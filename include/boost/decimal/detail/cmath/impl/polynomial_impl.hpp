// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_POLYNOMIAL_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_POLYNOMIAL_IMPL_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/remove_trailing_zeros.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace poly_detail {

using u128 = boost::int128::uint128_t;

// R_(k+1) = (a_k x + b_k) R_k - c_k R_(k-1) from R_0 = 1: hermite (a, b, c) = (2, 0, 2k) and P = R_n, legendre of order m
// (2k + 2m + 1, 0, k (k + 2m)) and P = (2m - 1)!! (1 - x^2)^(m/2) R_n / n!, laguerre (-1, 2k + m + 1, k (k + m)) and P = R_n / n!
enum class poly_kind { hermite, legendre, laguerre };

struct poly_step
{
    std::uint64_t a;
    std::uint64_t b;
    std::uint64_t c;
    bool a_neg;
};

constexpr auto poly_coeffs(const poly_kind kind, const std::uint64_t k, const std::uint64_t m) noexcept -> poly_step
{
    return kind == poly_kind::hermite ? poly_step {2U, 0U, 2U * k, false}
         : kind == poly_kind::legendre ? poly_step {2U * k + 2U * m + 1U, 0U, k * (k + 2U * m), false}
         : poly_step {1U, 2U * k + m + 1U, k * (k + m), true};
}

// m * 2^e with m in [1, 2), or zero
template <int N>
struct poly_float
{
    fx<N> m;
    int e;
    bool neg;
};

// As m is in [1, 2) or zero, its top word is zero only for zero
template <int N>
constexpr auto pf_is_zero(const poly_float<N>& a) noexcept -> bool
{
    return a.m.w[N - 1] == 0U;
}

template <int N>
constexpr auto pf_norm(const fx<N>& m, const int e, const bool neg) noexcept -> poly_float<N>
{
    const int l {fx_countl_zero(m)};
    return l == 64 * N ? poly_float<N> {m, 0, neg} : poly_float<N> {fx_shift(m, l - 1), e - (l - 1), neg};
}

// v as a float, normalized and truncated to N words
template <int N>
constexpr auto pf_from(const u128& v, const bool neg = false) noexcept -> poly_float<N>
{
    poly_float<N> r {fx<N> {}, 0, neg};
    if (v != 0U)
    {
        // t has its leading one at bit 126, and for z = 0 the last bit of v goes to the next word
        const int z {v.high != 0U ? countl_zero(v.high) : 64 + countl_zero(v.low)};
        const auto t {z == 0 ? v >> 1U : v << static_cast<unsigned>(z - 1)};
        r.m.w[N - 1] = t.high;
        BOOST_DECIMAL_IF_CONSTEXPR (N > 1)
        {
            r.m.w[N > 1 ? N - 2 : 0] = t.low;
        }
        BOOST_DECIMAL_IF_CONSTEXPR (N > 2)
        {
            r.m.w[N > 2 ? N - 3 : 0] = z == 0 ? v.low << 63U : 0U;
        }
        r.e = 127 - z;
    }
    return r;
}

// m * 2^e for m in [1, 4), with m brought into [1, 2)
template <int N>
constexpr auto pf_half(const fx<N>& m, const int e, const bool neg) noexcept -> poly_float<N>
{
    return m.w[N - 1] >= (UINT64_C(1) << 63U) ? poly_float<N> {fx_shr(m, 1U), e + 1, neg} : poly_float<N> {m, e, neg};
}

template <int N>
constexpr auto pf_mul(const poly_float<N>& a, const poly_float<N>& b) noexcept -> poly_float<N>
{
    return pf_is_zero(a) || pf_is_zero(b) ? poly_float<N> {} : pf_half(fx_mul(a.m, b.m), a.e + b.e, a.neg != b.neg);
}

// a * c for 0 < c < 2^63: a * c / 2^s is in [1/2, 2) for c in [2^(s - 1), 2^s)
template <int N>
constexpr auto pf_mul_u(const poly_float<N>& a, const std::uint64_t c) noexcept -> poly_float<N>
{
    if (pf_is_zero(a) || c == 0U)
    {
        return poly_float<N> {};
    }
    const auto s {static_cast<unsigned>(64 - countl_zero(c))};
    if (s == 1U)
    {
        return a;
    }
    const auto r {fx_mul_shr(a.m, c, s)};
    return r.w[N - 1] < (UINT64_C(1) << 62U) ? poly_float<N> {fx_shl(r, 1U), a.e + static_cast<int>(s) - 1, a.neg} : poly_float<N> {r, a.e + static_cast<int>(s), a.neg};
}

// a * k for k = first, first + step, ... up to last, with the factors gathered in one word while their product is
// below 2^63
template <int N>
constexpr auto pf_mul_range(poly_float<N> a, const std::uint64_t first, const std::uint64_t last, const std::uint64_t step) noexcept -> poly_float<N>
{
    std::uint64_t c {1U};
    for (std::uint64_t k {first}; k <= last; k += step)
    {
        if (128 - countl_zero(c) - countl_zero(k) > 63)
        {
            a = pf_mul_u(a, c);
            c = 1U;
        }
        c *= k;
    }
    return pf_mul_u(a, c);
}

// c as a float, for 0 < c < 2^63
template <int N>
constexpr auto pf_u(const std::uint64_t c) noexcept -> poly_float<N>
{
    poly_float<N> r {};
    const int b {64 - countl_zero(c)};
    r.m.w[N - 1] = c << (63 - b);
    r.e = b - 1;
    return r;
}

// a + b, or a - b for sub
template <int N>
constexpr auto pf_add(poly_float<N> a, poly_float<N> b, const bool sub) noexcept -> poly_float<N>
{
    b.neg = b.neg != sub;
    if (pf_is_zero(b))
    {
        return a;
    }
    if (pf_is_zero(a))
    {
        return b;
    }
    if (a.e < b.e)
    {
        const auto t {a};
        a = b;
        b = t;
    }
    const auto bm {fx_shift(b.m, b.e - a.e)};
    if (a.neg == b.neg)
    {
        return pf_half(fx_add(a.m, bm), a.e, a.neg);
    }
    return fx_less(a.m, bm) ? pf_norm(fx_sub(bm, a.m), a.e, b.neg) : pf_norm(fx_sub(a.m, bm), a.e, a.neg);
}

template <int N>
constexpr auto pf_recip(const poly_float<N>& a) noexcept -> poly_float<N>
{
    return pf_norm(fx_recip(a.m), -a.e, a.neg);
}

// r * b^n by squaring
template <int N>
constexpr auto pf_mul_pow(poly_float<N> r, poly_float<N> b, unsigned n) noexcept -> poly_float<N>
{
    while (n != 0U)
    {
        if ((n & 1U) != 0U)
        {
            r = pf_mul(r, b);
        }
        n >>= 1U;
        if (n != 0U)
        {
            b = pf_mul(b, b);
        }
    }
    return r;
}

// 10^k for |k| <= 38 from the table, with an error below one unit, and e = floor(k log2(10)) = floor(1701k / 512)
template <int N>
constexpr auto pf_pow10_small(const int k) noexcept -> std::enable_if_t<(N <= 3), poly_float<N>>
{
    return poly_float<N> {fx_top<N>(fx_table::pow10_bin[k + 38]), k >= 0 ? (k * 1701) >> 9 : -((511 - k * 1701) >> 9), false};
}

// The same in more words from 10^|k|, and its reciprocal for k < 0
template <int N>
constexpr auto pf_pow10_small(const int k) noexcept -> std::enable_if_t<(N > 3), poly_float<N>>
{
    const auto t {pf_from<N>(pow10(static_cast<u128>(k < 0 ? -k : k)))};
    return k < 0 ? pf_recip(t) : t;
}

// 10^k, by squaring for |k| > 38, with a relative error below 2 log2(|k|) + 2 units
template <int N>
constexpr auto pf_pow10(const int k) noexcept -> poly_float<N>
{
    const auto n {static_cast<unsigned>(k < 0 ? -k : k)};
    if (n <= 38U)
    {
        return pf_pow10_small<N>(k);
    }
    const auto r {pf_mul_pow(poly_float<N> {fx_one<N>(), 0, false}, pf_from<N>(u128 {10U}), n)};
    return k < 0 ? pf_recip(r) : r;
}

// (2m - 1)!! (1 - x^2)^(m/2) with a relative error below 16m + 32 units, from (10^q - sig) (10^q + sig) / 10^2q
// for x = sig * 10^-q where 10^q fits, or else 1 - x^2 with x below 10^-4
template <int N>
constexpr auto pf_legendre_scale(const unsigned m, const u128& sig, const int e, const poly_float<N>& x) noexcept -> poly_float<N>
{
    const poly_float<N> one {fx_one<N>(), 0, false};
    const int q {e < 0 ? -e : 0};
    const auto t {pow10(static_cast<u128>(q <= 38 ? q : 0))};
    const auto y {q <= 38 ? pf_mul(pf_mul(pf_from<N>(t - sig), pf_from<N>(t + sig)), pf_pow10<N>(-2 * q)) : pf_add(one, pf_mul(x, x), true)};
    // (2m - 1)!! in one word up to m = 16
    std::uint64_t w {1U};
    for (std::uint64_t k {3}; k < 2U * std::uint64_t {m} && k < 32U; k += 2U)
    {
        w *= k;
    }
    auto r {m <= 16U ? pf_u<N>(w) : pf_mul_range(one, 1U, 2U * std::uint64_t {m}, 2U)};
    r = pf_mul_pow(r, y, m / 2U);
    if ((m & 1U) != 0U)
    {
        const bool odd {(y.e & 1) != 0};
        r = pf_mul(r, pf_norm(fx_sqrt(odd ? fx_shl(y.m, 1U) : y.m), (y.e - (odd ? 1 : 0)) / 2, false));
    }
    return r;
}

// An upper bound m * 2^e
struct poly_bound
{
    std::uint64_t m;
    int e;
};

// v * 2^e rounded up to a mantissa of at most 2^31, so that products and sums fit in 64 bits
constexpr auto pb_make(const std::uint64_t v, const int e) noexcept -> poly_bound
{
    const int bits {64 - countl_zero(v)};
    if (bits <= 31)
    {
        return {v, e};
    }
    const int s {bits - 31};
    return {(v >> s) + ((v & ((UINT64_C(1) << s) - 1U)) != 0U ? 1U : 0U), e + s};
}

constexpr auto pb_mul(const poly_bound& a, const poly_bound& b) noexcept -> poly_bound
{
    return pb_make(a.m * b.m, a.e + b.e);
}

constexpr auto pb_add(poly_bound a, poly_bound b) noexcept -> poly_bound
{
    if (a.m == 0U || (b.m != 0U && a.e < b.e))
    {
        const auto t {a};
        a = b;
        b = t;
    }
    if (b.m == 0U)
    {
        return a;
    }
    const int d {a.e - b.e};
    return d >= 32 ? pb_make(a.m + 1U, a.e) : pb_make((a.m << d) + b.m, b.e);
}

// a1 b1 + a2 b2 rounded up, for products below 2^63
constexpr auto pb_dot(const poly_bound& a1, const poly_bound& b1, const poly_bound& a2, const poly_bound& b2) noexcept -> poly_bound
{
    auto p {a1.m * b1.m};
    auto e {a1.e + b1.e};
    auto q {a2.m * b2.m};
    auto f {a2.e + b2.e};
    if (p == 0U || (q != 0U && e < f))
    {
        const auto t {p};
        p = q;
        q = t;
        const int g {e};
        e = f;
        f = g;
    }
    if (q == 0U)
    {
        return pb_make(p, e);
    }
    return pb_make(p + (e - f >= 63 ? 1U : (q >> (e - f)) + 1U), e);
}

// A coefficient as a bound, kept whole below 2^32 for pb_dot
constexpr auto pb_coef(const std::uint64_t v) noexcept -> poly_bound
{
    return v < (UINT64_C(1) << 32U) ? poly_bound {v, 0} : pb_make(v, 0);
}

template <int N>
constexpr auto pb_of(const poly_float<N>& a) noexcept -> poly_bound
{
    return pf_is_zero(a) ? poly_bound {0U, 0} : poly_bound {(a.m.w[N - 1] >> 32U) + 1U, a.e - 30};
}

// A bound as an integer, at most 2^126
constexpr auto pb_int(const poly_bound& a) noexcept -> u128
{
    if (a.e >= 0)
    {
        return a.e > 94 ? u128 {1U} << 126U : static_cast<u128>(a.m) << a.e;
    }
    return a.e <= -32 ? u128 {1U} : static_cast<u128>((a.m >> -a.e) + 1U);
}

// a, or -a for neg, modulo 2^(64N)
template <int N>
constexpr auto poly_cneg(fx<N> a, const bool neg) noexcept -> fx<N>
{
    const std::uint64_t mask {neg ? UINT64_MAX : 0U};
    std::uint64_t c {neg ? 1U : 0U};
    BOOST_DECIMAL_FX_UNROLL
    for (int i {}; i < N; ++i)
    {
        const std::uint64_t v {(a.w[i] ^ mask) + c};
        c = v < c ? 1U : 0U;
        a.w[i] = v;
    }
    return a;
}

// a * n * 2^-s for 64 <= s, truncated, out of line with its arguments by value so that they stay in registers
template <int N>
constexpr auto poly_mul_shift_far(const fx<N> a, const std::uint64_t n, const int s) noexcept -> fx<N>
{
    fx<N + 1> p {};
    p.w[N] = fx_mul_add(fx_make_span<N>(p.w), fx_make_span(a.w), n);
    p = fx_shift(p, -s);
    fx<N> r {};
    for (int i {}; i < N; ++i)
    {
        r.w[i] = p.w[i];
    }
    return r;
}

// a * n * 2^-s for a in one word, truncated, for s >= 0 and a result below 2
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_mul_shift1(const std::uint64_t a, const std::uint64_t n, const int s) noexcept -> std::uint64_t
{
    std::uint64_t hi {};
    const auto lo {boost::int128::detail::umul(a, n, hi)};
    if (s >= 64)
    {
        return s >= 128 ? 0U : hi >> static_cast<unsigned>(s - 64);
    }
    const auto b {static_cast<unsigned>(s)};
    return (lo >> b) | ((hi << 1U) << (63U - b));
}

// 1/n! for n <= 20, normalized from the table, where its leading one is in the top word
template <int N>
constexpr auto pf_inv_factorial(const unsigned n) noexcept -> poly_float<N>
{
    const auto& t {fx_table::inv_factorial[n]};
    const int z {countl_zero(t.w[2]) - 1};
    return {fx_top<N>(fx_extend<(N > 3 ? N : 3)>(z == 0 ? t : fx_shl(t, static_cast<unsigned>(z)))), -z, false};
}

// The factor f that turns R_n into P: 1 / n! and, for legendre of order m, times (2m - 1)!! (1 - x^2)^(m/2), with
// a relative error below rho u
template <int N>
struct poly_factor
{
    poly_float<N> f;
    std::uint64_t rho;
};

// n! in one word, for n <= 20
constexpr auto poly_factorial(const unsigned n) noexcept -> std::uint64_t
{
    std::uint64_t w {1U};
    for (std::uint64_t k {2}; k <= n; ++k)
    {
        w *= k;
    }
    return w;
}

template <int N, bool Assoc>
constexpr auto poly_factor_of(const unsigned n, const unsigned m, const u128& sig, const int e, const poly_float<N>& x) noexcept -> poly_factor<N>
{
    // In one or two words 1/n! comes from the table, which keeps at least 128 bits of it up to n = 20
    poly_factor<N> r {N <= 2 && n <= 20U ? pf_inv_factorial<N>(n)
                                         : pf_recip(n <= 20U ? pf_u<N>(poly_factorial(n)) : pf_mul_range(poly_float<N> {fx_one<N>(), 0, false}, 2U, n, 1U)),
                      std::uint64_t {n} + 4U};
    BOOST_DECIMAL_IF_CONSTEXPR (Assoc)
    {
        r.f = pf_mul(r.f, pf_legendre_scale(m, sig, e, x));
        r.rho += 16U * std::uint64_t {m} + 32U;
    }
    return r;
}

// The digits of a result: d * 10^e10 with 38 or 39 digits and an error below units. amb is set where the error can reach
// a rounding boundary, and bound is the boundary next to d.
struct poly_out
{
    u128 d;
    u128 bound;
    u128 units;
    int e10;
    bool neg;
    bool amb;
};

// P(x) for x = sig * 10^e in N words with a bound on the error of each step: the error of R_(k+1) is below |a x + b|
// times that of R_k, plus c times that of R_(k-1), plus the rounding of the step
template <int N, bool Assoc>
constexpr auto poly_eval(const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg, const int p) noexcept -> poly_out
{
    const poly_bound u {1U, -(64 * N - 8)};
    const auto x {pf_mul(pf_from<N>(sig, xneg), pf_pow10<N>(e))};
    const auto bx {pb_of(x)};
    fx<N> f0 {};
    auto f1 {fx_shr(fx_one<N>(), 1U)};
    int d0 {};
    int d1 {1};
    poly_bound bnd0 {0U, 0};
    poly_bound bnd1 {1U, 0};
    // R_k = S_k 2^d_k with 2^d_k <= 8 B_k, so with big = |a x| + |b| and a margin for x, the rounding of a step is below
    // u (44 big B_k + 3 c B_(k-1)), and the error of R_k below 47 k u B_k for B_(k+1) = big B_k + c B_(k-1), which bounds |R_k|
    const auto bxm {pb_add(bx, pb_mul(bx, poly_bound {1U, -40}))};
    for (std::uint64_t k {}; k < n; ++k)
    {
        const auto s {poly_coeffs(kind, k, m)};
        auto coef {pf_mul_u(x, s.a)};
        coef.neg = coef.neg != s.a_neg;
        if (s.b != 0U)
        {
            coef = pf_add(coef, pf_u<N>(s.b), false);
        }
        const auto big {pb_dot(pb_coef(s.a), bxm, pb_coef(s.b), poly_bound {1U, 0})};
        const auto bnd2 {pb_dot(big, bnd1, pb_coef(s.c), bnd0)};
        const bool n1 {(f1.w[N - 1] >> 63U) != 0U};
        const bool n0 {(f0.w[N - 1] >> 63U) != 0U};
        const auto a1 {poly_cneg(f1, n1)};
        const auto a0 {poly_cneg(f0, n0)};
        const auto sc {static_cast<unsigned>(65 - countl_zero(s.c))};
        // Each term is below 2^(d2 - 1), and a zero takes no part
        const int l1 {fx_countl_zero(a1)};
        const int l0 {fx_countl_zero(a0)};
        const int e1 {l1 == 64 * N ? (std::numeric_limits<int>::min)() / 4 : coef.e + d1 + 3 - l1};
        const int e0 {l0 == 64 * N ? (std::numeric_limits<int>::min)() / 4 : static_cast<int>(sc) + d0 + 1 - l0};
        const int d2 {(e1 > e0 ? e1 : e0) + 1};
        const auto t1 {fx_shift(fx_mul(coef.m, a1), d1 + coef.e - d2)};
        const auto t2 {fx_shift(fx_mul_shr(a0, s.c, sc), d0 + static_cast<int>(sc) - d2)};
        const auto f2 {fx_add(poly_cneg(t1, n1 != coef.neg), poly_cneg(t2, !n0))};
        f0 = f1;
        f1 = f2;
        d0 = d1;
        d1 = d2;
        bnd0 = bnd1;
        bnd1 = bnd2;
    }
    const bool fneg {(f1.w[N - 1] >> 63U) != 0U};
    auto value {pf_norm(poly_cneg(f1, fneg), d1, fneg)};
    auto err {pb_mul(pb_mul(u, pb_make(47U * std::uint64_t {n}, 0)), bnd1)};
    if (kind != poly_kind::hermite)
    {
        const auto q {poly_factor_of<N, Assoc>(n, m, sig, e, x)};
        value = pf_mul(value, q.f);
        err = pb_add(pb_mul(pb_mul(err, pb_of(q.f)), pb_make((UINT64_C(1) << 40U) + 1U, -40)), pb_mul(pb_mul(u, pb_make(q.rho, 0)), pb_of(value)));
    }

    // A zero value leaves |P| at most err < 10^(t + 1): the units saturate, and the exact integers decide from there
    poly_out out {0U, 0U, 0U, 0, value.neg, true};
    if (pf_is_zero(value))
    {
        out.e10 = fx_floor_log10_2(err.e + 64 - countl_zero(err.m)) - 35;
        out.units = u128 {1U} << 126U;
        return out;
    }

    // |value| in [2^E, 2^(E + 1)) with 10^t <= 2^E < 10^(t + 1), or t one less: |value| * 10^(37 - t) is in
    // [10^37, 2.1 * 10^38)
    const int t {fx_floor_log10_2(value.e)};
    const int k10 {37 - t};
    const auto ten {pf_pow10<N>(k10)};
    const auto v {pf_mul(value, ten)};
    u128 hi {};
    hi.high = v.m.w[N - 1];
    hi.low = v.m.w[N - 2];
    // v = hi * 2^(v.e - 126), with the words below hi left out
    out.d = v.e >= 126 ? hi << (v.e - 126) : hi >> (126 - v.e);
    out.e10 = -k10;

    // err * 10^k10 + |v| (2 log2 |k10| + 8 <= 40) units, and 4 for the words left out
    auto total {pb_add(pb_mul(err, pb_of(ten)), pb_mul(pb_mul(u, pb_make(40U, 0)), pb_of(v)))};
    const auto units {pb_int(total) + 4U};
    // units saturates at 2^126, which is at least h
    const auto h {pow10(static_cast<u128>(37 - p)) * 5U};
    const auto r {out.d % h};
    out.units = units;
    out.amb = units >= h || r <= units || r + units + 1U >= h;
    out.bound = r <= units ? out.d - r : out.d - r + h;
    return out;
}

// An upper bound of 16 log2(v) for v >= 1: 16 e + t + 3 for v in [2^e, 2^(e + 1)) with the next four bits t, from the
// tangent of log2(1 + u) at u = 1/2
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_lg16(const std::uint64_t v) noexcept -> int
{
    const int z {countl_zero(v)};
    return 16 * (63 - z) + static_cast<int>(((v << static_cast<unsigned>(z)) >> 59U) & 15U) + 3;
}

// An upper bound of 16 log2(2^(a/16) + 2^(b/16)): 16 log2(1 + 2^-y) <= 23 2^-y, and 2^-r <= 1 - r/2 for r in [0, 1].
// 23 * 32 < 2^10, thus a shift by 10 already gives 0.
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_add16(const int a, const int b) noexcept -> int
{
    const int d {a > b ? a - b : b - a};
    const int q {(d >> 4) < 5 ? d >> 4 : 5};
    const int t {((23 * (32 - (d & 15))) >> (5 + q)) + 1};
    return (a > b ? a : b) + (t < 16 ? t : 16);
}

// An upper bound of 16 log2 |a| as poly_lg16 for a nonzero float
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto pf_lg16(const poly_float<N>& a) noexcept -> int
{
    return 16 * a.e + static_cast<int>((a.m.w[N - 1] >> 58U) & 15U) + 3;
}

// The 16 log2 bound of zero
constexpr int poly_none {(std::numeric_limits<int>::min)() / 8};

// g_(k+1) as an upper bound of 16 log2 B_(k+1) from g_k and g_(k-1) for |x| < 2^(lx / 16), and d_(k+1) = ceil(g_(k+1) / 16) + 1
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_grow16(const poly_step& s, const int lx, const bool xz, int& g0, int& g1) noexcept -> int
{
    const int la {xz ? poly_none : poly_lg16(s.a) + lx};
    const int l {s.b == 0U ? la : la == poly_none ? poly_lg16(s.b) : poly_add16(la, poly_lg16(s.b))};
    const int g2 {s.c == 0U ? l + g1 : poly_add16(l + g1, poly_lg16(s.c) + g0)};
    g0 = g1;
    g1 = g2;
    return (g2 >= 0 ? (g2 + 15) / 16 : -(-g2 / 16)) + 1;
}

// The loop of poly_try: R_n = f 2^d with B_n <= bnd, or B_n < 2^(d - 1) in two words or more
template <int N>
struct poly_run
{
    fx<N> f;
    poly_bound bnd;
    int d;
};

// The loop of poly_try in one word for each kind, with plain words
template <poly_kind kind>
constexpr auto poly_loop1(const unsigned n, const unsigned m, const poly_float<1>& x) noexcept -> poly_run<1>
{
    // |x| < bxm also covers a margin of 2^-40
    const auto bx {pb_of(x)};
    const auto bxm {pb_add(bx, pb_mul(bx, poly_bound {1U, -40}))};
    const bool xz {pf_is_zero(x)};
    const std::uint64_t xm {x.m.w[0]};
    poly_bound bnd0 {0U, 0};
    poly_bound bnd1 {1U, 0};
    std::uint64_t f0 {};
    std::uint64_t f1 {fx_shr(fx_one<1>(), 1U).w[0]};
    int d0 {};
    int d1 {1};
    for (std::uint64_t k {}; k < n; ++k)
    {
        const auto s {poly_coeffs(kind, k, m)};
        const auto big {pb_dot(pb_coef(s.a), bxm, pb_coef(s.b), poly_bound {1U, 0})};
        const auto bnd2 {pb_dot(big, bnd1, pb_coef(s.c), bnd0)};
        const int d2 {bnd2.m == 0U ? d1 : bnd2.e + 65 - countl_zero(bnd2.m)};
        bnd0 = bnd1;
        bnd1 = bnd2;
        const bool n1 {(f1 >> 63U) != 0U};
        const bool n0 {(f0 >> 63U) != 0U};
        const std::uint64_t a1 {n1 ? 0U - f1 : f1};
        const std::uint64_t a0 {n0 ? 0U - f0 : f0};
        // a x R_k, b R_k and c R_(k-1) in units of 2^d2, each with the factor 2 in n and 1 more in the shift
        std::uint64_t t {};
        if (!xz)
        {
            std::uint64_t hi {};
            const auto lo {boost::int128::detail::umul(xm, a1, hi)};
            t = poly_mul_shift1((lo >> 62U) + (hi << 2U), 2U * s.a, d2 - d1 - x.e + 1);
            t = n1 != (x.neg != s.a_neg) ? 0U - t : t;
        }
        if (s.b != 0U)
        {
            const auto v {poly_mul_shift1(a1, 2U * s.b, d2 - d1 + 1)};
            t += n1 ? 0U - v : v;
        }
        if (s.c != 0U)
        {
            const auto v {poly_mul_shift1(a0, 2U * s.c, d2 - d0 + 1)};
            t += n0 ? v : 0U - v;
        }
        f0 = f1;
        f1 = t;
        d0 = d1;
        d1 = d2;
    }
    fx<1> f {};
    f.w[0] = f1;
    return {f, bnd1, d1};
}

constexpr auto poly_loop(const poly_kind kind, const unsigned n, const unsigned m, const poly_float<1>& x) noexcept -> poly_run<1>
{
    return kind == poly_kind::hermite ? poly_loop1<poly_kind::hermite>(n, m, x)
         : kind == poly_kind::legendre ? poly_loop1<poly_kind::legendre>(n, m, x)
         : poly_loop1<poly_kind::laguerre>(n, m, x);
}

// p = a * n for a in two words, and p * 2^-s truncated to two words, for s >= 0
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_mul_shift2(const std::uint64_t lo, const std::uint64_t hi, const std::uint64_t n, const int s, std::uint64_t& rl, std::uint64_t& rh) noexcept -> void
{
    if (s >= 64)
    {
        fx<2> a {};
        a.w[0] = lo;
        a.w[1] = hi;
        const auto r {poly_mul_shift_far<2>(a, n, s)};
        rl = r.w[0];
        rh = r.w[1];
        return;
    }
    std::uint64_t c {};
    const auto p0 {mul_add(lo, n, 0U, c)};
    const auto p1 {mul_add(hi, n, 0U, c)};
    const auto b {static_cast<unsigned>(s)};
    rl = (p0 >> b) | ((p1 << 1U) << (63U - b));
    rh = (p1 >> b) | ((c << 1U) << (63U - b));
}

// lo, hi negated modulo 2^128 for neg
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_cneg2(std::uint64_t& lo, std::uint64_t& hi, const bool neg) noexcept -> void
{
    const std::uint64_t mask {neg ? UINT64_MAX : 0U};
    lo = (lo ^ mask) + (neg ? 1U : 0U);
    hi = (hi ^ mask) + (neg && lo == 0U ? 1U : 0U);
}

// lo, hi plus vl, vh modulo 2^128
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_add2(std::uint64_t& lo, std::uint64_t& hi, const std::uint64_t vl, const std::uint64_t vh) noexcept -> void
{
    lo += vl;
    hi += vh + (lo < vl ? 1U : 0U);
}

// The loop above in two words for each kind, with plain words that GCC keeps in registers and not the arrays of fx<2>
template <poly_kind kind>
constexpr auto poly_loop2(const unsigned n, const unsigned m, const poly_float<2>& x) noexcept -> poly_run<2>
{
    // In two words or more, B_k = 2^(g_k / 16) with g_k from upper bounds of 16 log2 is enough, and 2^(d_k - 1) >= B_k.
    // |x| < 2^(lx / 16) also covers a margin of 2^-40.
    const int lx {pf_lg16(x)};
    const bool xz {pf_is_zero(x)};
    const std::uint64_t xl {x.m.w[0]};
    const std::uint64_t xh {x.m.w[1]};
    int g0 {poly_none};
    int g1 {};
    std::uint64_t f0l {};
    std::uint64_t f0h {};
    std::uint64_t f1l {};
    std::uint64_t f1h {fx_shr(fx_one<2>(), 1U).w[1]};
    int d0 {};
    int d1 {1};
    for (std::uint64_t k {}; k < n; ++k)
    {
        const auto s {poly_coeffs(kind, k, m)};
        const int d2 {poly_grow16(s, lx, xz, g0, g1)};
        const bool n1 {(f1h >> 63U) != 0U};
        const bool n0 {(f0h >> 63U) != 0U};
        std::uint64_t a1l {f1l};
        std::uint64_t a1h {f1h};
        std::uint64_t a0l {f0l};
        std::uint64_t a0h {f0h};
        poly_cneg2(a1l, a1h, n1);
        poly_cneg2(a0l, a0h, n0);
        std::uint64_t tl {};
        std::uint64_t th {};
        if (!xz)
        {
            // x * |R_k| as in fx_mul_body
            std::uint64_t c {};
            const auto p1 {mul_add(xl, a1h, 0U, c)};
            std::uint64_t p2 {c};
            c = 0U;
            const auto q1 {mul_add(xh, a1l, p1, c)};
            p2 = mul_add(xh, a1h, p2, c);
            poly_mul_shift2((q1 >> 62U) | (p2 << 2U), (p2 >> 62U) | (c << 2U), 2U * s.a, d2 - d1 - x.e + 1, tl, th);
            poly_cneg2(tl, th, n1 != (x.neg != s.a_neg));
        }
        if (s.b != 0U)
        {
            std::uint64_t vl {};
            std::uint64_t vh {};
            poly_mul_shift2(a1l, a1h, 2U * s.b, d2 - d1 + 1, vl, vh);
            poly_cneg2(vl, vh, n1);
            poly_add2(tl, th, vl, vh);
        }
        if (s.c != 0U)
        {
            std::uint64_t vl {};
            std::uint64_t vh {};
            poly_mul_shift2(a0l, a0h, 2U * s.c, d2 - d0 + 1, vl, vh);
            poly_cneg2(vl, vh, !n0);
            poly_add2(tl, th, vl, vh);
        }
        f0l = f1l;
        f0h = f1h;
        f1l = tl;
        f1h = th;
        d0 = d1;
        d1 = d2;
    }
    fx<2> f {};
    f.w[0] = f1l;
    f.w[1] = f1h;
    return {f, poly_bound {1U, 0}, d1};
}

// Two words take poly_loop2 as one word takes poly_loop1
constexpr auto poly_loop(const poly_kind kind, const unsigned n, const unsigned m, const poly_float<2>& x) noexcept -> poly_run<2>
{
    return kind == poly_kind::hermite ? poly_loop2<poly_kind::hermite>(n, m, x)
         : kind == poly_kind::legendre ? poly_loop2<poly_kind::legendre>(n, m, x)
         : poly_loop2<poly_kind::laguerre>(n, m, x);
}

// p = a * n for a in three words, and p * 2^-s truncated to three words, for s >= 0
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_mul_shift3(const std::uint64_t a0, const std::uint64_t a1, const std::uint64_t a2, const std::uint64_t n, const int s,
                                                          std::uint64_t& r0, std::uint64_t& r1, std::uint64_t& r2) noexcept -> void
{
    if (s >= 64)
    {
        fx<3> a {};
        a.w[0] = a0;
        a.w[1] = a1;
        a.w[2] = a2;
        const auto r {poly_mul_shift_far<3>(a, n, s)};
        r0 = r.w[0];
        r1 = r.w[1];
        r2 = r.w[2];
        return;
    }
    std::uint64_t c {};
    const auto p0 {mul_add(a0, n, 0U, c)};
    const auto p1 {mul_add(a1, n, 0U, c)};
    const auto p2 {mul_add(a2, n, 0U, c)};
    const auto b {static_cast<unsigned>(s)};
    r0 = (p0 >> b) | ((p1 << 1U) << (63U - b));
    r1 = (p1 >> b) | ((p2 << 1U) << (63U - b));
    r2 = (p2 >> b) | ((c << 1U) << (63U - b));
}

// a0, a1, a2 negated modulo 2^192 for neg
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_cneg3(std::uint64_t& a0, std::uint64_t& a1, std::uint64_t& a2, const bool neg) noexcept -> void
{
    const std::uint64_t mask {neg ? UINT64_MAX : 0U};
    a0 = (a0 ^ mask) + (neg ? 1U : 0U);
    const std::uint64_t c1 {neg && a0 == 0U ? 1U : 0U};
    a1 = (a1 ^ mask) + c1;
    a2 = (a2 ^ mask) + (c1 != 0U && a1 == 0U ? 1U : 0U);
}

// a0, a1, a2 plus v0, v1, v2 modulo 2^192
BOOST_DECIMAL_FORCE_INLINE constexpr auto poly_add3(std::uint64_t& a0, std::uint64_t& a1, std::uint64_t& a2, const std::uint64_t v0, const std::uint64_t v1, const std::uint64_t v2) noexcept -> void
{
    const u128 v {v1, v0};
    const auto r {u128 {a1, a0} + v};
    a2 += v2 + (r < v ? 1U : 0U);
    a1 = r.high;
    a0 = r.low;
}

// The loop in three words for each kind, as poly_loop2
template <poly_kind kind>
constexpr auto poly_loop3(const unsigned n, const unsigned m, const poly_float<3>& x) noexcept -> poly_run<3>
{
    const int lx {pf_lg16(x)};
    const bool xz {pf_is_zero(x)};
    const std::uint64_t x0 {x.m.w[0]};
    const std::uint64_t x1 {x.m.w[1]};
    const std::uint64_t x2 {x.m.w[2]};
    int g0 {poly_none};
    int g1 {};
    std::uint64_t f00 {};
    std::uint64_t f01 {};
    std::uint64_t f02 {};
    std::uint64_t f10 {};
    std::uint64_t f11 {};
    std::uint64_t f12 {fx_shr(fx_one<3>(), 1U).w[2]};
    int d0 {};
    int d1 {1};
    for (std::uint64_t k {}; k < n; ++k)
    {
        const auto s {poly_coeffs(kind, k, m)};
        const int d2 {poly_grow16(s, lx, xz, g0, g1)};
        const bool n1 {(f12 >> 63U) != 0U};
        const bool n0 {(f02 >> 63U) != 0U};
        std::uint64_t a10 {f10};
        std::uint64_t a11 {f11};
        std::uint64_t a12 {f12};
        std::uint64_t a00 {f00};
        std::uint64_t a01 {f01};
        std::uint64_t a02 {f02};
        poly_cneg3(a10, a11, a12, n1);
        poly_cneg3(a00, a01, a02, n0);
        std::uint64_t t0 {};
        std::uint64_t t1 {};
        std::uint64_t t2 {};
        if (!xz)
        {
            // x * |R_k| as in fx_mul_body
            std::uint64_t c {};
            std::uint64_t p2 {mul_add(x0, a12, 0U, c)};
            std::uint64_t p3 {c};
            c = 0U;
            p2 = mul_add(x1, a11, p2, c);
            p3 = mul_add(x1, a12, p3, c);
            std::uint64_t p4 {c};
            c = 0U;
            p2 = mul_add(x2, a10, p2, c);
            p3 = mul_add(x2, a11, p3, c);
            p4 = mul_add(x2, a12, p4, c);
            poly_mul_shift3((p2 >> 62U) | (p3 << 2U), (p3 >> 62U) | (p4 << 2U), (p4 >> 62U) | (c << 2U), 2U * s.a, d2 - d1 - x.e + 1, t0, t1, t2);
            poly_cneg3(t0, t1, t2, n1 != (x.neg != s.a_neg));
        }
        if (s.b != 0U)
        {
            std::uint64_t v0 {};
            std::uint64_t v1 {};
            std::uint64_t v2 {};
            poly_mul_shift3(a10, a11, a12, 2U * s.b, d2 - d1 + 1, v0, v1, v2);
            poly_cneg3(v0, v1, v2, n1);
            poly_add3(t0, t1, t2, v0, v1, v2);
        }
        if (s.c != 0U)
        {
            std::uint64_t v0 {};
            std::uint64_t v1 {};
            std::uint64_t v2 {};
            poly_mul_shift3(a00, a01, a02, 2U * s.c, d2 - d0 + 1, v0, v1, v2);
            poly_cneg3(v0, v1, v2, !n0);
            poly_add3(t0, t1, t2, v0, v1, v2);
        }
        f00 = f10;
        f01 = f11;
        f02 = f12;
        f10 = t0;
        f11 = t1;
        f12 = t2;
        d0 = d1;
        d1 = d2;
    }
    fx<3> f {};
    f.w[0] = f10;
    f.w[1] = f11;
    f.w[2] = f12;
    return {f, poly_bound {1U, 0}, d1};
}

// Three words take poly_loop3 as two words take poly_loop2
constexpr auto poly_loop(const poly_kind kind, const unsigned n, const unsigned m, const poly_float<3>& x) noexcept -> poly_run<3>
{
    return kind == poly_kind::hermite ? poly_loop3<poly_kind::hermite>(n, m, x)
         : kind == poly_kind::legendre ? poly_loop3<poly_kind::legendre>(n, m, x)
         : poly_loop3<poly_kind::laguerre>(n, m, x);
}

// The digits d * 10^e10 of poly_try, to round once with fx_round_kept where ok
struct poly_try_out
{
    u128 d;
    int e10;
    bool neg;
    bool ok;
};

// P(x) as above with the scale of each step from a bound B_k of |R_k| and not from the leading bits of the values. Where
// B_k is far above |R_k|, the value keeps few bits.
template <int N, bool Assoc>
constexpr auto poly_try(const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg, const int p, const int lo10, const int hi10) noexcept -> poly_try_out
{
    const poly_bound u {1U, -(64 * N - 8)};
    const auto x {pf_mul(pf_from<N>(sig, xneg), pf_pow10<N>(e))};
    const auto r {poly_loop(kind, n, m, x)};
    const auto& f1 {r.f};
    const int d1 {r.d};
    const bool fneg {(f1.w[N - 1] >> 63U) != 0U};
    auto value {pf_norm(poly_cneg(f1, fneg), d1, fneg)};
    // The error is below 2^(l / 16) with l from upper bounds of 16 log2, and 47 n u 2^(d1 - 1) for R_n
    int l {(N == 1 ? poly_lg16(r.bnd.m) + 16 * r.bnd.e : 16 * (d1 - 1)) + poly_lg16(47U * std::uint64_t {n > 0U ? n : 1U}) + 16 * u.e};
    if (kind != poly_kind::hermite)
    {
        // err |f| (1 + 2^-30) + rho u |value|
        const auto q {poly_factor_of<N, Assoc>(n, m, sig, e, x)};
        value = pf_mul(value, q.f);
        l = poly_add16(l + pf_lg16(q.f) + 1, poly_lg16(q.rho) + 16 * u.e + pf_lg16(value));
    }
    poly_try_out out {0U, 0, value.neg, false};
    if (!pf_is_zero(value))
    {
        // |value| 10^(p - t) = D + F / 2^64 with D of p + 1 or p + 2 digits, to round once in lo10 <= e10 <= hi10 where
        // the error interval holds no integer: the digits past D are then not all zero
        const int k10 {p - fx_floor_log10_2(value.e)};
        const auto ten {pf_pow10<N>(k10)};
        const auto v {pf_mul(value, ten)};
        // err |ten| + 40 u |v| below 2^(lt / 16) <= 2^-3
        const int lt {poly_add16(l + pf_lg16(ten), poly_lg16(40U) + 16 * u.e + pf_lg16(v))};
        const int shift {64 * N - 66 - v.e};
        if ((N == 1 || shift >= 0) && lt <= -48)
        {
            // In one word D and F share the word, as D has few digits
            fx<N> w {};
            BOOST_DECIMAL_IF_CONSTEXPR (N == 1)
            {
                w.w[0] = v.m.w[0] << static_cast<unsigned>(v.e + 2);
            }
            else
            {
                w = fx_shift(v.m, -shift);
            }
            // 2^(lt / 16 + 64) <= 2^(h + 64) (1 + r / 16) for lt = 16 h + r, as 2^y <= 1 + y for y in [0, 1]
            const int h {-((15 - lt) / 16)};
            const std::uint64_t eps {h >= -60 ? (static_cast<std::uint64_t>(16 + lt - 16 * h) << static_cast<unsigned>(h + 60)) + 3U : 35U};
            u128 dd {N == 1 ? v.m.w[0] >> static_cast<unsigned>(62 - v.e) : w.w[N > 1 ? 1 : 0]};
            BOOST_DECIMAL_IF_CONSTEXPR (N > 2)
            {
                dd.high = w.w[2];
            }
            out.e10 = -k10;
            if (dd >= pow10(static_cast<u128>(p + 1)))
            {
                // One word where it fits, as a 128-bit division takes a loop
                dd = dd.high == 0U ? u128 {dd.low / 10U} : dd / 10U;
                ++out.e10;
            }
            if (((w.w[0] > eps && w.w[0] < ~eps) || out.e10 > hi10) && out.e10 + 1 >= lo10)
            {
                out.d = dd;
                out.ok = true;
            }
        }
    }
    return out;
}

// A signed integer of up to poly_big_words words, low word first, where ok turns false once a result does not fit
constexpr int poly_big_words {576};

struct poly_big
{
    std::uint64_t w[poly_big_words];
    int len;
    bool neg;
    bool ok;
};

constexpr auto big_trim(poly_big& r) noexcept -> void
{
    while (r.len > 0 && r.w[r.len - 1] == 0U)
    {
        --r.len;
    }
}

// a in r, by the words that a uses
constexpr auto big_copy(const poly_big& a, poly_big& r) noexcept -> void
{
    for (int i {}; i < a.len; ++i)
    {
        r.w[i] = a.w[i];
    }
    r.len = a.len;
    r.neg = a.neg;
    r.ok = a.ok;
}

constexpr auto big_mul_u(poly_big& a, const std::uint64_t f) noexcept -> void
{
    std::uint64_t carry {};
    for (int i {}; i < a.len; ++i)
    {
        const auto t {static_cast<u128>(a.w[i]) * f + carry};
        a.w[i] = t.low;
        carry = t.high;
    }
    a.ok = a.ok && (carry == 0U || a.len < poly_big_words);
    if (carry != 0U && a.ok)
    {
        a.w[a.len++] = carry;
    }
}

constexpr auto big_mul_pow10(poly_big& a, int j) noexcept -> void
{
    for (; j > 0 && a.ok; j -= 19)
    {
        big_mul_u(a, pow10(static_cast<std::uint64_t>(j < 19 ? j : 19)));
    }
}

// v * 10^j in r
constexpr auto big_set(poly_big& r, const u128& v, const int j) noexcept -> void
{
    r.w[0] = v.low;
    r.w[1] = v.high;
    r.len = 2;
    r.neg = false;
    r.ok = true;
    big_trim(r);
    big_mul_pow10(r, j);
}

// a * b in r, which is neither
constexpr auto big_mul(const poly_big& a, const poly_big& b, poly_big& r) noexcept -> void
{
    r.ok = a.ok && b.ok && a.len + b.len <= poly_big_words;
    r.len = r.ok ? a.len + b.len : 0;
    r.neg = a.neg != b.neg;
    for (int i {}; i < r.len; ++i)
    {
        r.w[i] = 0U;
    }
    for (int i {}; i < (r.ok ? a.len : 0); ++i)
    {
        std::uint64_t carry {};
        for (int j {}; j < b.len; ++j)
        {
            const auto t {static_cast<u128>(a.w[i]) * b.w[j] + r.w[i + j] + carry};
            r.w[i + j] = t.low;
            carry = t.high;
        }
        r.w[i + b.len] = carry;
    }
    big_trim(r);
}

constexpr auto big_cmp(const poly_big& a, const poly_big& b) noexcept -> int
{
    if (a.len != b.len)
    {
        return a.len < b.len ? -1 : 1;
    }
    for (int i {a.len - 1}; i >= 0; --i)
    {
        if (a.w[i] != b.w[i])
        {
            return a.w[i] < b.w[i] ? -1 : 1;
        }
    }
    return 0;
}

// a + b, or a - b for sub, in r, which may be a
constexpr auto big_add(const poly_big& a, const poly_big& b, const bool sub, poly_big& r) noexcept -> void
{
    const bool bneg {b.neg != sub};
    const bool fits {a.len < poly_big_words && b.len < poly_big_words};
    const int len {!fits ? 0 : a.len > b.len ? a.len : b.len};
    const bool ok {a.ok && b.ok && fits};
    const bool swap {a.neg != bneg && big_cmp(a, b) < 0};
    const auto& x {swap ? b : a};
    const auto& y {swap ? a : b};
    const bool neg {swap ? bneg : a.neg};
    std::uint64_t carry {};
    for (int i {}; i < len; ++i)
    {
        const std::uint64_t xi {i < x.len ? x.w[i] : 0U};
        const std::uint64_t yi {i < y.len ? y.w[i] : 0U};
        if (a.neg == bneg)
        {
            const auto t {static_cast<u128>(xi) + yi + carry};
            r.w[i] = t.low;
            carry = t.high;
        }
        else
        {
            r.w[i] = xi - yi - carry;
            carry = (xi < yi || (xi == yi && carry != 0U)) ? 1U : 0U;
        }
    }
    r.w[len] = a.neg == bneg ? carry : 0U;
    r.len = len + 1;
    big_trim(r);
    r.neg = neg;
    r.ok = ok;
}

// The results of an exact comparison besides the sign -1, 0 or 1
constexpr int poly_unknown {2};
constexpr int poly_is_zero {3};

// P(x) in exact integers G_k = R_k 10^(qk) for x = sig * 10^e and q = max(-e, 0). a is |G_n|, or for legendre of order m
// |L| with L = G_n (2m - 1)!! M^(m/2) and M = 10^2q - s^2, or L^2 for an odd m (square). sign is that of P or poly_is_zero.
struct poly_exact
{
    poly_big g[3];
    poly_big t[2];
    poly_big* a;
    poly_big* w;
    int qn;
    int sign;
    bool square;
};

// G_0 = f * fu in place of 1 multiplies P by f * fu
template <bool Assoc>
constexpr auto poly_exact_init(poly_exact& st,const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg,
                               const u128& f = 1U, const std::uint64_t fu = 1U) noexcept -> void
{
    const int q {e < 0 ? -e : 0};
    auto& g {st.g};
    auto& t {st.t};
    big_set(g[0], 0U, 0);
    big_set(g[1], f, 0);
    big_mul_u(g[1], fu);
    // G_k in g[(k + 1) % 3], and t for the terms of a step
    unsigned steps {};
    for (unsigned k {}; k < n && g[(k + 1) % 3].ok; steps = ++k)
    {
        // G_(k+1) = (a_k s + b_k 10^q) G_k - c_k 10^(2q) G_(k-1) for s = sig * 10^max(e, 0)
        const auto s {poly_coeffs(kind, k, m)};
        auto& r {g[(k + 2) % 3]};
        big_set(r, 0U, 0);
        for (int j {}; j < 3; ++j)
        {
            big_set(t[0], j == 0 ? sig * s.a : u128 {j == 1 ? s.b : s.c}, j == 0 ? (e > 0 ? e : 0) : j * q);
            t[0].neg = j == 0 ? s.a_neg != xneg : j == 2;
            big_mul(t[0], g[(k + (j == 2 ? 0U : 1U)) % 3], t[1]);
            big_add(r, t[1], false, r);
        }
    }
    auto& gn {g[(n + 1) % 3]};
    st.a = &gn;
    // w is n!, or 1 for hermite, in t[0] below
    st.w = &t[0];
    // q (n + m), capped far above the most that the integers can hold
    const std::int64_t qn {q * (std::int64_t {n} + (Assoc ? std::int64_t {m} : 0))};
    st.qn = static_cast<int>(qn < (std::int64_t {1} << 24) ? qn : (std::int64_t {1} << 24));
    st.square = false;
    // poly_unknown where the integers do not fit
    st.sign = steps < n || !gn.ok ? poly_unknown : gn.len == 0 ? poly_is_zero : gn.neg ? -1 : 1;
    gn.neg = false;
    big_set(t[0], 1U, 0);
    for (std::uint64_t k {2}; kind != poly_kind::hermite && k <= n; ++k)
    {
        big_mul_u(t[0], k);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (Assoc)
    {
        // G_(n-1) and G_(n-2) are free now: b reuses them
        poly_big* b {&g[n % 3]};
        big_set(t[1], sig, 0);
        big_mul(t[1], t[1], *b);
        big_set(t[1], 1U, 2 * q);
        big_add(t[1], *b, true, t[1]);
        for (unsigned k {1}; k <= m; ++k)
        {
            big_mul_u(*st.a, 2U * k - 1U);
        }
        for (unsigned k {}; k < m / 2U; ++k)
        {
            big_mul(*st.a, t[1], *b);
            const auto c {st.a};
            st.a = b;
            b = c;
        }
        if ((m & 1U) != 0U)
        {
            big_mul(*st.a, *st.a, *b);
            big_mul(*b, t[1], *st.a);
            st.square = true;
        }
    }
}

// The sign of |P| * 10^k10 - bound for bound > 0: a 10^j against bound 10^-j w for j = k10 - q(n + m), or both squared
constexpr auto poly_exact_cmp(poly_exact& st,const u128& bound, const int k10) noexcept -> int
{
    // Three of the five integers are free
    poly_big* all[5] {&st.g[0], &st.g[1], &st.g[2], &st.t[0], &st.t[1]};
    poly_big* f[3] {};
    int nf {};
    for (int i {}; i < 5 && nf < 3; ++i)
    {
        if (all[i] != st.a && all[i] != st.w)
        {
            f[nf++] = all[i];
        }
    }
    const int j {k10 - st.qn};
    big_copy(*st.a, *f[0]);
    big_mul_pow10(*f[0], st.square ? 2 * j : j);
    big_set(*f[1], bound, -j);
    big_mul(*f[1], *st.w, *f[2]);
    if (st.square)
    {
        big_mul(*f[2], *f[2], *f[1]);
    }
    const auto& rhs {st.square ? *f[1] : *f[2]};
    return st.sign != poly_unknown && f[0]->ok && rhs.ok ? big_cmp(*f[0], rhs) : poly_unknown;
}

// The decimal digits of v, for the bounds of poly_exact_setup
constexpr auto poly_digits(std::uint64_t v) noexcept -> int
{
    int d {};
    for (; v != 0U; v /= 10U)
    {
        ++d;
    }
    return d;
}

// st for P(x), or where x = sig * 10^e is tiny, for the lowest term L: P(0), or P'(0) x for an odd P. Returns the side
// of P where L 10^k10 is the bound, or 0 for P itself, and the shift of k10.
template <bool Assoc>
constexpr auto poly_exact_setup(poly_exact& st, const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg,
                                const u128& bound, const int k10, const int p, int& shift) noexcept -> int
{
    shift = 0;
    const int q {-e};
    const bool odd {kind != poly_kind::laguerre && (n & 1U) != 0U};
    const int skew {k10 < 0 ? -k10 : 0};
    // |P(x) - L| < 4 rho |L| for rho = (n + m + 2)^4 x^2, or n |x| for laguerre, and L decides from k10 on where this times
    // 10^k10 is below one over 2^(n + m) 10^(q odd + max(-k10, 0)), as the denominator of L 10^k10 and |L| 10^k10 < 3 * 10^38
    const int need {kind == poly_kind::laguerre ? 40 + poly_digits(n) + p + skew
                    : 40 + 4 * poly_digits(std::uint64_t {n} + m + 2U) + 2 * p + static_cast<int>((std::uint64_t {n} + m) * 302U / 1000U) + 1 + skew};
    if (e >= 0 || (odd || kind == poly_kind::laguerre ? q : 2 * q) < need)
    {
        // The integers of P do not fit where their comparison with bound does not
        poly_exact_init<Assoc>(st, kind, n, m, sig, e, xneg);
        st.sign = st.sign != poly_is_zero && poly_exact_cmp(st, bound, k10) == poly_unknown ? poly_unknown : st.sign;
        return 0;
    }
    if (odd)
    {
        // P'(0) = f P_(n-1)(0) for f = 2n for hermite, and l + m for legendre of degree l = n + m
        poly_exact_init<Assoc>(st, kind, n - 1U, m, 0U, 0, false, sig, kind == poly_kind::hermite ? 2U * std::uint64_t {n} : std::uint64_t {n} + 2U * m);
        shift = -q;
    }
    else
    {
        poly_exact_init<Assoc>(st, kind, n, m, 0U, 0, false);
    }
    st.sign = st.sign == poly_is_zero ? poly_unknown : st.sign;
    // By the differential equation of P at 0, the next term moves P toward zero, or away from zero for laguerre with x < 0
    return kind == poly_kind::laguerre && xneg ? 1 : -1;
}

// P(x) rounded once from an approximation r whose error interval holds a rounding boundary, out of line so that the
// exact integers take no stack in the common path
template <typename T, bool Assoc>
constexpr auto poly_resolve(const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg, const poly_out& r) noexcept -> T
{
    constexpr int p {precision_v<T>};
    const auto h {pow10(static_cast<u128>(37 - p)) * 5U};
    const auto p37 {pow10(static_cast<u128>(37))};
    // The exact integers with x without its trailing zeros
    const auto tz {remove_trailing_zeros(sig)};
    const auto sx {static_cast<u128>(tz.trimmed_number)};
    const int ex {e + static_cast<int>(tz.number_of_removed_zeros)};
    int k10 {-r.e10};
    int shift {};
    poly_exact st {};
    const int side {poly_exact_setup<Assoc>(st, kind, n, m, sx, ex, xneg, r.bound, k10, p, shift)};
    if (st.sign == poly_is_zero)
    {
        return T {0};
    }

    // The boundaries are the multiples of h from lo to hi, and v the last one below |P| * 10^k10, or lo - 1 where it is
    // below all. A value below 10^37 has the boundaries 10 times closer, thus lo - h would not do.
    auto lo {r.bound};
    auto hi {r.bound};
    bool neg {r.neg};
    if (r.units >= h / 2U)
    {
        hi = (r.d + r.units) / h * h;
        lo = r.d >= p37 + r.units ? (r.d - r.units + h - 1U) / h * h : p37;
        if (r.d < p37 + r.units)
        {
            // Next to a root the sign and the magnitude come from the exact integers: the least k10 with
            // |P| * 10^k10 >= 10^37
            if (st.sign == poly_unknown)
            {
                return fx_round<T>(r.d, r.e10, r.neg);
            }
            neg = st.sign < 0;
            int k_lo {k10};
            int k_hi {40 - etiny_v<T>};
            while (k_lo < k_hi)
            {
                const int k {k_lo + (k_hi - k_lo) / 2};
                const int ck {poly_exact_cmp(st, p37, k + shift)};
                k_hi = ck >= 0 ? k : k_hi;
                k_lo = ck >= 0 ? k_lo : k + 1;
            }
            hi = k_lo == k10 ? hi : pow10(static_cast<u128>(38));
            k10 = k_lo;
        }
    }
    auto v {lo - 1U};
    while (lo <= hi)
    {
        const auto mid {lo + (hi - lo) / (2U * h) * h};
        int c {poly_exact_cmp(st, mid, k10 + shift)};
        if (c == poly_unknown)
        {
            // The integers do not decide, so the approximation picks the side
            c = r.d < mid ? -1 : 1;
        }
        else if (c == 0 && side != 0)
        {
            c = side;
        }
        if (c == 0)
        {
            return T {mid, -k10, neg};
        }
        v = c > 0 ? mid : v;
        lo = c > 0 ? mid + h : lo;
        hi = c > 0 ? hi : mid - h;
    }
    // |P| * 10^k10 is above v, and below the next boundary
    return fx_round<T>(v, -k10, neg);
}

// a * b in r, which may be a, or false where it can reach 2^126
constexpr auto poly_mul_fits(const u128& a, const u128& b, u128& r) noexcept -> bool
{
    const bool fits {256 - countl_zero(a) - countl_zero(b) <= 126};
    r = a * b;
    return fits;
}

// m1 - m2 or m1 + m2 by the signs n1 and n2, with the sign of the result in n1
constexpr auto poly_add_signed(u128& m1, bool& n1, const u128& m2, const bool n2) noexcept -> void
{
    const bool swap {n1 != n2 && m1 < m2};
    m1 = n1 == n2 ? m1 + m2 : swap ? m2 - m1 : m1 - m2;
    n1 = swap ? n2 : n1;
}

// P(x) rounded once in out where its exact integers G_k of poly_exact fit in 128 bits, as for a short x and a low degree:
// P = G_n / (n! 10^(qn)), from p + 1 digits of the quotient or more, and a last digit for the rest
template <bool Assoc, typename T, typename Sig>
constexpr auto poly_exact_small(const poly_kind kind, const unsigned n, const unsigned m, const Sig sig, const int e, const bool xneg, T& out) noexcept -> bool
{
    constexpr int p {precision_v<T>};
    const auto tz {remove_trailing_zeros(sig)};
    const int ex {e + static_cast<int>(tz.number_of_removed_zeros)};
    const int q {ex < 0 ? -ex : 0};
    u128 s {tz.trimmed_number};
    const auto tq {pow10(static_cast<u128>(q < 38 ? q : 38))};
    const unsigned nm {Assoc ? n + m : n};
    if ((kind != poly_kind::hermite && n > 20U) || ex > 38 || q > 19 || (Assoc && m % 2U != 0U) ||
        (ex > 0 && !poly_mul_fits(s, pow10(static_cast<u128>(ex)), s)))
    {
        return false;
    }
    const std::uint64_t w {kind == poly_kind::hermite ? 1U : poly_factorial(n)};
    // G_n has near n times the bits of s or 10^q, and those of n!, of 4^n for hermite and of (2m - 1)!!
    const auto bits {std::uint64_t {nm} * static_cast<std::uint64_t>(128 - countl_zero(s > tq ? s : tq)) +
                     (kind == poly_kind::hermite ? 2U * std::uint64_t {n} : static_cast<std::uint64_t>(64 - countl_zero(w))) +
                     (Assoc ? std::uint64_t {m} * static_cast<std::uint64_t>(64 - countl_zero(std::uint64_t {2U * m})) : 0U)};
    if (bits > 136U)
    {
        return false;
    }
    const auto tq2 {tq * tq};
    u128 g0 {};
    u128 g1 {1U};
    bool n0 {};
    bool n1 {};
    for (unsigned k {}; k < n; ++k)
    {
        // G_(k+1) = a s G_k + b 10^q G_k - c 10^2q G_(k-1)
        const auto c {poly_coeffs(kind, k, m)};
        u128 u {};
        u128 t1 {};
        u128 t2 {};
        u128 t3 {};
        if (!poly_mul_fits(s, u128 {c.a}, u) || !poly_mul_fits(u, g1, t1) || !poly_mul_fits(tq, u128 {c.b}, u) ||
            !poly_mul_fits(u, g1, t2) || !poly_mul_fits(tq2, u128 {c.c}, u) || !poly_mul_fits(u, g0, t3))
        {
            return false;
        }
        bool nr {(c.a_neg != xneg) != n1};
        poly_add_signed(t1, nr, t2, n1);
        poly_add_signed(t1, nr, t3, !n0);
        g0 = g1;
        n0 = n1;
        g1 = t1;
        n1 = nr;
    }
    for (unsigned k {1}; Assoc && k <= m; ++k)
    {
        // G_n takes (2m - 1)!! (10^2q - s^2)^(m/2) for an even m, and qn is q (n + m): 2k - 1 for each k, and for each
        // other k 10^2q - s^2 with s < 10^q as |x| < 1
        if (!poly_mul_fits(g1, u128 {2U * k - 1U}, g1) || ((k & 1U) == 0U && !poly_mul_fits(g1, tq2 - s * s, g1)))
        {
            return false;
        }
    }
    if (g1 == 0U)
    {
        out = T {0};
        return true;
    }
    // G_n / w to p + 1 digits or more by long division, 19 digits at a time: as the rest r is below w, r 10^j / w is
    // below 2^64
    std::uint64_t r {};
    const auto low {boost::int128::detail::udiv_2by1(g1.high % w, g1.low, w, r)};
    u128 d {g1.high / w, low};
    int k10 {};
    for (int nd {num_digits(d)}; nd <= p; nd = num_digits(d))
    {
        const int j {p + 1 - nd < 19 ? p + 1 - nd : 19};
        const auto t {u128 {r} * pow10(static_cast<std::uint64_t>(j))};
        d = d * pow10(static_cast<u128>(j)) + boost::int128::detail::udiv_2by1(t.high, t.low, w, r);
        k10 += j;
    }
    // The constructor sees a last 0 as exact and a last 5 as a tie, so a rest makes it 1 or 6
    if (r != 0U && (d.high % 5U + d.low % 5U) % 5U == 0U)
    {
        d += 1U;
    }
    out = T {d, -k10 - q * static_cast<int>(nm), n1};
    return true;
}

// P(x) rounded once in the current mode. Where the error interval covers several rounding boundaries, the approximation
// runs again in two more words, and next to the boundaries a search with the exact integers decides.
template <typename T, bool Assoc = false>
constexpr auto poly_fixed(const poly_kind kind, const unsigned n, const unsigned m, const T x) noexcept -> T
{
    constexpr int p {precision_v<T>};
    constexpr int words {p <= 16 ? 2 : 3};
    int e {};
    const auto sig_x {frexp10(x, &e)};
    const auto sig {static_cast<u128>(sig_x)};
    const bool xneg {signbit(x)};
    // An odd hermite or legendre, and assoc_legendre for an odd l - m, is +0 at x = 0
    if (sig == 0U && kind != poly_kind::laguerre && (n & 1U) != 0U)
    {
        return T {0};
    }
    // A first try where the error interval holds no rounding boundary, below a degree where it seldom decides, and in one
    // word first for decimal32_t
    constexpr int lo10 {etiny_v<T>};
    constexpr int hi10 {std::numeric_limits<T>::max_exponent10 - p};
    using sig_type = typename T::significand_type;
    BOOST_DECIMAL_IF_CONSTEXPR (p <= 7)
    {
        if (n + m <= (kind == poly_kind::laguerre ? 11U : 22U))
        {
            const auto q1 {poly_try<1, Assoc>(kind, n, m, sig, e, xneg, p, lo10, hi10)};
            if (q1.ok)
            {
                return fx_round_kept<T>(static_cast<sig_type>(q1.d), q1.e10, q1.neg);
            }
        }
    }
    T out {};
    if (poly_exact_small<Assoc>(kind, n, m, sig_x, e, xneg, out))
    {
        return out;
    }
    if (n + m <= (kind == poly_kind::laguerre ? 22U : kind == poly_kind::hermite ? 64U : 48U))
    {
        const auto q {poly_try<words, Assoc>(kind, n, m, sig, e, xneg, p, lo10, hi10)};
        if (q.ok)
        {
            return fx_round_kept<T>(static_cast<sig_type>(q.d), q.e10, q.neg);
        }
    }
    const auto h {pow10(static_cast<u128>(37 - p)) * 5U};
    auto r {poly_eval<words, Assoc>(kind, n, m, sig, e, xneg, p)};
    if (r.amb && r.units >= h / 2U)
    {
        r = poly_eval<words + 2, Assoc>(kind, n, m, sig, e, xneg, p);
    }
    // |P| is below 10^(e10 + 39), and at least 10^(e10 + 37) where d - units is: past the range of T, the constructor
    // overflows or underflows in the mode
    const auto p37 {pow10(static_cast<u128>(37))};
    if ((r.e10 + 37 > std::numeric_limits<T>::max_exponent10 && r.d >= p37 + r.units) || r.e10 + 40 < etiny_v<T>)
    {
        return T {1U, (r.e10 > 0 ? 2 : -2) * std::numeric_limits<T>::max_exponent10, r.neg};
    }
    return r.amb ? poly_resolve<T, Assoc>(kind, n, m, sig, e, xneg, r) : fx_round<T>(r.d, r.e10, r.neg);
}

} // namespace poly_detail
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_POLYNOMIAL_IMPL_HPP
