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

// R_(k+1) = (a_k x + b_k) R_k - c_k R_(k-1) from R_0 = 1: hermite (a, b, c) = (2, 0, 2k) and P = R_n, laguerre (-1, 2k + m + 1,
// k (k + m)) and P = R_n / n!, legendre of order m (2k + 2m + 1, 0, k (k + 2m)) and P = (2m - 1)!! (1 - x^2)^(m/2) R_n / n!
enum class poly_kind { hermite, legendre, laguerre };

struct poly_step
{
    std::uint64_t a;
    bool a_neg;
    std::uint64_t b;
    std::uint64_t c;
};

constexpr auto poly_coeffs(const poly_kind kind, const std::uint64_t k, const std::uint64_t m) noexcept -> poly_step
{
    return kind == poly_kind::hermite ? poly_step {2U, false, 0U, 2U * k}
         : kind == poly_kind::legendre ? poly_step {2U * k + 2U * m + 1U, false, 0U, k * (k + 2U * m)}
         : poly_step {1U, true, 2U * k + m + 1U, k * (k + m)};
}

// m * 2^e with m in [1, 2), or zero
template <int N>
struct poly_float
{
    fx<N> m;
    int e;
    bool neg;
};

template <int N>
constexpr auto pf_zero(const poly_float<N>& a) noexcept -> bool
{
    return fx_lead(a.m) == 64 * N;
}

template <int N>
constexpr auto pf_norm(const fx<N>& m, const int e, const bool neg) noexcept -> poly_float<N>
{
    const int l {fx_lead(m)};
    return l == 64 * N ? poly_float<N> {m, 0, neg} : poly_float<N> {fx_shift(m, l - 1), e - (l - 1), neg};
}

// v as a float: v * 2^-126 in N + 2 words, normalized and truncated to N
template <int N>
constexpr auto pf_from(const u128& v, const bool neg = false) noexcept -> poly_float<N>
{
    fx<N + 2> t {};
    t.w[N] = v.low;
    t.w[N + 1] = v.high;
    const int l {fx_lead(t)};
    return l == 64 * (N + 2) ? poly_float<N> {fx<N> {}, 0, neg} : poly_float<N> {fx_top<N>(fx_shift(t, l - 1)), 126 - (l - 1), neg};
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
    return pf_zero(a) || pf_zero(b) ? poly_float<N> {} : pf_half(fx_mul(a.m, b.m), a.e + b.e, a.neg != b.neg);
}

// a * c for c > 0: a * c / 2^s is in [1/2, 2) for c in [2^(s - 1), 2^s)
template <int N>
constexpr auto pf_mul_u(const poly_float<N>& a, const std::uint64_t c) noexcept -> poly_float<N>
{
    if (pf_zero(a) || c == 0U)
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
    if (pf_zero(b))
    {
        return a;
    }
    if (pf_zero(a))
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
        b = pf_mul(b, b);
        n >>= 1U;
    }
    return r;
}

// 10^k by squaring, with a relative error below 2 log2(|k|) + 2 units
template <int N>
constexpr auto pf_pow10(const int k) noexcept -> poly_float<N>
{
    const auto n {static_cast<unsigned>(k < 0 ? -k : k)};
    if (n <= 38U)
    {
        const auto t {pf_from<N>(pow10(static_cast<u128>(n)))};
        return k < 0 ? pf_recip(t) : t;
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
    auto y {pf_add(one, pf_mul(x, x), true)};
    if (q <= 38)
    {
        const auto t {pow10(static_cast<u128>(q))};
        y = pf_mul(pf_mul(pf_from<N>(t - sig), pf_from<N>(t + sig)), pf_pow10<N>(-2 * q));
    }
    auto r {one};
    for (unsigned k {1}; k <= m; ++k)
    {
        r = pf_mul_u(r, 2U * k - 1U);
    }
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

template <int N>
constexpr auto pb_of(const poly_float<N>& a) noexcept -> poly_bound
{
    return pf_zero(a) ? poly_bound {0U, 0} : poly_bound {(a.m.w[N - 1] >> 32U) + 1U, a.e - 30};
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

// The digits of a result: d * 10^e10 with 38 or 39 digits. amb is set where the error can reach a rounding boundary,
// bound is the boundary next to d, and zero is set for a zero value.
struct poly_out
{
    u128 d;
    int e10;
    bool neg;
    bool amb;
    u128 bound;
    bool zero;
    u128 units;
};

// P(x) for x = sig * 10^e in N words with a bound on the error of each step: the error of R_(k+1) is below |a x + b|
// times that of R_k, plus c times that of R_(k-1), plus the rounding of the step
template <int N, bool Assoc>
constexpr auto poly_eval(const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg, const int p) noexcept -> poly_out
{
    const poly_bound u {1U, -(64 * N - 8)};
    const auto x {pf_mul(pf_from<N>(sig, xneg), pf_pow10<N>(e))};
    const auto bx {pb_of(x)};
    poly_float<N> r0 {};
    poly_float<N> r1 {fx_one<N>(), 0, false};
    poly_bound bnd0 {0U, 0};
    poly_bound bnd1 {1U, 0};
    // With big = |a x| + |b| and a margin for the error of x, the rounding of a step is below u (44 big |R_k| + 3 c |R_(k-1)|),
    // and the error of R_k below 47 k u B_k for B_(k+1) = big B_k + c B_(k-1), which bounds |R_k|
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
        const auto t1 {pf_mul(coef, r1)};
        const auto t2 {pf_mul_u(r0, s.c)};
        const auto r2 {pf_add(t1, t2, true)};

        const auto big {pb_add(pb_mul(pb_make(s.a, 0), bxm), pb_make(s.b, 0))};
        const auto bnd2 {pb_add(pb_mul(big, bnd1), pb_mul(pb_make(s.c, 0), bnd0))};
        r0 = r1;
        r1 = r2;
        bnd0 = bnd1;
        bnd1 = bnd2;
    }
    auto value {r1};
    auto err {pb_mul(pb_mul(u, pb_make(47U * n, 0)), bnd1)};
    if (kind != poly_kind::hermite)
    {
        poly_float<N> f {fx_one<N>(), 0, false};
        for (std::uint64_t k {2}; k <= n; ++k)
        {
            f = pf_mul_u(f, k);
        }
        auto inv {pf_recip(f)};
        std::uint64_t rho {n + 4U};
        BOOST_DECIMAL_IF_CONSTEXPR (Assoc)
        {
            inv = pf_mul(inv, pf_legendre_scale(m, sig, e, x));
            rho += 16U * m + 32U;
        }
        value = pf_mul(r1, inv);
        err = pb_add(pb_mul(pb_mul(err, pb_of(inv)), pb_make((UINT64_C(1) << 40U) + 1U, -40)), pb_mul(pb_mul(u, pb_make(rho, 0)), pb_of(value)));
    }

    poly_out out {0U, 0, value.neg, true, 0U, false, 0U};
    if (pf_zero(value))
    {
        out.zero = true;
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

// P(x) in exact integers: with x = sig * 10^e, q = max(-e, 0) and s = sig * 10^max(e, 0), G_k = R_k 10^(qk) gives
// G_(k+1) = (a_k s + b_k 10^q) G_k - c_k 10^(2q) G_(k-1). a is |G_n|, or for legendre of order m, |L| for
// L = G_n (2m - 1)!! M^(m/2) with M = 10^2q - s^2, or L^2 M for an odd m (square). w is n!, or 1 for hermite.
// sign is that of P, or poly_is_zero, or poly_unknown where the integers do not fit.
template <bool Assoc>
struct poly_exact
{
    poly_big g[3];
    poly_big t[2];
    poly_big* a;
    poly_big* w;
    int qn;
    bool square;
    int sign;
};

// G_0 = f * fu in place of 1 multiplies P by f * fu
template <bool Assoc>
constexpr auto poly_exact_init(poly_exact<Assoc>& st, const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg,
                               const u128& f = 1U, const std::uint64_t fu = 1U) noexcept -> void
{
    const int q {e < 0 ? -e : 0};
    auto& g {st.g};
    auto& t {st.t};
    big_set(g[0], 0U, 0);
    big_set(g[1], f, 0);
    big_mul_u(g[1], fu);
    // G_k in g[(k + 1) % 3], and t for the terms of a step
    for (unsigned k {}; k < n && g[(k + 1) % 3].ok; ++k)
    {
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
    st.w = &t[0];
    st.qn = q * static_cast<int>(Assoc ? n + m : n);
    st.square = false;
    st.sign = !gn.ok ? poly_unknown : gn.len == 0 ? poly_is_zero : gn.neg ? -1 : 1;
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
template <bool Assoc>
constexpr auto poly_exact_cmp(poly_exact<Assoc>& st, const u128& bound, const int k10) noexcept -> int
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
    *f[0] = *st.a;
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

// st for P(x), or where its integers do not fit and x = sig * 10^e is tiny, for the lowest term L: P(0), or P'(0) x for
// an odd P, with P'(0) = f P_(n-1)(0) for f = 2n for hermite, and l + m for legendre of degree l = n + m. The ratios of
// the next coefficients bound |P(x) - L| below 4 rho |L| for rho = (n + m + 2)^4 x^2, or n |x| for laguerre. As |L| 10^k10
// < 3 * 10^38 and the denominator of L 10^k10 is below 2^(n + m) 10^(q odd + max(-k10, 0)), 4 rho |L| 10^k10 below one over
// it lets L decide, from k10 on. Where L 10^k10 is the bound, the next term moves P toward zero, or away from zero for
// laguerre with x < 0, by the differential equation of P at 0. Returns that side, or 0 for P itself, and the shift of k10.
// The integers of P do not fit where their comparison with bound does not.
template <bool Assoc>
constexpr auto poly_exact_setup(poly_exact<Assoc>& st, const poly_kind kind, const unsigned n, const unsigned m, const u128& sig, const int e, const bool xneg,
                                const u128& bound, const int k10, const int p, int& shift) noexcept -> int
{
    poly_exact_init(st, kind, n, m, sig, e, xneg);
    st.sign = st.sign != poly_is_zero && poly_exact_cmp(st, bound, k10) == poly_unknown ? poly_unknown : st.sign;
    shift = 0;
    const int q {-e};
    const bool odd {kind != poly_kind::laguerre && (n & 1U) != 0U};
    const int skew {k10 < 0 ? -k10 : 0};
    const int need {kind == poly_kind::laguerre ? 40 + poly_digits(n) + p + skew
                    : 40 + 4 * poly_digits(std::uint64_t {n} + m + 2U) + 2 * p + static_cast<int>((std::uint64_t {n} + m) * 302U / 1000U) + 1 + skew};
    if (st.sign != poly_unknown || e >= 0 || (odd ? q : 2 * q) < need)
    {
        return 0;
    }
    if (odd)
    {
        poly_exact_init(st, kind, n - 1U, m, 0U, 0, false, sig, kind == poly_kind::hermite ? 2U * std::uint64_t {n} : std::uint64_t {n} + 2U * m);
        shift = -q;
    }
    else
    {
        poly_exact_init(st, kind, n, m, 0U, 0, false);
    }
    st.sign = st.sign == poly_is_zero ? poly_unknown : st.sign;
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
    poly_exact<Assoc> st {};
    const int rule {poly_exact_setup(st, kind, n, m, sx, ex, xneg, r.bound, k10, p, shift)};
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
        else if (c == 0 && rule != 0)
        {
            c = rule;
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


// P(x) rounded once in the current mode. Where the error interval of the approximation covers several rounding boundaries,
// it runs again in two more words. Next to the boundaries a search with the exact integers decides. Where neither the
// integers nor the tiny x decide, for a high degree, the approximation picks the side.
template <typename T, bool Assoc = false>
constexpr auto poly_fixed(const poly_kind kind, const unsigned n, const unsigned m, const T x) noexcept -> T
{
    constexpr int p {precision_v<T>};
    constexpr int words {p <= 16 ? 2 : 3};
    int e {};
    const auto sig {static_cast<u128>(frexp10(x, &e))};
    const bool xneg {signbit(x)};
    const auto h {pow10(static_cast<u128>(37 - p)) * 5U};
    auto r {poly_eval<words, Assoc>(kind, n, m, sig, e, xneg, p)};
    if (r.amb && r.units >= h / 2U)
    {
        r = poly_eval<words + 2, Assoc>(kind, n, m, sig, e, xneg, p);
    }
    if (r.zero)
    {
        return T {0};
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
