// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/impl/gamma_impl.hpp>

namespace boost {
namespace decimal {

namespace detail {

namespace beta_detail {

using gamma_detail::gamma_float;
using gamma_detail::gamma_exp;
using gamma_detail::gamma_parts;

// x + y < 1 for x = sx * 10^ex and y = sy * 10^ey with p digits, as v * 2^n2 * 10^q with q = max(ex, ey) + p: the
// halves of x / 10^q and y / 10^q add up below 1, and the smaller loses only its digits below 2^-64N of the sum
template <int N>
constexpr auto beta_small_sum(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const int p) noexcept -> gamma_float<N>
{
    const int q {(ex > ey ? ex : ey) + p};
    const auto h {fx_add(fx_shr(fx_from_decimal<N>(sx, q - ex), 1U), fx_shr(fx_from_decimal<N>(sy, q - ey), 1U))};
    unsigned s {};
    const auto v {fx_normalize(h, s)};
    return {v, 1 - static_cast<int>(s), q};
}

// acc * f, or false where it does not fit in 256 bits
constexpr auto beta_mul(u256& acc, const boost::int128::uint128_t& f) noexcept -> bool
{
    std::uint64_t r[6] {};
    const std::uint64_t fw[2] {f.low, f.high};
    for (int j {}; j < 2; ++j)
    {
        std::uint64_t carry {};
        for (int i {}; i < 4; ++i)
        {
            const auto t {boost::int128::uint128_t {acc[static_cast<std::size_t>(i)]} * fw[j] + r[i + j] + carry};
            r[i + j] = t.low;
            carry = t.high;
        }
        r[4 + j] = carry;
    }
    acc = u256 {r[3], r[2], r[1], r[0]};
    return (r[4] | r[5]) == 0U;
}

// acc * 10^d, or false where it does not fit in 256 bits
constexpr auto beta_mul_pow10(u256& acc, int d) noexcept -> bool
{
    for (; d > 0; d -= 38)
    {
        if (!beta_mul(acc, pow10(static_cast<boost::int128::uint128_t>(d < 38 ? d : 38))))
        {
            return false;
        }
    }
    return true;
}

// The other argument z = a * 10^-k of beta(n, z) for an integer n in [2, 35], or n = 0 for none
struct beta_integer
{
    boost::int128::uint128_t a;
    int k;
    unsigned n;
};

// Whether beta(n, z) = (n - 1)! / (z (z + 1) ... (z + n - 1)) is r * 10^q exactly: r * prod (a + i 10^k) * 10^q =
// (n - 1)! * 10^(kn), in 256 bits. Where it does not fit, the result is not exact. From k = 36 on, the parts of the
// factors prime to 10 multiply past (n - 1)!, and below, the factors fit in 128 bits.
constexpr auto beta_is_exact(const beta_integer& b, const boost::int128::uint128_t& r, const int q) noexcept -> bool
{
    if (b.k >= 36)
    {
        return false;
    }
    u256 lhs {r};
    u256 rhs {gamma_detail::gamma_at(b.n)};
    const auto step {pow10(static_cast<boost::int128::uint128_t>(b.k))};
    for (unsigned i {}; i < b.n; ++i)
    {
        if (!beta_mul(lhs, b.a + step * i))
        {
            return false;
        }
    }
    const int d {q - b.k * static_cast<int>(b.n)};
    return (d >= 0 ? beta_mul_pow10(lhs, d) : beta_mul_pow10(rhs, -d)) && lhs == rhs;
}

// Whether digits with a relative error below 2^-bits can round another way than the exact value: the boundaries
// are the multiples of half a unit in the last of p digits
template <typename digits_type>
constexpr auto beta_near_boundary(const digits_type& digits, const int p, const int bits) noexcept -> bool
{
    const auto m {pow10(static_cast<digits_type>(num_digits(digits) - p - 1)) * 5U};
    const auto r {digits % m};
    const auto w {(digits >> bits) + 1U};
    return r < w || r > m - w;
}

// beta(m, M) = 1 / m + d for a tiny m = sm * 10^em, with d > 0 for M < 1, d < 0 for M > 1 and |d| < 10^bound, in out. An
// exact 1 / m moves to the side of d, else 1 / m + d rounds as 1 / m where |d| is below its gap to the boundaries, or false.
template <typename T>
constexpr auto beta_tiny(const boost::int128::uint128_t& sm, const int em, const bool above, const int bound, T& out) noexcept -> bool
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto r {impl::div_mod(umul256(pow10(static_cast<boost::int128::uint128_t>(p)), pow10(static_cast<boost::int128::uint128_t>(p))), sm)};
    const auto q {static_cast<boost::int128::uint128_t>(r.quotient)};
    const auto rem {static_cast<boost::int128::uint128_t>(r.remainder)};
    const int e {-em - 2 * p};
    if (rem == 0U)
    {
        out = T {above ? q * 10U + 1U : q * 10U - 1U, e - 1};
        return bound < e;
    }
    // |d| < 10^bound <= 10^e * gap / sm, with gap = min(rem, sm - rem) >= 10^(digits of gap - 1)
    const auto gap {rem < sm - rem ? rem : sm - rem};
    out = T {q * 10U + 5U, e - 1};
    return bound + num_digits(sm) <= e + num_digits(gap) - 1;
}

// beta(m, M) = gamma(m) M^-m e^-(F + d) for M = sM * 10^eM >= 10^15 and 0 < m = sm * 10^em < 10^4, with r = m / M:
// F = (M + m - 1/2) ln(1 + r) - m = r m A - ln(1 + r) / 2 for A = sum (-r)^(j - 1) / (j (j + 1)), d = -r (1 - r) / (12 M)
template <typename T, int N>
constexpr auto beta_large(const boost::int128::uint128_t& sm, const int em, const boost::int128::uint128_t& sM, const int eM) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto gm {gamma_detail::gamma_positive<N>(sm, em, p)};
    const auto mw {gamma_detail::gamma_wide<N>(sm, em)};
    // m < 10^4 < 2^14 and ln(M) < 2^14, so m / 2^14 and ln(M) / 2^14 fit below 1
    const auto m14 {fx_wide_shr(mw, 14U)};
    const int nd {num_digits(sM)};
    const auto ln_m {log_detail::log_wide_of<N>(sM, nd, eM, log_detail::log_base::e)};

    // m ln(M) = m * (ln(M) / 2^14) * 2^14
    const auto m_ln {fx_wide_shl(fx_wide_mul(mw, fx_wide_shr(ln_m, 14U)), 14U)};

    // 1 / M, r = m / M and a = r m
    const auto inv {gamma_detail::gamma_recip<N>(sM, nd, eM)};
    const auto r {fx_shl(fx_mul(inv, m14), 14U)};
    const auto a {fx_shl(fx_mul(r, m14), 14U)};

    // A and ln(1 + r) / 2 as alternating sums, which stay between 0 and their first terms
    fx<N> sum_a {};
    fx<N> half_ln {};
    auto pw {fx_one<N>()};
    for (std::uint64_t j {1}; fx_countl_zero(pw) < 64 * N; ++j)
    {
        const auto ta {fx_div_small(pw, j * (j + 1U))};
        pw = fx_mul(pw, r);
        const auto tb {fx_div_small(pw, 2U * j)};
        sum_a = (j & 1U) != 0U ? fx_add(sum_a, ta) : fx_sub(sum_a, ta);
        half_ln = (j & 1U) != 0U ? fx_add(half_ln, tb) : fx_sub(half_ln, tb);
    }

    // st = st(m) - m ln(M) - F - d = st(m) - m ln(M) - a A - q r + ln(1 + r) / 2 + q, with q = r / (12 M)
    const auto q {fx_div_small(fx_mul(r, inv), 12U)};
    auto st {fx_wide_add(gm.st, m_ln, true)};
    st = fx_wide_add(st, fx_widen(fx_add(fx_mul(a, sum_a), fx_mul(q, r))), true);
    st = fx_wide_add(st, fx_widen(fx_add(half_ln, q)));
    const auto res {gamma_detail::gamma_exp_digits<T>(st, gamma_detail::gamma_frecip(gm.c))};
    if (res.kind != exp_detail::exp_kind::digits)
    {
        return exp_detail::exp10_beyond<T, N>(res.kind, false);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N < 3)
    {
        if (beta_near_boundary(res.digits, p, 72))
        {
            return beta_large<T, (N < 3 ? N + 1 : 3)>(sm, em, sM, eM);
        }
    }
    return fx_round<T>(res.digits, res.exp10, false);
}

// beta(x, y) = c * e^st for x + y < 20, with gamma(z) = c(z) * e^st(z) from gamma_poly_up
template <int N>
constexpr auto beta_poly(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const int p) noexcept
    -> std::enable_if_t<(N < 3), gamma_exp<N>>
{
    const auto wx {gamma_detail::gamma_wide<N>(sx, ex)};
    const auto wy {gamma_detail::gamma_wide<N>(sy, ey)};
    const auto w {fx_wide_add(wx, wy)};
    const auto gx {gamma_detail::gamma_poly_up(wx, gamma_detail::gamma_whole(wx) == 0U ? gamma_detail::gamma_fdec<N>(sx, p, ex) : gamma_float<N> {})};
    const auto gy {gamma_detail::gamma_poly_up(wy, gamma_detail::gamma_whole(wy) == 0U ? gamma_detail::gamma_fdec<N>(sy, p, ey) : gamma_float<N> {})};
    const auto gs {gamma_detail::gamma_poly_up(w, gamma_detail::gamma_whole(w) == 0U ? beta_small_sum<N>(sx, ex, sy, ey, p) : gamma_float<N> {})};
    const auto c {gamma_detail::gamma_fmul(gamma_detail::gamma_fmul(gx.c, gy.c), gamma_detail::gamma_frecip(gs.c))};
    return {fx_wide_add(fx_widen(fx_add(gx.st, gy.st)), fx_widen(gs.st), true), c};
}

template <int N>
constexpr auto beta_poly(const boost::int128::uint128_t&, int, const boost::int128::uint128_t&, int, int) noexcept
    -> std::enable_if_t<(N == 3), gamma_exp<N>>
{
    return {};
}

// beta(x, y) for x + y < 20 for decimal64_t in one word, rounded once in out where it cannot be next to a rounding
// boundary: the two products, the reciprocal and gamma_float_digits add 8 units of 2^-62 to the errors of the gammas
template <typename T>
constexpr auto beta_try(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, T& out) noexcept -> bool
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto wx {gamma_detail::gamma_wide<1>(sx, ex)};
    const auto wy {gamma_detail::gamma_wide<1>(sy, ey)};
    const auto w {fx_wide_add(wx, wy)};
    const auto nx {gamma_detail::gamma_whole(wx)};
    const auto ny {gamma_detail::gamma_whole(wy)};
    const auto n {gamma_detail::gamma_whole(w)};
    const auto gx {gamma_detail::gamma_poly_one(wx, sx, ex)};
    const auto gy {gamma_detail::gamma_poly_one(wy, sy, ey)};
    const auto gs {gamma_detail::gamma_poly_up<1, 16>(w, n == 0U ? gamma_detail::gamma_ftop(beta_small_sum<2>(sx, ex, sy, ey, p)) : gamma_float<1> {})};
    const auto r {gamma_detail::gamma_float_digits(gamma_detail::gamma_fmul(gamma_detail::gamma_fmul(gx, gy), gamma_detail::gamma_frecip(gs.c)))};
    return fx_try_round(r.digits, r.exp10, false, ((r.digits >> 58U) + 1U) * (gamma_detail::gamma_poly_error(nx) + gamma_detail::gamma_poly_error(ny) + gamma_detail::gamma_poly_error(n) + 8U) / 16U + 2U, out);
}

// ln gamma(z) added to st for z = sig * 10^e >= 10, else gamma(z) multiplied into c, and the error of the step in err
constexpr auto beta_add_gamma(const fx_wide<2>& z, const boost::int128::uint128_t& sig, const int e, fx_wide<2>& st, gamma_float<1>& c, std::uint64_t& err) noexcept -> void
{
    const auto n {gamma_detail::gamma_whole(z)};
    if (n >= 10U)
    {
        st = fx_wide_add(st, gamma_detail::gamma_stirling_one(z));
        err += 2U;
    }
    else
    {
        c = gamma_detail::gamma_fmul(c, gamma_detail::gamma_poly_one(fx_wide<1> {{z.w[1], z.w[2]}}, sig, e));
        err += gamma_detail::gamma_poly_error(n) + 2U;
    }
}

// beta(x, y) for x + y >= 20 for decimal64_t in one word as in tgamma, rounded once in out where it cannot be next to a
// rounding boundary: its error is below 48 units of 2^-62, 2 more per other Stirling tail and those of gamma_poly_one
template <typename T>
constexpr auto beta_try_large(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, T& out) noexcept -> bool
{
    const auto wx {gamma_detail::gamma_wide<2>(sx, ex)};
    const auto wy {gamma_detail::gamma_wide<2>(sy, ey)};
    auto st {fx_wide_negate(gamma_detail::gamma_stirling_one(fx_wide_add(wx, wy)))};
    gamma_float<1> c {fx_one<1>(), 0, 0};
    std::uint64_t err {48U};
    beta_add_gamma(wx, sx, ex, st, c, err);
    beta_add_gamma(wy, sy, ey, st, c, err);
    const auto r {gamma_detail::gamma_exp_digits<T, 2, false, 1>(st, gamma_float<2> {fx_extend<2>(c.v), c.n2, c.n10})};
    if (r.kind != exp_detail::exp_kind::digits)
    {
        // t in two words is far past the bound where the result underflows or overflows, as in beta_words
        out = exp_detail::exp10_beyond<T, 2>(r.kind, false);
        return true;
    }
    return fx_try_round(r.digits, r.exp10, false, (r.digits >> 57U) * err / 32U + 2U, out);
}

// x = a + eps for an integer a >= 1 and |eps| < 10^-25, with |eps| * 10^25 = s * 10^-k, for decimal128_t
struct beta_near
{
    boost::int128::uint128_t s;
    unsigned a;
    int k;
    bool neg;
    bool ok;
};

constexpr auto beta_near_of(const boost::int128::uint128_t& sig, const int e) noexcept -> beta_near
{
    const int d {-e};
    if (d < 26 || d > 34)
    {
        return {0U, 0U, 0, false, false};
    }
    const auto ten {pow10(static_cast<boost::int128::uint128_t>(d))};
    const auto a {(sig + ten / 2U) / ten};
    const auto whole {a * ten};
    const bool neg {sig < whole};
    const auto s {neg ? whole - sig : sig - whole};
    const bool ok {a != 0U && a < 30U && s < pow10(static_cast<boost::int128::uint128_t>(d - 25))};
    return {s, static_cast<unsigned>(a), d - 25, neg, ok};
}

// floor(m * g * 2^40)
constexpr auto beta_scale40(const boost::int128::uint128_t& m, const fx<3>& g) noexcept -> u256
{
    std::uint64_t r[5] {};
    const std::uint64_t mw[2] {m.low, m.high};
    for (int i {}; i < 2; ++i)
    {
        std::uint64_t carry {};
        for (int j {}; j < 3; ++j)
        {
            const auto t {static_cast<boost::int128::uint128_t>(mw[i]) * g.w[j] + r[i + j] + carry};
            r[i + j] = t.low;
            carry = t.high;
        }
        r[i + 3] = carry;
    }
    return u256 {0U, r[4] >> 22U, (r[4] << 42U) | (r[3] >> 22U), (r[3] << 42U) | (r[2] >> 22U)};
}

// beta(a + eps, b + db) = P / Q e^S with S = -(A eps + B db) + (A2 eps^2 + B2 db^2) / 2 - eps db psi'(a + b), for
// A = sum 1 / i, A2 = sum 1 / i^2 from i = a, B and B2 from b, to a + b - 1. The rest is far below one unit.
template <typename T>
constexpr auto beta_near_integers(const beta_near& x, const beta_near& y) noexcept -> T
{
    using u128 = boost::int128::uint128_t;

    const auto one {fx_one<3>()};
    u128 num {1U};
    u128 den {1U};
    for (unsigned i {2U}; i < x.a; ++i)
    {
        num *= i;
    }
    for (unsigned i {2U}; i < y.a; ++i)
    {
        num *= i;
    }
    for (unsigned i {2U}; i < x.a + y.a; ++i)
    {
        den *= i;
    }

    // The sums from a and from b, A and B in quarters
    fx<3> sa {};
    fx<3> sb {};
    fx<3> sa2 {};
    fx<3> sb2 {};
    fx<3> h2 {};
    for (unsigned i {1U}; i < x.a + y.a; ++i)
    {
        const auto r1 {fx_div_small(one, 4U * i)};
        const auto r2 {fx_div_small(one, static_cast<std::uint64_t>(i) * i)};
        if (i >= x.a)
        {
            sa = fx_add(sa, r1);
            sa2 = fx_add(sa2, r2);
        }
        if (i >= y.a)
        {
            sb = fx_add(sb, r1);
            sb2 = fx_add(sb2, r2);
        }
        h2 = fx_add(h2, r2);
    }
    // zeta(2) = 1 + 2 (zeta(2) - 1) / 2
    const auto psi1 {fx_sub(fx_add(one, fx_shl(fx_top<3>(fx_table::lgamma_zeta[0]), 1U)), h2)};

    // For eps = e 10^-25 and db = d 10^-25, the terms P1 10^-25 and P2 10^-50 of S give g = expm1(S) 10^25 / 4 =
    // P1 / 4 + 10^-25 (P2 / 4 + 2 (P1 / 4)^2)
    const auto e {fx_from_decimal<3>(x.s, x.k)};
    const auto d {fx_from_decimal<3>(y.s, y.k)};
    const auto p1 {fx_signed_add<3>({fx_mul(sa, e), !x.neg}, {fx_mul(sb, d), !y.neg})};
    const auto sq {fx_shr(fx_add(fx_mul(sa2, fx_mul(e, e)), fx_mul(sb2, fx_mul(d, d))), 3U)};
    const auto p2 {fx_signed_add<3>({sq, false}, {fx_shr(fx_mul(fx_mul(e, d), psi1), 2U), x.neg == y.neg})};
    const auto inner {fx_signed_add<3>(p2, {fx_shl(fx_mul(p1.m, p1.m), 1U), false})};
    const auto g {fx_signed_add<3>(p1, {fx_mul(inner.m, fx_from_decimal<3>(1U, 25)), inner.neg})};

    // q0 + r0 / Q = P 10^K / Q with 38 digits in q0
    int k {37 - (num_digits(num) - num_digits(den))};
    const auto scaled {[&](const int n) { return umul256(num * pow10(static_cast<u128>(n - 38 > 0 ? n - 38 : 0)), pow10(static_cast<u128>(n < 38 ? n : 38))); }};
    auto qr {impl::div_mod(scaled(k), den)};
    if (static_cast<u128>(qr.quotient) < pow10(static_cast<u128>(37)))
    {
        ++k;
        qr = impl::div_mod(scaled(k), den);
    }
    const auto q0 {static_cast<u128>(qr.quotient)};
    const auto r0 {static_cast<u128>(qr.remainder)};

    // floor(q0 + r0 / Q + (q0 + r0 / Q) E) for E = expm1(S) = 4 g 10^-25: q0 E = Z / D2 for Z = q0 g 2^40 and
    // D2 = 10^25 / 4 * 2^40, and r0 E D2 = r0 g 2^40
    const auto d2 {(pow10(static_cast<u128>(25)) / 4U) << 40U};
    const auto yq {impl::div_mod(beta_scale40(q0, g.m), d2)};
    const auto yi {static_cast<u128>(yq.quotient)};
    const auto lhs {umul256(static_cast<u128>(yq.remainder), den) + beta_scale40(r0, g.m)};
    const auto base {umul256(r0, d2)};
    u128 digits {};
    if (g.neg)
    {
        digits = q0 - yi - (base < lhs ? 1U : 0U);
    }
    else
    {
        const auto dq {umul256(d2, den)};
        const auto n {lhs + base};
        digits = q0 + yi + (n >= dq ? 1U : 0U) + (n >= dq + dq ? 1U : 0U);
    }
    return fx_round<T>(digits, -k, false);
}

// beta(x, y) = e^st * c in N words: from beta_poly for x + y < 20, else with gamma(z) = e^st(z) / c(z) from gamma_recur
// and gamma(w + 1) / w for w = x + y < 1. One or two words compute it again in one more near a rounding boundary.
template <typename T, int N>
constexpr auto beta_words(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const beta_integer& b) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto w {fx_wide_add(gamma_detail::gamma_wide<N>(sx, ex), gamma_detail::gamma_wide<N>(sy, ey))};
    gamma_exp<N> g {};
    // In one word the fit of gamma_d32 holds at any y, inside the 48 bits that the check below allows
    if (N < 3 && gamma_detail::gamma_whole(w) < 20U)
    {
        BOOST_DECIMAL_IF_CONSTEXPR (p == 16 && N == 2)
        {
            T r {};
            if (beta_try(sx, ex, sy, ey, r))
            {
                return r;
            }
        }
        g = beta_poly<N>(sx, ex, sy, ey, p);
    }
    else
    {
        BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
        {
            // For x + y < 20, a first try in two words with the precise series of exp: its error is below 640 units of
            // the digits
            if (gamma_detail::gamma_whole(w) < 20U)
            {
                const auto g2 {beta_poly<2>(sx, ex, sy, ey, p)};
                const auto r2 {gamma_detail::gamma_exp_digits<T, 2, true>(g2.st, g2.c)};
                if (r2.kind == exp_detail::exp_kind::digits && !fx_near_boundary(r2.digits, 640U))
                {
                    return fx_round<T>(r2.digits, r2.exp10, false);
                }
            }
        }
        BOOST_DECIMAL_IF_CONSTEXPR (p == 16 && N == 2)
        {
            T r {};
            if (beta_try_large(sx, ex, sy, ey, r))
            {
                return r;
            }
        }
        const auto gx {gamma_detail::gamma_positive<N>(sx, ex, p)};
        const auto gy {gamma_detail::gamma_positive<N>(sy, ey, p)};
        gamma_parts<N> gs {};
        if (gamma_detail::gamma_whole(w) >= 1U)
        {
            gs = gamma_detail::gamma_recur(w);
        }
        else
        {
            gs = gamma_detail::gamma_recur(fx_wide_add(w, fx_wide_int<N>(1U)));
            gs.c = gamma_detail::gamma_fmul(gs.c, beta_small_sum<N>(sx, ex, sy, ey, p));
        }
        g.c = gamma_detail::gamma_fmul(gs.c, gamma_detail::gamma_frecip(gamma_detail::gamma_fmul(gx.c, gy.c)));
        g.st = fx_wide_add(fx_wide_add(gx.st, gy.st), gs.st, true);
    }
    const auto r {gamma_detail::gamma_exp_digits<T>(g.st, g.c)};
    if (r.kind != exp_detail::exp_kind::digits)
    {
        return exp_detail::exp10_beyond<T, N>(r.kind, false);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N < 3)
    {
        if (beta_near_boundary(r.digits, p, N == 1 ? 48 : 72))
        {
            return beta_words<T, (N < 3 ? N + 1 : 3)>(sx, ex, sy, ey, b);
        }
    }
    if (b.n != 0U)
    {
        // The value next to the digits, without its trailing zeros
        const int drop {num_digits(r.digits) - p};
        const auto unit {pow10(static_cast<boost::int128::uint128_t>(drop))};
        auto sig {(r.digits + unit / 2U) / unit};
        const auto t {remove_trailing_zeros(sig)};
        if (beta_is_exact(b, t.trimmed_number, r.exp10 + drop + static_cast<int>(t.number_of_removed_zeros)))
        {
            return T {sig, r.exp10 + drop};
        }
    }
    return fx_round<T>(r.digits, r.exp10, false);
}

} // namespace beta_detail

template <typename T>
constexpr auto beta_impl(const T x, const T y) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    // beta is defined for x > 0 and y > 0, and other arguments are a domain error
    if (isnan(x) || isnan(y) || x <= T {0} || y <= T {0})
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (isinf(x) || isinf(y))
    {
        return T {0};
    }
    #else
    if (x <= T {0} || y <= T {0})
    {
        return T {0};
    }
    #endif

    // beta(x, 1) = 1 / x
    if (x == T {1} || y == T {1})
    {
        return T {1} / (x == T {1} ? y : x);
    }
    constexpr int p {std::numeric_limits<T>::digits10};
    int ex {};
    int ey {};
    const auto sx {static_cast<boost::int128::uint128_t>(frexp10(x, &ex))};
    const auto sy {static_cast<boost::int128::uint128_t>(frexp10(y, &ey))};

    BOOST_DECIMAL_IF_CONSTEXPR (p > 16)
    {
        // Next to integers a and b, beta(a, b) is rational and the first order terms are exact decimals
        const auto nx {beta_detail::beta_near_of(sx, ex)};
        const auto ny {beta_detail::beta_near_of(sy, ey)};
        if (nx.ok && ny.ok && nx.a + ny.a <= 30U && (nx.s != 0U || ny.s != 0U))
        {
            return beta_detail::beta_near_integers<T>(nx, ny);
        }
    }

    // A tiny m = min(x, y) next to M = max(x, y), with m < 10^mag_m and M < 10^mag_M
    const bool x_min {x < y};
    const auto s_m {x_min ? sx : sy};
    const int e_m {x_min ? ex : ey};
    const auto s_M {x_min ? sy : sx};
    const int e_M {x_min ? ey : ex};
    const int mag_m {e_m + p};
    const int mag_M {e_M + p};
    if (mag_m <= -20 && mag_m - mag_M <= -20)
    {
        const int bound {mag_M <= 0 ? 2 - mag_M : num_digits(static_cast<unsigned>(3 + 3 * mag_M)) + 1};
        T r {};
        if (beta_detail::beta_tiny<T>(s_m, e_m, x_min ? y < T {1} : x < T {1}, bound, r))
        {
            return r;
        }
    }

    constexpr int words {p <= 16 ? 2 : 3};
    if (ex + p <= 15 && ey + p <= 15)
    {
        // beta(n, z) can be exact for an integer n
        const bool ix {gamma_detail::gamma_is_integer(sx, ex)};
        const bool iy {gamma_detail::gamma_is_integer(sy, ey)};
        const bool use_x {ix && (!iy || x < y)};
        beta_detail::beta_integer b {0U, 0, 0U};
        if ((ix || iy) && (use_x ? x : y) <= T {35})
        {
            const auto sn {use_x ? sx : sy};
            const int en {use_x ? ex : ey};
            const auto tz {remove_trailing_zeros(use_x ? sy : sx)};
            const int ez {(use_x ? ey : ex) + static_cast<int>(tz.number_of_removed_zeros)};
            b.n = static_cast<unsigned>(gamma_detail::gamma_integer_of(sn, en));
            b.a = ez >= 0 ? tz.trimmed_number * pow10(static_cast<boost::int128::uint128_t>(ez)) : tz.trimmed_number;
            b.k = ez >= 0 ? 0 : -ez;
        }
        // One word for decimal32_t below 1000, else two words, or three for decimal128_t
        if (p <= 7 && ex + p <= 3 && ey + p <= 3)
        {
            return beta_detail::beta_words<T, (p <= 7 ? 1 : words)>(sx, ex, sy, ey, b);
        }
        return beta_detail::beta_words<T, words>(sx, ex, sy, ey, b);
    }

    // M >= 10^15, where m >= 10^4 gives beta(m, M) <= beta(10^4, 10^15) < e^-250000
    if (mag_m <= 4)
    {
        return beta_detail::beta_large<T, words>(s_m, e_m, s_M, e_M);
    }
    return T {1U, -2 * std::numeric_limits<T>::max_exponent10};
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto beta(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::beta_impl(static_cast<evaluation_type>(y), static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP
