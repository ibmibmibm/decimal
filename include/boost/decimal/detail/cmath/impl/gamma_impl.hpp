// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_GAMMA_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_GAMMA_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/u256.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/cmath/impl/exp_impl.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace gamma_detail {

template <int N>
using wide = fx_wide<N>;

// The Stirling series in N words converges for y >= X0, with the terms in its table
constexpr auto gamma_x0(const int n) noexcept -> std::uint64_t
{
    return n == 1 ? 10U : n == 2 ? 20U : 30U;
}

// The terms of the Stirling series in w, with the later ones in fewer words. In one word, the terms
// after the first 7 are below 2^-64 for w <= 1/4.
constexpr auto gamma_stirling_sum(const fx<1>& w) noexcept -> fx<1>
{
    const auto c {fx_make_span(fx_table::stirling_x10)};
    std::size_t k {w.w[0] <= (UINT64_C(1) << 60U) ? 7U : c.size()};
    fx<1> s {c[k - 1U]};
    while (--k > 0U)
    {
        s = fx_sub(c[k - 1U], fx_mul(w, s));
    }
    return s;
}

constexpr auto gamma_stirling_sum(const fx<2>& w) noexcept -> fx<2>
{
    return fx_alternating_split<fx_table::stirling_x20_shift>(w, fx_table::stirling_x20_head, fx_table::stirling_x20_tail);
}

constexpr auto gamma_stirling_sum(const fx<3>& w) noexcept -> fx<3>
{
    return fx_alternating_split<fx_table::stirling_x30_mid_shift, fx_table::stirling_x30_tail_shift>(w, fx_table::stirling_x30_head, fx_table::stirling_x30_mid, fx_table::stirling_x30_tail);
}

// The integer part of a >= 0
template <int N>
constexpr auto gamma_whole(const wide<N>& a) noexcept -> std::uint64_t
{
    return (a.w[N] << 2U) | (a.w[N - 1] >> 62U);
}

// v * 2^n2 * 10^n10 with v in [1, 2)
template <int N>
struct gamma_float
{
    fx<N> v;
    int n2;
    int n10;
};

// a with v in [1, 4) as v in [1, 2)
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto gamma_fnorm(gamma_float<N> a) noexcept -> gamma_float<N>
{
    if (a.v.w[N - 1] >= (UINT64_C(1) << 63U))
    {
        a.v = fx_shr(a.v, 1U);
        ++a.n2;
    }
    return a;
}

template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto gamma_fmul(const gamma_float<N>& a, const gamma_float<N>& b) noexcept -> gamma_float<N>
{
    return gamma_fnorm(gamma_float<N> {fx_mul(a.v, b.v), a.n2 + b.n2, a.n10 + b.n10});
}

// 1 / a, with 2 / v in (1, 2]
template <int N>
constexpr auto gamma_frecip(const gamma_float<N>& a) noexcept -> gamma_float<N>
{
    return gamma_fnorm(gamma_float<N> {fx_shl(fx_recip(a.v), 1U), -a.n2 - 1, -a.n10});
}

// sig * 10^e for a sig of nd digits
template <int N>
constexpr auto gamma_fdec(const boost::int128::uint128_t& sig, const int nd, const int e) noexcept -> gamma_float<N>
{
    unsigned s {};
    const auto v {fx_normalize(fx_from_decimal<N>(sig, nd), s)};
    return {v, -static_cast<int>(s), e + nd};
}

// a / 2^s as fixed point, for a wide a >= 0 below 2^(s + 2)
template <int N>
constexpr auto gamma_shr(const wide<N>& a, const unsigned s) noexcept -> fx<N>
{
    fx<N> r {};
    fx_shr_words(fx_make_span(r.w), fx_make_span<N>(a.w), a.w[N], s);
    return r;
}

// a in one word, truncated
constexpr auto gamma_ftop(const gamma_float<2>& a) noexcept -> gamma_float<1>
{
    return {fx_top<1>(a.v), a.n2, a.n10};
}

// a >= 1 as v * 2^n2 with v in [1, 2)
template <int N>
constexpr auto gamma_normalize(const wide<N>& a) noexcept -> gamma_float<N>
{
    const auto s {static_cast<unsigned>(63 - countl_zero(gamma_whole(a)))};
    gamma_float<N> r {{}, static_cast<int>(s), 0};
    if (s == 0U)
    {
        r.v = fx_low<N>(a.w);
    }
    else
    {
        r.v = gamma_shr(a, s);
    }
    return r;
}

// x = sig * 10^e < 10^16 as wide fixed point
template <int N>
constexpr auto gamma_wide(const boost::int128::uint128_t& sig, const int e) noexcept -> wide<N>
{
    if (e >= 0)
    {
        return fx_wide_int<N>(static_cast<std::uint64_t>(sig * pow10(static_cast<boost::int128::uint128_t>(e))));
    }
    const int k {-e};
    const auto whole {static_cast<std::uint64_t>(exp_detail::exp_whole(sig, k))};
    return fx_wide_add(fx_widen(exp_detail::exp_frac<N>(sig, k, whole)), fx_wide_int<N>(whole));
}

// a * 2^s for 0 < s < 64, where the result fits
template <int N>
constexpr auto gamma_shl(wide<N> a, const unsigned s) noexcept -> wide<N>
{
    for (int i {N}; i > 0; --i)
    {
        a.w[i] = (a.w[i] << s) | (a.w[i - 1] >> (64U - s));
    }
    a.w[0] <<= s;
    return a;
}

// 1 / y for y = sig * 10^e with nd digits: 1 / v * 2^s * 10^-(e + nd), with y = u * 10^(e + nd) and u = v * 2^-s in
// [0.1, 1)
template <int N>
constexpr auto gamma_inv(const boost::int128::uint128_t& sig, const int nd, const int e) noexcept -> fx<N>
{
    unsigned s {};
    const auto v {fx_normalize(fx_from_decimal<N>(sig, nd), s)};
    return fx_div(fx_shl(fx_from_decimal<N>(1U, e + nd), s), v);
}

// (y - 1/2) ln(y) - y + ln(2 pi) / 2 for y = f, the head of the Stirling series
template <int N>
constexpr auto gamma_stirling_head(const wide<N>& y, const gamma_float<N>& f) noexcept -> wide<N>
{
    const auto ln_y {log_detail::log_wide_of(f.v, f.n2, 0, log_detail::log_base::e)};
    wide<N> half {};
    half.w[N - 1] = UINT64_C(1) << 61U;

    // (y - 1/2) * ln(y) / 16 * 16, with ln(y) < 37
    const auto r {gamma_shl(fx_wide_mul(fx_wide_add(y, half, true), gamma_shr(ln_y, 4U)), 4U)};
    return fx_wide_add(fx_wide_add(r, y, true), fx_widen(fx_top<N>(fx_table::ln_sqrt_2pi)));
}

// The rest z * s(w) of the Stirling series in M words for y = f, with z = 1 / y and the terms of s alternating in
// sign in w = (X0 z)^2 <= 1
template <int M, int N>
constexpr auto gamma_stirling_tail(const gamma_float<N>& f) noexcept -> fx<M>
{
    const auto z {fx_shr(fx_recip(fx_top<M>(f.v)), static_cast<unsigned>(f.n2))};
    const auto xz {fx_mul_small(z, gamma_x0(M))};
    return fx_mul(z, gamma_stirling_sum(fx_mul(xz, xz)));
}

// ln gamma(y) for X0 <= y < 10^16 by the Stirling series
template <int N>
constexpr auto gamma_stirling(const wide<N>& y) noexcept -> wide<N>
{
    const auto f {gamma_normalize(y)};
    return fx_wide_add(gamma_stirling_head(y, f), fx_widen(gamma_stirling_tail<N>(f)));
}

// ln gamma(y) for 10 <= y < 10^16 with the head of the Stirling series in two words and its tail in one
constexpr auto gamma_stirling_one(const wide<2>& y) noexcept -> wide<2>
{
    const auto f {gamma_normalize(y)};
    return fx_wide_add(gamma_stirling_head(y, f), fx_widen(fx_extend<2>(gamma_stirling_tail<1>(f))));
}

// ln|gamma(x)| = st - ln(c), or ln(c) - st with neg_st, and the sign neg of gamma(x)
template <int N>
struct gamma_parts
{
    wide<N> st;
    gamma_float<N> c;
    bool neg_st;
    bool neg;
};

// ln gamma(y) for 1 <= y < 10^16. Below X0, ln gamma(y) = ln gamma(y + n) - ln(y (y + 1) ... (y + n - 1)) with y + n >= X0.
template <int N>
constexpr auto gamma_shift(wide<N> y) noexcept -> gamma_parts<N>
{
    gamma_parts<N> r {{}, {fx_one<N>(), 0, 0}, false, false};
    const std::uint64_t whole {gamma_whole(y)};
    if (whole < gamma_x0(N))
    {
        // With v = 2y + n - 1, the pairs (y + i)(y + n - 1 - i) are (v^2 - k^2) / 4 for k = n - 1 - 2i, and an odd n
        // has the middle factor v / 2: one product for two factors
        const std::uint64_t n {gamma_x0(N) - whole};
        auto vf {gamma_normalize(fx_wide_add(fx_wide_add(y, y), fx_wide_int<N>(n - 1U)))};
        const auto v2 {gamma_shl(fx_widen(fx_mul(vf.v, vf.v)), 2U * static_cast<unsigned>(vf.n2))};
        // Two pairs (v^2 - a)(v^2 - b) = v^4 - (a + b) v^2 + ab: one product for four factors, with -(a + b) v^2 as
        // (a + b) times -v^2 modulo 2^(64 (N + 1))
        const auto v2f {gamma_normalize(v2)};
        const auto v4 {gamma_shl(fx_widen(fx_mul(v2f.v, v2f.v)), 2U * static_cast<unsigned>(v2f.n2))};
        const auto neg_v2 {fx_wide_negate(v2)};
        std::uint64_t k {n - 1U};
        for (; k > 2U && k < n; k -= 4U)
        {
            const std::uint64_t a {k * k};
            const std::uint64_t b {(k - 2U) * (k - 2U)};
            auto d {fx_wide_add(v4, fx_wide_int<N>(a * b))};
            fx_mul_add(fx_make_span(d.w), fx_make_span(neg_v2.w), a + b);
            auto f {gamma_normalize(d)};
            f.n2 -= 4;
            r.c = gamma_fmul(r.c, f);
        }
        if (k > 0U && k < n)
        {
            auto d {v2};
            const std::uint64_t low {(k * k) << 62U};
            d.w[N] -= ((k * k) >> 2U) + (d.w[N - 1] < low ? 1U : 0U);
            d.w[N - 1] -= low;
            auto f {gamma_normalize(d)};
            f.n2 -= 2;
            r.c = gamma_fmul(r.c, f);
        }
        if (n % 2U == 1U)
        {
            --vf.n2;
            r.c = gamma_fmul(r.c, vf);
        }
        y = fx_wide_add(y, fx_wide_int<N>(n));
    }
    r.st = gamma_stirling(y);
    return r;
}

// ln gamma(x) for 0 < x = sig * 10^e < 10^16 with p digits in sig, with ln gamma(x) = ln gamma(x + 1) - ln(x) below 1
template <int N>
constexpr auto gamma_positive(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> gamma_parts<N>
{
    const auto y {gamma_wide<N>(sig, e)};
    if (e + p > 0)
    {
        return gamma_shift(y);
    }
    auto r {gamma_shift(fx_wide_add(y, fx_wide_int<N>(1U)))};
    r.c = gamma_fmul(r.c, gamma_fdec<N>(sig, p, e));
    return r;
}

// F(z) = (1 - euler) + sum (-1)^k (zeta(k) - 1) / k * z^(k - 1) for k >= 2 and |z| < 1/16, with ln gamma(2 + z) = z F(z).
// The terms of the sum alternate in sign for z > 0, and are all negative for z < 0.
template <int N>
constexpr auto gamma_near_two(const fx<N>& az, const bool neg_z) noexcept -> fx<N>
{
    // The first term left out, at k = 16, 29 or 42, is below 2^-79, 2^-145 and 2^-211
    const auto zeta {fx_make_span<N == 1 ? 14 : N == 2 ? 27 : 40>(fx_table::lgamma_zeta)};
    auto s {fx_top<N>(zeta[zeta.size() - 1U])};
    for (std::size_t k {zeta.size() - 1U}; k-- > 0U;)
    {
        s = neg_z ? fx_add(fx_top<N>(zeta[k]), fx_mul(az, s)) : fx_sub(fx_top<N>(zeta[k]), fx_mul(az, s));
    }
    const auto c {fx_top<N>(fx_table::one_minus_euler)};
    return neg_z ? fx_sub(c, fx_mul(az, s)) : fx_add(c, fx_mul(az, s));
}

// ln gamma(2 + z), or ln gamma(1 + z) with one and |z| < 1/256, for z = dz * 10^e != 0 with the sign neg_z and |z| < 1/16.
// With ln(1 + z) = z * l(z) as in log_near, ln gamma(1 + z) = z F(z) - ln(1 + z) = -z (l(z) - F(z)).
template <int N>
constexpr auto gamma_near(const boost::int128::uint128_t& dz, const int e, const bool neg_z, const bool one) noexcept -> fx_wide_digits_of<N>
{
    using digits_type = fx_digits_t<N>;

    const auto az {fx_from_decimal<N>(dz, -e)};
    auto m {gamma_near_two(az, neg_z)};
    bool neg {neg_z};
    if (one)
    {
        const auto l {log_detail::log_ratio(az, neg_z)};
        m = fx_sub(l, m);
        neg = !neg_z;
    }

    // m < 1, thus m * dz scaled to len digits has len - 1 or len digits
    constexpr int len {N == 1 ? 19 : 38};
    const int nd {num_digits(dz)};
    return {static_cast<digits_type>(fx_scale(dz * pow10(static_cast<boost::int128::uint128_t>(len - nd)), m)), e - (len - nd), neg};
}

// Whether digits * 10^exp10 from two words, with an error below 2^-112 * 10^(digits of the integer part), can round
// another way than the exact value
template <typename digits_type>
constexpr auto gamma_near_rounding(const digits_type& digits, const int exp10, const int p) noexcept -> bool
{
    constexpr int scale10 {fx_scale10(2)};
    constexpr int kept {33};
    const int len {digits >= pow10(static_cast<digits_type>(scale10)) ? scale10 + 1 : scale10};
    const int err {-exp10 - kept + (exp10 + len > 1 ? exp10 + len : 0)};
    if (err >= len - p - 1)
    {
        return true;
    }
    const auto m {pow10(static_cast<digits_type>(len - p - 1)) * 5U};
    const auto r {digits % m};
    const auto w {pow10(static_cast<digits_type>(err > 0 ? err + 1 : 1))};
    return r < w || r > m - w;
}

// a + z s / 2^30, or a - z s / 2^30 for neg: one step of the series of gamma_root
template <int N>
constexpr auto gamma_root_step(const fx<N>& a, const fx<N>& z, const fx<N>& s, const bool neg) noexcept -> fx<N>
{
    const auto m {fx_shr(fx_mul(z, s), 30U)};
    return neg ? fx_sub(a, m) : fx_add(a, m);
}

// ln|gamma(x)| for decimal128_t x = -sig * 10^e (e = -32 above ten, else -33) next to a root r, rounded once in out. For
// ||x| - |r|| = (n + f) 10^e and v = 4 (n + f) 10^-k10 < 2^-30, it is v S(v), or -v S(v) with S alternating for neg.
template <typename T>
constexpr auto gamma_root(const boost::int128::uint128_t& sig, const bool above_ten, T& out) noexcept -> bool
{
    // The roots below 10 come first, each decade by |r|
    constexpr int count {static_cast<int>(sizeof(fx_table::lgamma_root_sig) / sizeof(fx_table::lgamma_root_sig[0]))};
    constexpr int split {16};
    const int first {above_ten ? split : 0};
    const int last {above_ten ? count : split};
    int i {first};
    int j {last};
    while (i < j)
    {
        const int mid {i + (j - i) / 2};
        i = fx_table::lgamma_root_sig[mid] < sig ? mid + 1 : i;
        j = fx_table::lgamma_root_sig[mid] < sig ? j : mid;
    }
    for (i = i > first ? i - 1 : i; i < last && i <= j; ++i)
    {
        const auto c {fx_table::lgamma_root_sig[i]};
        const bool above {sig > c};
        const auto d {above ? sig - c : c - sig};
        const int k10 {fx_table::lgamma_root_k10[i]};
        if (d >= (pow10(static_cast<boost::int128::uint128_t>(k10)) >> 32U))
        {
            continue;
        }

        // n = d - 1 and f = 1 - g above c, else n = d and f = g, and z = v * 2^30
        const auto a {fx_make_span<3>(fx_table::lgamma_root, i * 3)};
        const auto n {above ? d - 1U : d};
        const auto f {above ? fx_sub(fx_one<3>(), a[0]) : a[0]};
        const auto z {fx_shl(fx_add(fx_from_decimal<3>(n, k10), fx_mul(f, fx_from_decimal<3>(1U, k10))), 32U)};
        const bool neg {above == (i % 2 == 0)};
        const auto t {fx_make_span<3>(fx_table::lgamma_root_tail, i * 3)};
        const auto m {fx_make_span<2>(fx_table::lgamma_root_mid, i * 2)};
        const auto z1 {fx_top<1>(z)};
        auto s1 {gamma_root_step(fx<1> {{t[1]}}, z1, fx<1> {{t[2]}}, neg)};
        s1 = gamma_root_step(fx<1> {{t[0]}}, z1, s1, neg);
        const auto z2 {fx_top<2>(z)};
        const auto s2 {gamma_root_step(m[0], z2, gamma_root_step(m[1], z2, fx_extend<2>(s1), neg), neg)};
        auto s {gamma_root_step(a[1], z, gamma_root_step(a[2], z, fx_extend<3>(s2), neg), neg)};

        // 0.24 < 4 S < 2.4 and (n + f) 10^(38 - nd) is in [8 * 10^35, 10^38]: the digits of v S fit
        s = fx_shl(s, 2U);
        const int nd {n == 0U ? 0 : num_digits(n)};
        const auto scale {pow10(static_cast<boost::int128::uint128_t>(38 - nd))};
        std::uint64_t tail {};
        const auto low {fx_scale_frac(scale, f, tail)};
        out = fx_round<T>(fx_scale(scale * n + low, tail, s), nd - 38 - k10, neg);
        return true;
    }
    return false;
}

// sin(pi f) / (pi f) for f = d * 10^-k <= 1/2 by the Taylor series in w = (2f)^2 <= 1, whose terms alternate in sign. One
// word takes the top words of the terms of two words, as the terms after them are below 2^-67.
constexpr auto gamma_sinc(const boost::int128::uint128_t& d, const int k, std::integral_constant<int, 1>) noexcept -> fx<1>
{
    const auto t {fx_from_decimal<1>(d * 2U, k)};
    const auto w {fx_mul(t, t)};
    const auto c {fx_make_span(fx_table::sinpi_x2_head)};
    auto s {fx_top<1>(c[c.size() - 1U])};
    for (std::size_t i {c.size() - 1U}; i-- > 0U;)
    {
        s = fx_sub(fx_top<1>(c[i]), fx_mul(w, s));
    }
    return s;
}

// The same in two or three words, with the later terms in fewer words
constexpr auto gamma_sinc(const boost::int128::uint128_t& d, const int k, std::integral_constant<int, 2>) noexcept -> fx<2>
{
    const auto t {fx_from_decimal<2>(d * 2U, k)};
    return fx_alternating_split<fx_table::sinpi_x2_shift>(fx_mul(t, t), fx_table::sinpi_x2_head, fx_table::sinpi_x2_tail);
}

constexpr auto gamma_sinc(const boost::int128::uint128_t& d, const int k, std::integral_constant<int, 3>) noexcept -> fx<3>
{
    const auto t {fx_from_decimal<3>(d * 2U, k)};
    return fx_alternating_split<fx_table::sinpi_x3_mid_shift, fx_table::sinpi_x3_tail_shift>(fx_mul(t, t), fx_table::sinpi_x3_head, fx_table::sinpi_x3_mid, fx_table::sinpi_x3_tail);
}

// |sin(pi x)| / pi, and the sign of gamma(x)
template <int N>
struct gamma_sine
{
    gamma_float<N> s;
    bool neg;
};

// For x = -a < 0 that is not an integer, a = sig * 10^e: f = x - round(x) = -(a - n) is exact, sin(pi x) =
// -(-1)^n sin(pi (a - n)), and |sin(pi x)| / pi = |f| sin(pi f) / (pi f). The sign of gamma(x) is that of sin(pi x).
template <int N>
constexpr auto gamma_reflect(const boost::int128::uint128_t& sig, const int e) noexcept -> gamma_sine<N>
{
    // |a - n| = d * 10^-k, for the integer n nearest to a
    const int k {-e};
    auto d {sig};
    bool odd {false};
    bool below {false};
    if (k < 39)
    {
        const auto ten {pow10(static_cast<boost::int128::uint128_t>(k))};
        const auto n {sig / ten};
        d = sig - n * ten;
        odd = (n.low & 1U) != 0U;
        if (d * 2U > ten)
        {
            d = ten - d;
            odd = !odd;
            below = true;
        }
    }
    const bool neg {odd == below};

    const auto s {gamma_sinc(d, k, std::integral_constant<int, N> {})};
    return {gamma_fmul(gamma_fdec<N>(d, num_digits(d), e), gamma_float<N> {fx_shl(fx_clamp_below_one(s), 1U), -1, 0}), neg};
}

// ln|gamma(x)| = ln(pi) - ln|sin(pi x)| - ln gamma(1 + a) for x = -a < 0 as in gamma_reflect, or gamma_positive for x > 0
template <int N>
constexpr auto gamma_split(const boost::int128::uint128_t& sig, const int e, const bool neg_x, const int p) noexcept -> gamma_parts<N>
{
    if (!neg_x)
    {
        return gamma_positive<N>(sig, e, p);
    }
    const auto g {gamma_reflect<N>(sig, e)};
    auto r {gamma_shift(fx_wide_add(gamma_wide<N>(sig, e), fx_wide_int<N>(1U)))};
    r.c = gamma_fmul(r.c, gamma_frecip(g.s));
    r.neg_st = true;
    r.neg = g.neg;
    return r;
}

// ln gamma(y), negated with neg, for y = sig * 10^e >= 10^16, with extra added to the part below y: ln gamma(y) =
// y * (ln(y) - 1 + (ln(2 pi) / 2 - ln(y) / 2 + extra) / y + 1 / (12 y^2)), where the next term is below 2^-200
template <typename T, int N>
constexpr auto gamma_large(const boost::int128::uint128_t& sig, const int e, const wide<N>& extra, const bool neg) noexcept -> T
{
    const int nd {num_digits(sig)};
    const auto ln_y {log_detail::log_wide_of<N>(sig, nd, e, log_detail::log_base::e)};

    const auto inv_y {gamma_inv<N>(sig, nd, e)};

    // r = ln(2 pi) / 2 - ln(y) / 2 + extra, and q = ln(y) - 1 + r / y + 1 / (12 y^2)
    wide<N> half_ln {};
    fx_shr_words(fx_make_span(half_ln.w), fx_make_span(ln_y.w), 0U, 1U);
    const auto r {fx_wide_add(fx_wide_add(fx_widen(fx_top<N>(fx_table::ln_sqrt_2pi)), half_ln, true), extra)};
    bool neg_r {};
    const auto abs_r {fx_wide_abs(r, neg_r)};
    auto q {fx_wide_add(ln_y, fx_wide_int<N>(1U), true)};
    q = fx_wide_add(q, fx_wide_mul(abs_r, inv_y), neg_r);
    q = fx_wide_add(q, fx_widen(fx_div_small(fx_mul(inv_y, inv_y), 12U)));

    // y * q = sig * 10^e * q with sig scaled so that the product has 37 to 39 digits
    const auto q_whole {gamma_whole(q)};
    const auto q_frac {fx_frac<N>(q.w)};
    const int extra_digits {38 - num_digits(q_whole) - nd};
    const int scale {extra_digits > 0 ? extra_digits : 0};
    if (e - scale + 39 > std::numeric_limits<T>::max_exponent10 + std::numeric_limits<T>::digits10 + 1)
    {
        return T {1U, 2 * std::numeric_limits<T>::max_exponent10, neg};
    }
    const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(scale))};
    return fx_round<T>(m * q_whole + fx_scale(m, q_frac), e - scale, neg);
}

// Whether x = sig * 10^e is an integer
constexpr auto gamma_is_integer(const boost::int128::uint128_t& sig, const int e) noexcept -> bool
{
    if (e < 0 && e > -20 && sig.high == 0U)
    {
        return sig.low % pow10(static_cast<std::uint64_t>(-e)) == 0U;
    }
    return e >= 0 || (e > -39 && sig % pow10(static_cast<boost::int128::uint128_t>(-e)) == 0U);
}

// The integer sig * 10^e, for gamma_is_integer(sig, e) and an integer that fits
constexpr auto gamma_integer_of(const boost::int128::uint128_t& sig, const int e) noexcept -> boost::int128::uint128_t
{
    return e >= 0 ? sig * pow10(static_cast<boost::int128::uint128_t>(e)) : sig / pow10(static_cast<boost::int128::uint128_t>(-e));
}

// gamma(n) = (n - 1)! for 1 <= n <= 35
constexpr auto gamma_at(const std::uint64_t n) noexcept -> boost::int128::uint128_t
{
    boost::int128::uint128_t f {1U};
    for (std::uint64_t i {2U}; i < n; ++i)
    {
        f *= i;
    }
    return f;
}

// a * b / 2^60 for signed a and b, truncated toward zero
constexpr auto gamma_smul(const std::int64_t a, const std::int64_t b) noexcept -> std::int64_t
{
    const auto ua {static_cast<std::uint64_t>(a < 0 ? -a : a)};
    const auto ub {static_cast<std::uint64_t>(b < 0 ? -b : b)};
    const auto m {static_cast<std::int64_t>(((static_cast<boost::int128::uint128_t>(ua) * ub) >> 60U).low)};
    return (a < 0) != (b < 0) ? -m : m;
}

// gamma(2.5 + w) for |w| <= 1/2 with the scale 2^60, from gamma_d32 with the even and the odd terms as two chains in w^2
constexpr auto gamma_poly(const std::int64_t w) noexcept -> std::int64_t
{
    const auto c {fx_make_span(fx_table::gamma_d32)};
    const auto v {gamma_smul(w, w)};
    std::size_t k {c.size() - 1U};
    auto s0 {c[k]};
    auto s1 {c[k - 1U]};
    while (k >= 3U)
    {
        k -= 2U;
        s0 = c[k] + gamma_smul(v, s0);
        s1 = c[k - 1U] + gamma_smul(v, s1);
    }
    return s1 + gamma_smul(w, s0);
}

// a * b / 2^62 for signed a and b, rounded, from the product in two's complement
constexpr auto gamma_rmul(const std::int64_t a, const std::int64_t b) noexcept -> std::int64_t
{
    const auto ua {static_cast<std::uint64_t>(a)};
    const auto ub {static_cast<std::uint64_t>(b)};
    auto m {static_cast<boost::int128::uint128_t>(ua) * ub};
    m.high -= (ub & (0U - (ua >> 63U))) + (ua & (0U - (ub >> 63U)));
    return static_cast<std::int64_t>(((m + (UINT64_C(1) << 61U)) >> 62U).low);
}

// gamma(2.5 + w) for |w| <= 1/2 with the scale 2^62 from gamma_d64 by Horner's rule. The sum can reach 2, so it is
// unsigned.
constexpr auto gamma_poly_d64(const std::int64_t w) noexcept -> std::uint64_t
{
    const auto c {fx_make_span(fx_table::gamma_d64)};
    std::size_t k {c.size() - 1U};
    auto s {c[k]};
    while (--k > 0U)
    {
        s = c[k] + gamma_rmul(w, s);
    }
    return static_cast<std::uint64_t>(c[0]) + static_cast<std::uint64_t>(gamma_rmul(w, s));
}

// The terms of lgamma_d64 per interval in two words, and in one
constexpr int gamma_lpoly_head {12};
constexpr int gamma_lpoly_tail {13};

// Whether the tail of each interval sums to below 2^-62, so that its partial sums fit in one word scaled by 2^64
constexpr auto gamma_lpoly_fits() noexcept -> bool
{
    for (int i {}; i < 8; ++i)
    {
        boost::int128::uint128_t sum {};
        for (int k {}; k < gamma_lpoly_tail; ++k)
        {
            sum += fx_table::lgamma_d64_tail[i * gamma_lpoly_tail + k] + UINT64_C(1);
        }
        if (sum.high != 0U)
        {
            return false;
        }
    }
    return true;
}

static_assert(gamma_lpoly_fits(), "the tails fit in one word");

// ln gamma(2 + u) for 0 <= u < 1 from lgamma_d64 at the center y0 of the eighth of u, with t = 16 |u + 2 - y0|. The
// terms from k = 2 on alternate in sign above y0 and are all positive below; the tail runs in one word.
constexpr auto gamma_lpoly(const fx<2>& u) noexcept -> fx<2>
{
    const auto i {static_cast<int>(u.w[1] >> 59U)};
    const fx<2> y0 {{0U, static_cast<std::uint64_t>(2 * i + 1) << 58U}};
    const bool below {fx_less(u, y0)};
    const auto t {fx_shl(below ? fx_sub(y0, u) : fx_sub(u, y0), 4U)};
    const auto a {fx_make_span<gamma_lpoly_head>(fx_table::lgamma_d64_head, i * gamma_lpoly_head)};
    const auto b {fx_make_span<gamma_lpoly_tail>(fx_table::lgamma_d64_tail, i * gamma_lpoly_tail)};

    const fx<1> t1 {{t.w[1]}};
    fx<1> s1 {{b[b.size() - 1U]}};
    for (std::size_t k {b.size() - 1U}; k-- > 0U;)
    {
        const auto m {fx_mul(t1, s1)};
        s1 = below ? fx_add(fx<1> {{b[k]}}, m) : fx_sub(fx<1> {{b[k]}}, m);
    }
    fx<2> s {{s1.w[0], 0U}};
    for (std::size_t k {a.size()}; k-- > 0U;)
    {
        const auto m {fx_mul(t, s)};
        s = below == (k > 1U) ? fx_add(a[k], m) : fx_sub(a[k], m);
    }
    // ln gamma(2 + u) >= 0, but next to u = 0 the error can take the sum below 0, which is near 4 as an unsigned value
    return (s.w[1] >> 63U) != 0U ? fx<2> {} : s;
}

// gamma(u + 2) = c * e^st for the fraction u of y: c from gamma_poly or gamma_poly_d64 for P = 7 or 16 digits and st = 0
// in one word, else c = 1 and st from gamma_lpoly
template <int N>
struct gamma_start
{
    gamma_float<N> c;
    fx<N> st;
};

template <int P>
constexpr auto gamma_poly_start(const wide<1>& z, std::integral_constant<int, P>) noexcept -> gamma_start<1>
{
    // gamma(u + 2) is in [1, 2], and the error of the fit can take it just below 1 or to 2
    gamma_float<1> c {{}, 0, 0};
    c.v.w[0] = P == 7 ? static_cast<std::uint64_t>(gamma_poly(static_cast<std::int64_t>(z.w[0] >> 2U) - (INT64_C(1) << 59U))) << 2U
                      : gamma_poly_d64(static_cast<std::int64_t>(z.w[0]) - (INT64_C(1) << 61U));
    if (c.v.w[0] < (UINT64_C(1) << 62U))
    {
        c.v = fx_shl(c.v, 1U);
        c.n2 = -1;
    }
    return {gamma_fnorm(c), {}};
}

constexpr auto gamma_poly_start(const wide<2>& z, std::integral_constant<int, 16>) noexcept -> gamma_start<2>
{
    return {{fx_one<2>(), 0, 0}, gamma_lpoly(fx_low<2>(z.w))};
}

// gamma(y) = c * e^st for 0 < y < 20 is (y - 1) ... (u + 2) gamma(u + 2) for the fraction u of y, else gamma(u + 2) / y
// below 2, or gamma(u + 2) / (y (y + 1)) below 1 with y as the float small
template <int N, int P = (N == 1 ? 7 : 16)>
constexpr auto gamma_poly_up(const wide<N>& y, const gamma_float<N>& small) noexcept -> gamma_start<N>
{
    const auto n {gamma_whole(y)};
    auto z {y};
    z.w[N] = 0U;
    z.w[N - 1] &= (UINT64_C(1) << 62U) - 1U;
    const auto start {gamma_poly_start(z, std::integral_constant<int, P> {})};
    auto c {start.c};
    if (n >= 2U)
    {
        z = fx_wide_add(z, fx_wide_int<N>(2U));
        for (std::uint64_t k {2U}; k < n; ++k)
        {
            c = gamma_fmul(c, gamma_normalize(z));
            z = fx_wide_add(z, fx_wide_int<N>(1U));
        }
    }
    else
    {
        const auto f {gamma_normalize(n == 1U ? y : fx_wide_add(y, fx_wide_int<N>(1U)))};
        c = gamma_fmul(c, gamma_frecip(n == 1U ? f : gamma_fmul(small, f)));
    }
    return {c, start.st};
}

// The error of c from gamma_poly_up<1, 16> in units of 2^-62, for the integer part n of y: 2 for each product, 3 from y
// in all of them, and 3 from the fit, or 9 with the reciprocals below 2
constexpr auto gamma_poly_error(const std::uint64_t n) noexcept -> std::uint64_t
{
    return n < 4U ? 9U : 2U * n + 2U;
}

// c from gamma_poly_up<1, 16> for y = x = sig * 10^e with 0 < x < 20
constexpr auto gamma_poly_one(const wide<1>& y, const boost::int128::uint128_t& sig, const int e) noexcept -> gamma_float<1>
{
    return gamma_poly_up<1, 16>(y, gamma_whole(y) == 0U ? gamma_ftop(gamma_fdec<2>(sig, 16, e)) : gamma_float<1> {}).c;
}

// ln|gamma(x)| for 0 < |x| = sig * 10^e < 20 for decimal64_t in one word, rounded once in out where it is at least 0.1
// and not next to a rounding boundary: the log adds 10 units of 2^-62 to the error of c, and the reflection 8 for x < 0
template <typename T>
constexpr auto lgamma_try(const boost::int128::uint128_t& sig, const int e, const bool neg_x, T& out) noexcept -> bool
{
    auto y {gamma_wide<1>(sig, e)};
    gamma_float<1> c {};
    if (neg_x)
    {
        // ln|gamma(x)| = -ln(|sin(pi x)| / pi * gamma(1 - x)), with the sine in two words as in tgamma_try
        y = fx_wide_add(y, fx_wide_int<1>(1U));
        c = gamma_fmul(gamma_ftop(gamma_reflect<2>(sig, e).s), gamma_poly_one(y, sig, e));
    }
    else
    {
        c = gamma_poly_one(y, sig, e);
    }
    // n2 ln(2) from two words of ln(2)
    const auto l2 {log_detail::log_mul_small<2>(fx_table::ln2, static_cast<std::uint64_t>(c.n2 < 0 ? -c.n2 : c.n2))};
    const fx_wide<1> t {{l2.w[1], l2.w[2]}};
    const auto r {fx_wide_add(log_detail::log_wide_of(c.v, 0, c.n10, log_detail::log_base::e), t, c.n2 < 0)};
    // The result for x < 0 can be near 0, so keep |r| >= 0.1
    const auto a {(r.w[1] >> 63U) != 0U ? fx_wide_negate(r) : r};
    if (a.w[1] == 0U && a.w[0] < UINT64_C(0x0666666666666666))
    {
        return false;
    }
    const auto d {fx_wide_digits(r, neg_x)};
    const auto err {gamma_poly_error(gamma_whole(y)) + (neg_x ? 18U : 10U)};
    const auto width {(static_cast<boost::int128::uint128_t>(pow10(static_cast<std::uint64_t>(-d.exp10))) * err) >> 62U};
    return fx_try_round(d.digits, d.exp10, d.neg, static_cast<std::uint64_t>(width) + 2U, out);
}

// The parts of gamma(x) for 0 < |x| = sig * 10^e < 12 (20 in two words): gamma(y) for y = x, or 1 - x by reflection,
// from gamma_poly_up
template <int N>
constexpr auto gamma_poly_parts(const boost::int128::uint128_t& sig, const int e, const bool neg_x, const int p) noexcept -> gamma_parts<N>
{
    auto y {gamma_wide<N>(sig, e)};
    gamma_sine<N> g {};
    if (neg_x)
    {
        g = gamma_reflect<N>(sig, e);
        y = fx_wide_add(y, fx_wide_int<N>(1U));
    }
    const auto start {gamma_poly_up(y, gamma_whole(y) == 0U ? gamma_fdec<N>(sig, p, e) : gamma_float<N> {})};
    auto c {start.c};
    const auto st {fx_widen(start.st)};
    if (neg_x)
    {
        c = gamma_frecip(gamma_fmul(g.s, c));
    }
    return {N == 1 || neg_x ? st : fx_wide_negate(st), c, true, g.neg};
}

// The parts of gamma(x): from gamma_poly_parts for |x| < 12 in one word and |x| < 20 in two, else from gamma_split
template <int N>
constexpr auto gamma_parts_of(const boost::int128::uint128_t& sig, const int e, const bool neg_x, const int p) noexcept
    -> std::enable_if_t<(N < 3), gamma_parts<N>>
{
    constexpr std::uint64_t bound {N == 1 ? 12U : 20U};
    const bool poly {e + p < 2 || (e + p == 2 && sig < pow10(static_cast<boost::int128::uint128_t>(p - 2)) * bound)};
    return poly ? gamma_poly_parts<N>(sig, e, neg_x, p) : gamma_split<N>(sig, e, neg_x, p);
}

template <int N>
constexpr auto gamma_parts_of(const boost::int128::uint128_t& sig, const int e, const bool neg_x, const int p) noexcept
    -> std::enable_if_t<(N == 3), gamma_parts<N>>
{
    return gamma_split<N>(sig, e, neg_x, p);
}

// ln|gamma(x)| from its parts, signed
template <int N>
constexpr auto lgamma_of(const gamma_parts<N>& g) noexcept -> wide<N>
{
    auto r {g.st};
    if (g.c.n2 != 0 || g.c.n10 != 0 || fx_less(fx_one<N>(), g.c.v))
    {
        r = fx_wide_add(r, log_detail::log_wide_of(g.c.v, g.c.n2, g.c.n10, log_detail::log_base::e), true);
    }
    return g.neg_st ? fx_wide_negate(r) : r;
}

// ln|gamma(x)| for x = sig * 10^e < 10^16 in N words, rounded once in the current mode. Two words compute it again in
// three where the error can change the rounding, and one word needs no more: a check of all decimal32_t values shows it.
template <typename T, int N>
constexpr auto lgamma_words(const boost::int128::uint128_t& sig, const int e, const bool neg_x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    BOOST_DECIMAL_IF_CONSTEXPR (p == 16 && N == 2)
    {
        T r {};
        if (!decimal_at_least(sig, e, p, 20U, 2) && lgamma_try(sig, e, neg_x, r))
        {
            return r;
        }
        // For |x| >= 20, a first try in one word by the Stirling series: for y = |x| or 1 - x below w < 2^n, its
        // error is below (n + 28) w + 64 units of 2^-62, and 128 more for x < 0
        if (decimal_at_least(sig, e, p, 20U, 2))
        {
            const auto w {static_cast<std::uint64_t>(e < 0 ? exp_detail::exp_whole(sig, -e) : sig) + (neg_x ? 2U : 1U)};
            const auto err {w * static_cast<std::uint64_t>(92 - countl_zero(w)) + (neg_x ? 192U : 64U)};
            const auto d1 {fx_wide_digits(lgamma_of(gamma_parts_of<1>(sig, e, neg_x, p)), false)};
            const auto width {(static_cast<boost::int128::uint128_t>(pow10(static_cast<std::uint64_t>(-d1.exp10))) * err) >> 62U};
            if (fx_try_round(d1.digits, d1.exp10, d1.neg, static_cast<std::uint64_t>(width) + 2U, r))
            {
                return r;
            }
        }
    }
    fx_wide_digits_of<N> d {};
    bool near {false};
    if (!neg_x && (e == -p || e == 1 - p))
    {
        // z = x - 1 or x - 2 is exact, for |z| < 1/256 near 1, where log_near holds, and |z| < 1/16 near 2
        const auto one {pow10(static_cast<boost::int128::uint128_t>(p - 1))};
        const bool two {e == 1 - p && sig * 2U >= one * 3U};
        const auto c {e == -p ? one * 10U : two ? one * 2U : one};
        const bool neg_z {sig < c};
        const auto dz {neg_z ? c - sig : sig - c};
        if (dz == 0U)
        {
            return T {0};
        }
        if (dz * (two ? 16U : 256U) < pow10(static_cast<boost::int128::uint128_t>(-e)))
        {
            d = gamma_near<N>(dz, e, neg_z, !two);
            near = true;
        }
    }
    if (!near)
    {
        BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
        {
            // A first try in two words: its error is below 2^8 units of 2^-126 for 0 < x < 20, 2^9 for -20 < x < 0, and
            // 2^(m + 7) for 20 <= |x| < 10^k < 2^m with m = ceil(10k / 3), by the Stirling series
            T r {};
            const unsigned bits {decimal_at_least(sig, e, p, 20U, 2) ? static_cast<unsigned>((10 * (e + p) + 2) / 3 + 7) : neg_x ? 9U : 8U};
            if (fx_wide_try_round(lgamma_of(gamma_parts_of<2>(sig, e, neg_x, p)), false, r, bits))
            {
                return r;
            }
        }
        d = fx_wide_digits(lgamma_of(gamma_parts_of<N>(sig, e, neg_x, p)), false);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N == 2)
    {
        if (gamma_near_rounding(d.digits, d.exp10, p))
        {
            return lgamma_words<T, 3>(sig, e, neg_x);
        }
    }
    return fx_round<T>(d.digits, d.exp10, d.neg);
}

// ln|gamma(x)| for x = sig * 10^e (negative with neg_x) not 0 or a negative integer, rounded once in the current mode.
// For |x| >= 10^16, ln gamma(|x|) = |x| * q(|x|).
template <typename T>
constexpr auto lgamma_fixed(const boost::int128::uint128_t& sig, const int e, const bool neg_x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int large_words {p <= 16 ? 2 : 3};

    if (e + p > 16)
    {
        if (!neg_x)
        {
            return gamma_large<T, large_words>(sig, e, wide<large_words> {}, false);
        }
        const auto c {gamma_fmul(gamma_fdec<large_words>(sig, p, e), gamma_reflect<large_words>(sig, e).s)};
        return gamma_large<T, large_words>(sig, e, log_detail::log_wide_of(c.v, c.n2, c.n10, log_detail::log_base::e), true);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // Three words lose the digits of the result next to the roots below -2
        T r {};
        if (neg_x && (e == 1 - p || e == 2 - p) && gamma_root(sig, e == 2 - p, r))
        {
            return r;
        }
    }
    return lgamma_words<T, words>(sig, e, neg_x);
}

// 10^t * v, or 10^-t * v with neg, for t >= 0 and 1 <= v < 2 as in exp10_wide: 10^f * v = 2^a * m * v in [1, 20),
// with the wide series of exp for Wide, and 10^f * v in the top M words
template <int N, bool Wide = false, int M = N>
constexpr auto gamma_exp10(const wide<N>& t, const bool neg, const std::uint64_t bound, const fx<N>& v) noexcept -> exp_detail::exp_digits<M>
{
    constexpr int scale10 {fx_scale10(M)};
    using digits_type = fx_digits_t<M>;

    const std::uint64_t whole {gamma_whole(t)};
    if (whole >= bound)
    {
        return {0U, 0, neg ? exp_detail::exp_kind::underflow : exp_detail::exp_kind::overflow};
    }
    auto f {fx_top<M>(fx_frac<N>(t.w))};
    auto exp10 {static_cast<int>(whole)};
    if (neg)
    {
        exp10 = -exp10 - 1;
        f = fx_sub(fx_one<M>(), f);
    }
    const auto r10 {exp_detail::exp10_frac<M, Wide>(f)};
    auto a {r10.a};
    auto m {fx_mul(r10.m, fx_top<M>(v))};
    if (m.w[M - 1] >= (UINT64_C(1) << 63U))
    {
        m = fx_shr(m, 1U);
        ++a;
    }

    // From 10 on, one digit less
    const int ten {a == 4U || (a == 3U && m.w[M - 1] >= (UINT64_C(5) << 60U)) ? 1 : 0};
    return {static_cast<digits_type>(exp_detail::exp10_digits(exp_detail::exp_mantissa<M> {m, a}, scale10 - ten)), exp10 - scale10 + ten,
            exp_detail::exp_kind::digits};
}

// e^st * c as digits for st = st_abs, or -st_abs with neg_st: 10^t * v with t = st * log10(e) + n2 * log10(2) + n10,
// for c = v * 2^n2 * 10^n10
template <typename T, int N, bool Wide = false, int M = N>
constexpr auto gamma_exp_digits(const wide<N>& st_abs, const bool neg_st, const gamma_float<N>& c) noexcept -> exp_detail::exp_digits<M>
{
    auto t {fx_wide_mul(st_abs, fx_top<N>(fx_table::log10_e))};
    t = fx_wide_add(neg_st ? fx_wide_negate(t) : t, log_detail::log_mul_small<N>(fx_table::log10_2, static_cast<std::uint64_t>(c.n2 < 0 ? -c.n2 : c.n2)), c.n2 < 0);
    t = fx_wide_add(t, fx_wide_int<N>(static_cast<std::uint64_t>(c.n10 < 0 ? -c.n10 : c.n10)), c.n10 < 0);
    bool neg_t {};
    t = fx_wide_abs(t, neg_t);
    return gamma_exp10<N, Wide, M>(t, neg_t, static_cast<std::uint64_t>(std::numeric_limits<T>::max_exponent10 + std::numeric_limits<T>::digits10 + 1), c.v);
}

// The same for a signed st
template <typename T, int N, bool Wide = false, int M = N>
constexpr auto gamma_exp_digits(const wide<N>& st, const gamma_float<N>& c) noexcept -> exp_detail::exp_digits<M>
{
    bool neg_st {};
    const auto st_abs {fx_wide_abs(st, neg_st)};
    return gamma_exp_digits<T, N, Wide, M>(st_abs, neg_st, c);
}

// A result of gamma_exp_digits, rounded once in the current mode with the sign neg
template <typename T, int N>
constexpr auto gamma_exp_round(const exp_detail::exp_digits<N>& r, const bool neg) noexcept -> T
{
    if (r.kind != exp_detail::exp_kind::digits)
    {
        return exp_detail::exp10_beyond<T, N>(r.kind, neg);
    }
    return fx_round<T>(r.digits, r.exp10, neg);
}

// c as digits: 2^n2 = 10^q * m * 2^s for 10^-q = m * 2^(s - n2) from pow10_bin, with 2^s * m in [1, 10)
constexpr auto gamma_float_digits(gamma_float<1> c) noexcept -> exp_detail::exp_digits<1>
{
    const int q {c.n2 >= 0 ? (c.n2 * 1233) >> 12 : -((4095 - c.n2 * 1233) >> 12)};
    const int s {c.n2 + (q <= 0 ? (-q * 1701) >> 9 : -((511 + q * 1701) >> 9))};
    c = gamma_fmul(gamma_float<1> {c.v, 0, c.n10 + q}, gamma_float<1> {fx_top<1>(fx_table::pow10_bin[38 - q]), s, 0});
    const int ten {c.n2 >= 4 || (c.n2 == 3 && c.v.w[0] >= (UINT64_C(5) << 60U)) ? 1 : 0};
    return {static_cast<std::uint64_t>(exp_detail::exp10_digits(exp_detail::exp_mantissa<1> {c.v, static_cast<unsigned>(c.n2)}, 18 - ten)), c.n10 - 18 + ten,
            exp_detail::exp_kind::digits};
}

// gamma(x) for 0 < |x| = sig * 10^e < 20 for decimal64_t in one word, rounded once in out where it cannot be next to a
// rounding boundary: gamma_float_digits adds 4 units of 2^-62 to the error of c, and for x < 0 the reflection 8 more
template <typename T>
constexpr auto tgamma_try(const boost::int128::uint128_t& sig, const int e, const bool neg_x, T& out) noexcept -> bool
{
    auto y {gamma_wide<1>(sig, e)};
    if (!neg_x)
    {
        const auto r {gamma_float_digits(gamma_poly_one(y, sig, e))};
        return fx_try_round(r.digits, r.exp10, false, ((r.digits >> 58U) + 1U) * (gamma_poly_error(gamma_whole(y)) + 4U) / 16U + 2U, out);
    }
    // gamma(x) = 1 / (|sin(pi x)| / pi * gamma(1 - x)), with the sine in two words
    const auto g {gamma_reflect<2>(sig, e)};
    y = fx_wide_add(y, fx_wide_int<1>(1U));
    const auto r {gamma_float_digits(gamma_frecip(gamma_fmul(gamma_ftop(g.s), gamma_poly_one(y, sig, e))))};
    return fx_try_round(r.digits, r.exp10, g.neg, ((r.digits >> 58U) + 1U) * (gamma_poly_error(gamma_whole(y)) + 12U) / 16U + 2U, out);
}

// gamma(x) for x = sig * 10^e < 10^16 in N words, rounded once in the current mode: exp(ln|gamma(x)|) with the sign, or
// for |x| < 12 (20) from gamma_poly_parts. A check of all decimal32_t values shows that one word rounds right.
template <typename T, int N>
constexpr auto tgamma_words(const boost::int128::uint128_t& sig, const int e, const bool neg_x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    BOOST_DECIMAL_IF_CONSTEXPR (p == 16)
    {
        T r {};
        if (!decimal_at_least(sig, e, p, 20U, 2) && tgamma_try(sig, e, neg_x, r))
        {
            return r;
        }
        // For |x| >= 20, a first try with the Stirling head in two words and its tail and 10^t in one: its error is below
        // 40 units of 2^-62 for x > 0, else 48, so below d1 * 2^-57 * (5 / 4 or 3 / 2) + 2 units of d1
        if (decimal_at_least(sig, e, p, 20U, 2))
        {
            auto y {gamma_wide<2>(sig, e)};
            gamma_sine<2> g {{fx_one<2>(), 0, 0}, false};
            if (neg_x)
            {
                g = gamma_reflect<2>(sig, e);
                y = fx_wide_add(y, fx_wide_int<2>(1U));
            }
            const auto f {gamma_normalize(y)};
            const auto head {gamma_stirling_head(y, f)};
            const auto st1 {fx_wide_add(head, fx_widen(fx_extend<2>(gamma_stirling_tail<1>(f))))};
            exp_detail::exp_digits<1> r1 {};
            auto c {g.s};
            if (neg_x)
            {
                // gamma(x) = e^-st * c for st = ln gamma(1 - x) and c = 1 / (|sin(pi x)| / pi)
                c = gamma_fmul(gamma_float<2> {fx_one<2>(), 0, 0}, gamma_frecip(g.s));
                r1 = gamma_exp_digits<T, 2, false, 1>(fx_wide_negate(st1), c);
            }
            else
            {
                const auto t {fx_wide_mul(st1, fx_top<2>(fx_table::log10_e))};
                const auto whole {gamma_whole(t)};
                r1 = {static_cast<std::uint64_t>(exp_detail::exp10_digits(exp_detail::exp10_frac(fx_top<1>(fx_frac<2>(t.w))), 18)), static_cast<int>(whole) - 18,
                      whole < static_cast<std::uint64_t>(std::numeric_limits<T>::max_exponent10 + p + 1) ? exp_detail::exp_kind::digits : exp_detail::exp_kind::overflow};
            }
            if (r1.kind == exp_detail::exp_kind::digits && fx_try_round(r1.digits, r1.exp10, g.neg, (r1.digits >> 57U) + (r1.digits >> (neg_x ? 58U : 59U)) + 2U, r))
            {
                return r;
            }
            // Else the tail in two words
            const auto st {fx_wide_add(head, fx_widen(gamma_stirling_tail<2>(f)))};
            return gamma_exp_round<T>(gamma_exp_digits<T>(neg_x ? fx_wide_negate(st) : st, neg_x ? c : gamma_frecip(c)), g.neg);
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
    {
        // For |x| < 20, a first try in two words with the wide series of exp: its error is below 220 units of the
        // digits for x > 0, else 420
        if (!decimal_at_least(sig, e, p, 20U, 2))
        {
            const auto g2 {gamma_poly_parts<2>(sig, e, neg_x, p)};
            const auto r2 {gamma_exp_digits<T, 2, true>(fx_wide_negate(g2.st), g2.c)};
            if (!fx_near_boundary(r2.digits, neg_x ? 420U : 220U))
            {
                return fx_round<T>(r2.digits, r2.exp10, g2.neg);
            }
        }
    }
    // |gamma(x)| = e^(+-st) * c with c or 1 / c
    const auto g {gamma_parts_of<N>(sig, e, neg_x, p)};
    const auto c {g.neg_st ? g.c : gamma_frecip(g.c)};
    // In two words, gamma_poly_parts gives st < 0 for x > 0
    const auto st {g.neg_st ? fx_wide_negate(g.st) : g.st};
    const auto r {N == 2 ? gamma_exp_digits<T>(st, c) : gamma_exp_digits<T>(g.st, g.neg_st, c)};
    return gamma_exp_round<T>(r, g.neg);
}

// gamma(x) = 1/x - euler + O(x) for a tiny x = sig * 10^e: with k = 2p + 1 and b = floor(euler sig 10^(k + e)), |10^k / x
// - euler| is in [D, D + 1] units 10^(-e - k) with D = floor((10^k - b - 1) / sig), or floor((10^k + b) / sig) for x < 0
template <typename T>
constexpr auto tgamma_tiny(const boost::int128::uint128_t& sig, const int e, const bool neg_x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int k {2 * p + 1};
    using digits_type = fx_digits_t<words>;
    using wide_type = std::conditional_t<words == 3, u256, boost::int128::uint128_t>;

    const auto euler {fx_sub(fx_one<words>(), fx_top<words>(fx_table::one_minus_euler))};
    const auto b {static_cast<wide_type>(fx_scale(sig, fx_mul(euler, fx_from_decimal<words>(1U, -k - e))))};
    const auto ten {pow10(static_cast<wide_type>(k))};
    const auto d {static_cast<boost::int128::uint128_t>((neg_x ? ten + b : ten - b - 1U) / static_cast<wide_type>(sig))};
    return fx_round<T>(static_cast<digits_type>(d * 10U + 1U), -e - k - 1, neg_x);
}

// gamma(x) for x = sig * 10^e (negative with neg_x) not 0 or a negative integer, rounded once in the current mode,
// with (x - 1)! for an integer x where it fits in 128 bits
template <typename T>
constexpr auto tgamma_fixed(const boost::int128::uint128_t& sig, const int e, const bool neg_x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    if (e + p > 16)
    {
        return neg_x ? T {1U, -2 * std::numeric_limits<T>::max_exponent10, gamma_reflect<words>(sig, e).neg}
                     : T {1U, 2 * std::numeric_limits<T>::max_exponent10};
    }
    if (e + p <= -(p + (words == 1 ? 10 : words == 2 ? 11 : 16)))
    {
        return tgamma_tiny<T>(sig, e, neg_x);
    }
    if (!neg_x && gamma_is_integer(sig, e))
    {
        const auto n {static_cast<std::uint64_t>(gamma_integer_of(sig, e))};
        if (n <= 37U)
        {
            // 34! has seven trailing zeros, so 35! and 36! fit in 128 bits without them
            return n <= 35U ? T {gamma_at(n), 0} : T {gamma_at(35U) / 10000000U * (n == 36U ? 35U : 35U * 36U), 7};
        }
    }
    return tgamma_words<T, words>(sig, e, neg_x);
}

} // namespace gamma_detail
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_GAMMA_IMPL_HPP
