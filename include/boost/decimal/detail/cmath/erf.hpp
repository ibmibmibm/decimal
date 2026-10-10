// Copyright John Maddock 2006.
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP

#include <boost/decimal/detail/cmath/impl/exp_impl.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace erf_detail {

// The polynomial of g on interval i for the given words: its first coefficient in the table, the number of
// coefficients, its center c = cw + 1/2 for half, log2 of its half width h, and the first coefficient of the tail
struct erf_interval
{
    int first;
    int count;
    unsigned cw;
    int log2_h;
    int split;
    bool half;
};

constexpr auto erf_interval_of(const int words, const unsigned i) noexcept -> erf_interval
{
    // [1, 2], [2, 4] and [4, 8] for all types, then [8, 10] for two words, and [8, 12] and [12, 16] for three words
    return words == 1 ? (i == 0U   ? erf_interval {0, 20, 1U, -1, 20, true}
                         : i == 1U ? erf_interval {20, 23, 3U, 0, 23, false}
                                   : erf_interval {43, 26, 6U, 1, 26, false})
         : words == 2 ? (i == 0U   ? erf_interval {0, 31, 1U, -1, 23, true}
                         : i == 1U ? erf_interval {31, 36, 3U, 0, 28, false}
                         : i == 2U ? erf_interval {67, 40, 6U, 1, 33, false}
                                   : erf_interval {107, 27, 9U, 0, 19, false})
                      : (i == 0U   ? erf_interval {0, 43, 1U, -1, 23, true}
                         : i == 1U ? erf_interval {43, 50, 3U, 0, 28, false}
                         : i == 2U ? erf_interval {93, 57, 6U, 1, 33, false}
                         : i == 3U ? erf_interval {150, 49, 10U, 1, 25, false}
                                   : erf_interval {199, 44, 14U, 1, 21, false});
}

// Whether each first is the sum of the counts before it, and the last interval ends at the size of the table
constexpr auto erf_intervals_fit(const int words, const unsigned last, const int size) noexcept -> bool
{
    for (unsigned i {1}; i <= last; ++i)
    {
        if (erf_interval_of(words, i).first != erf_interval_of(words, i - 1U).first + erf_interval_of(words, i - 1U).count)
        {
            return false;
        }
    }
    return erf_interval_of(words, last).first + erf_interval_of(words, last).count == size;
}

static_assert(erf_intervals_fit(1, 2U, static_cast<int>(sizeof(fx_table::erfc_d32) / sizeof(fx_table::erfc_d32[0]))), "the intervals cover the table");
static_assert(erf_intervals_fit(2, 3U, static_cast<int>(sizeof(fx_table::erfc_d64) / sizeof(fx_table::erfc_d64[0]))), "the intervals cover the table");
static_assert(erf_intervals_fit(3, 4U, static_cast<int>(sizeof(fx_table::erfc_d128) / sizeof(fx_table::erfc_d128[0]))), "the intervals cover the table");

// Whether the tail of each interval, the coefficients from its split on, sums to below 2^-62: for |u| <= 1 its
// partial sums then fit in N - 1 words scaled by 2^64
template <int N, std::size_t M>
constexpr auto erf_splits_fit(const fx<N> (&c)[M], const unsigned last) noexcept -> bool
{
    for (unsigned i {}; i <= last; ++i)
    {
        const auto iv {erf_interval_of(N, i)};
        boost::int128::uint128_t sum {};
        for (int k {iv.first + iv.split}; k < iv.first + iv.count; ++k)
        {
            if (c[k].w[N - 1] != 0U)
            {
                return false;
            }
            sum += c[k].w[N - 2] + UINT64_C(1);
        }
        if (sum.high != 0U)
        {
            return false;
        }
    }
    return true;
}

static_assert(erf_splits_fit(fx_table::erfc_d64, 3U) && erf_splits_fit(fx_table::erfc_d128, 4U), "the tails fit in fewer words");

// The coefficients of the polynomial of the interval iv in the table of N words
constexpr auto erf_poly(const fx<1>&, const erf_interval& iv) noexcept -> fx_span<const fx<1>> { return fx_make_span(fx_table::erfc_d32, iv.first, iv.count); }
constexpr auto erf_poly(const fx<2>&, const erf_interval& iv) noexcept -> fx_span<const fx<2>> { return fx_make_span(fx_table::erfc_d64, iv.first, iv.count); }
constexpr auto erf_poly(const fx<3>&, const erf_interval& iv) noexcept -> fx_span<const fx<3>> { return fx_make_span(fx_table::erfc_d128, iv.first, iv.count); }

// erf(x) / x for z = x^2 <= 1 by the Taylor series, whose terms alternate in sign. The terms from z^20, z^33 and z^45
// on are below 2^-64, 2^-128 and 2^-192. As z can be near 1, all terms take N words.
constexpr auto erf_series(const fx<1>& z) noexcept -> fx<1>
{
    return fx_even_odd(z, fx_make_span<20>(fx_table::erf_taylor), false);
}

template <int N>
constexpr auto erf_series(const fx<N>& z) noexcept -> fx<N>
{
    constexpr std::size_t terms {N == 2 ? 33U : 45U};
    return fx_series<true, terms, terms, terms>(z, fx_table::erf_taylor);
}

// x = whole + frac as x / 4 in fixed point, for x < 16
template <int N>
constexpr auto erf_quarter(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    auto r {fx_shr(frac, 2U)};
    r.w[N - 1] += whole << 60U;
    return r;
}

// The polynomial c in u >= 0, alternating where above, else with positive terms
constexpr auto erf_g_sum(const fx<1>& u, const fx_span<const fx<1>> c, const bool above, std::size_t) noexcept -> fx<1>
{
    return fx_even_odd(u, c, !above);
}

// The tail from c[split] on runs in N - 1 words scaled by 2^64, the low words of c, with the same last bit
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto erf_g_sum(const fx<N>& u, const fx_span<const fx<N>> c, const bool above, const std::size_t split) noexcept -> fx<N>
{
    const auto ut {fx_top<N - 1>(u)};
    auto t {fx_low<N - 1>(c[c.size() - 1U].w)};
    for (std::size_t k {c.size() - 1U}; k-- > split;)
    {
        const auto ck {fx_low<N - 1>(c[k].w)};
        t = above ? fx_sub(ck, fx_mul(ut, t)) : fx_add(ck, fx_mul(ut, t));
    }
    fx<N> s {};
    for (int i {}; i < N - 1; ++i)
    {
        s.w[i] = t.w[i];
    }
    for (std::size_t k {split}; k-- > 0U;)
    {
        s = above ? fx_sub(c[k], fx_mul(u, s)) : fx_add(c[k], fx_mul(u, s));
    }
    return s;
}

// g(x) = erfc(x) exp(x^2) for 1 <= x = whole + frac below the asymptotic range: the polynomial of the interval in
// u = (x - c) / h, which is alternating for u >= 0 and has positive terms in |u| for u < 0
template <int N>
constexpr auto erf_g_table(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    const unsigned i {whole < 2U ? 0U : whole < 4U ? 1U : whole < 8U ? 2U : N == 3 && whole >= 12U ? 4U : 3U};
    const auto iv {erf_interval_of(N, i)};
    const auto x4 {erf_quarter(whole, frac)};
    fx<N> c4 {};
    c4.w[N - 1] = (static_cast<std::uint64_t>(iv.cw) << 60U) | (iv.half ? UINT64_C(1) << 59U : 0U);
    const bool above {!fx_less(x4, c4)};
    const auto u {fx_shl(above ? fx_sub(x4, c4) : fx_sub(c4, x4), static_cast<unsigned>(2 - iv.log2_h))};

    return erf_g_sum(u, erf_poly(u, iv), above, static_cast<std::size_t>(iv.split));
}

// g(x) for x >= 8, 10 or 16 in one, two or three words by the asymptotic series 1 / (x sqrt(pi)) * sum (-1)^n
// (2n - 1)!! / (2x^2)^n, where the first term left out is below 2^-(64N - 2)
template <int N>
constexpr auto erf_g_asymptotic(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    // 1 / x = 1 / d / 2^s with x = d * 2^s and d in [1, 2)
    const auto s {static_cast<unsigned>(63 - countl_zero(whole))};
    auto d {fx_shr(frac, s)};
    d.w[N - 1] += whole << (62U - s);
    const auto inv_x {fx_shr(fx_recip(d), s)};
    const auto v {fx_shr(fx_mul(inv_x, inv_x), 1U)};

    // The sum stops where a term is below 4 units of the last bit; 200 terms bound the loop
    auto term {fx_one<N>()};
    auto pos {term};
    fx<N> neg {};
    for (std::uint64_t n {1U}; n < 200U; ++n)
    {
        const auto t {fx_mul(term, v)};
        term = fx_mul_small(t, 2U * n - 1U);
        bool small {term.w[0] < 4U};
        for (int i {1}; i < N; ++i)
        {
            small = small && term.w[i] == 0U;
        }
        if (small)
        {
            break;
        }
        if ((n & 1U) != 0U)
        {
            neg = fx_add(neg, term);
        }
        else
        {
            pos = fx_add(pos, term);
        }
    }
    return fx_mul(fx_mul(fx_sub(pos, neg), inv_x), fx_top<N>(fx_table::inv_sqrt_pi));
}

// erf(x), or erfc(x) with complement, as digits * 10^exp10 for 0 < |x| < 1: erf(x) = x * s(x^2) by the Taylor series
// keeps all digits, and erfc(x) = 1 - erf(x)
template <typename T, int N>
constexpr auto erf_below_one(const boost::int128::uint128_t& sig, const int e, const bool complement, const bool neg) noexcept -> exp_detail::exp_digits<N>
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int scale10 {fx_scale10(N)};
    using digits_type = fx_digits_t<N>;
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};

    const auto ax {fx_from_decimal<N>(sig, -e)};
    const auto s {erf_series(fx_mul(ax, ax))};
    if (!complement)
    {
        // m has scale10 + 1 digits where s < 1, else scale10, and one more where m * s has a digit fewer
        int len {s.w[N - 1] >= (UINT64_C(1) << 62U) ? scale10 : scale10 + 1};
        const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(len - p))};
        auto d {fx_scale(m, s)};
        if (d < unit)
        {
            ++len;
            d = fx_scale(m * 10U, s);
        }
        return {static_cast<digits_type>(d), e - (len - p), exp_detail::exp_kind::digits};
    }
    const auto er {fx_mul(ax, s)};
    if (neg)
    {
        return {static_cast<digits_type>(fx_scale(unit * 2U, fx_shr(fx_add(fx_one<N>(), er), 1U))), -scale10, exp_detail::exp_kind::digits};
    }
    return {static_cast<digits_type>(fx_scale(unit * 10U, fx_sub(fx_one<N>(), er))), -scale10 - 1, exp_detail::exp_kind::digits};
}

// erf(x), or erfc(x) with complement, in N words for |x| = whole + frac >= 1: erfc(|x|) = exp(-x^2) * g(|x|), erf(x) =
// 1 - erfc(x) and erfc(-|x|) = 2 - erfc(|x|). Kind below_one where erfc(|x|) cannot change the rounding of 1 or 2.
template <typename T, int N>
constexpr auto erf_above_one(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole, const bool complement, const bool neg) noexcept -> exp_detail::exp_digits<N>
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int scale10 {fx_scale10(N)};
    using digits_type = fx_digits_t<N>;
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};
    const auto frac {exp_detail::exp_frac<N>(sig, k, whole)};

    // x^2 = whole^2 + 2 whole frac + frac^2, and exp(-x^2) = 10^(-n - 1) * 10^(1 - f) for x^2 log10(e) = n + f
    std::uint64_t sq[static_cast<std::size_t>(N + 1)] {};
    sq[N] = fx_mul_add(fx_make_span<N>(sq), fx_make_span(frac.w), 2U * whole);
    const auto f2 {fx_mul(frac, frac)};
    sq[N] += fx_add_to(fx_make_span<N>(sq), fx_make_span(f2.w));
    const auto sq_frac {fx_frac<N>(sq)};
    const auto y {exp_detail::exp_log10(whole * whole + ((sq[N] << 2U) | (sq[N - 1] >> 62U)), sq_frac, fx_table::log10_e)};
    int exp10 {-y.exp10 - 1};
    // A first try in two words for decimal128_t takes the wide series of exp(-z), with an error below 2^-124
    const auto r10 {exp_detail::exp10_frac<N, (N == 2 && fx_words(p) == 3)>(fx_sub(fx_one<N>(), y.f))};

    // The end of the last interval of erf_interval_of, where g becomes asymptotic
    constexpr std::uint64_t asymptotic {N == 1 ? 8U : N == 2 ? 10U : 16U};
    const auto g {whole < asymptotic ? erf_g_table(whole, frac) : erf_g_asymptotic(whole, frac)};

    // erfc(|x|) = 16 h * 10^exp10 with h in [1/16, 10/16)
    auto h {fx_shr(fx_mul(r10.m, g), 4U - r10.a)};
    while (h.w[N - 1] < (UINT64_C(1) << 58U))
    {
        h = fx_mul_small(h, 10U);
        --exp10;
    }
    if (complement && !neg)
    {
        return {static_cast<digits_type>(fx_scale(unit * 16U, h)), exp10 - scale10, exp_detail::exp_kind::digits};
    }

    // erfc(|x|) in fixed point. Below 10^-(p + 1) it cannot change the rounding of 1 or 2.
    if (exp10 < -(p + 1))
    {
        return {0U, 0, exp_detail::exp_kind::below_one};
    }
    const auto c {fx_mul(fx_shl(h, 2U), fx_from_decimal<N>(4U, -exp10))};
    if (complement)
    {
        return {static_cast<digits_type>(fx_scale(unit * 2U, fx_shr(fx_sub(fx_int<N>(2U), c), 1U))), -scale10, exp_detail::exp_kind::digits};
    }
    return {static_cast<digits_type>(fx_scale(unit * 10U, fx_sub(fx_one<N>(), c))), -scale10 - 1, exp_detail::exp_kind::digits};
}

// A value in (2 - 10^-scale10, 2), rounded once in the current mode
template <typename T, typename digits_type>
constexpr auto below_two(const int scale10) noexcept -> T
{
    return fx_round<T>(static_cast<digits_type>(pow10(static_cast<boost::int128::uint128_t>(scale10)) * 2U - 1U), -scale10, false);
}

// erf(x), or erfc(x) with complement, for a finite x != 0, rounded once in the current mode. For |x| >= 1, erfc(|x|) =
// exp(-x^2) * g(|x|), and erf(x) = 1 - erfc(x) and erfc(-|x|) = 2 - erfc(|x|) in fixed point.
template <typename T>
constexpr auto erf_fixed(const T x, const bool complement) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};
    const int k {-e};
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};

    if (k >= p)
    {
        // For |x| < 10^-(p + 1), erfc(x) rounds as a value next to 1
        if (complement && k - p >= p + 1)
        {
            return neg ? fx_round<T>(static_cast<digits_type>(unit), -scale10, false) : exp_detail::below_one<T, digits_type>(scale10, false);
        }
        BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
        {
            // A first try in one word: the error of the sum is below 22 units of 2^-62, so below d1 * 2^-57 + 2 units
            // of d1 for erf(x), and 54 or 8 units for erfc(x) with x > 0 or x < 0
            const auto r1 {erf_below_one<T, 1>(sig, e, complement, neg)};
            T r {};
            if (fx_try_round(r1.digits, r1.exp10, !complement && neg, complement ? (neg ? 8U : 54U) : (r1.digits >> 57U) + 2U, r))
            {
                return r;
            }
        }
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // A first try in two words: the error of the sum is below 210 units of 2^-126, so below d2 * 2^-118 + 2
            // units of d2 for erf(x), and 250 units for erfc(x)
            const auto r2 {erf_below_one<T, 2>(sig, e, complement, neg)};
            if (!fx_near_boundary(r2.digits, complement ? 250U : static_cast<std::uint64_t>(r2.digits >> 118U) + 2U))
            {
                return fx_round<T>(r2.digits, r2.exp10, !complement && neg);
            }
        }
        const auto r {erf_below_one<T, words>(sig, e, complement, neg)};
        return fx_round<T>(r.digits, r.exp10, !complement && neg);
    }

    // Above these bounds, erfc(x) is below half of the smallest subnormal value, and erf(x) rounds as a value next to 1
    const auto wide_whole {exp_detail::exp_whole(sig, k)};
    constexpr std::uint64_t tiny_bound {words == 1 ? 16U : words == 2 ? 31U : 120U};
    constexpr std::uint64_t one_bound {words == 1 ? 7U : 10U};
    if (wide_whole >= (complement && !neg ? tiny_bound : one_bound))
    {
        if (!complement)
        {
            return exp_detail::below_one<T, digits_type>(scale10, neg);
        }
        return neg ? below_two<T, digits_type>(scale10) : T {1U, -2 * std::numeric_limits<T>::max_exponent10};
    }
    const auto whole {static_cast<std::uint64_t>(wide_whole)};
    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        // For erf(x) and erfc(-|x|), a first try in one word: the error of erfc(|x|) is below 20 units of 2^-62, so
        // below 46 or 8 units of d1. Where it is too small in one word, it is below 10^-(p + 1) * (1 + 2^-50).
        if (!complement || neg)
        {
            const auto r1 {erf_above_one<T, 1>(sig, k, whole, complement, neg)};
            if (r1.kind != exp_detail::exp_kind::digits)
            {
                return complement ? below_two<T, digits_type>(scale10) : exp_detail::below_one<T, digits_type>(scale10, neg);
            }
            T r {};
            if (fx_try_round(r1.digits, r1.exp10, !complement && neg, complement ? 8U : 46U, r))
            {
                return r;
            }
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // The same in two words: the error of erfc(|x|) is below 128 units of 2^-126, so below 152 or 17 units of d2
        if (!complement || neg)
        {
            const auto r2 {erf_above_one<T, 2>(sig, k, whole, complement, neg)};
            if (r2.kind != exp_detail::exp_kind::digits)
            {
                return complement ? below_two<T, digits_type>(scale10) : exp_detail::below_one<T, digits_type>(scale10, neg);
            }
            if (!fx_near_boundary(r2.digits, complement ? 17U : 152U))
            {
                return fx_round<T>(r2.digits, r2.exp10, !complement && neg);
            }
        }
    }
    const auto r {erf_above_one<T, words>(sig, k, whole, complement, neg)};
    if (r.kind != exp_detail::exp_kind::digits)
    {
        return complement ? below_two<T, digits_type>(scale10) : exp_detail::below_one<T, digits_type>(scale10, neg);
    }
    return fx_round<T>(r.digits, r.exp10, !complement && neg);
}

} // namespace erf_detail

template <typename T>
constexpr auto erf_impl(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // Edge cases
    const auto fp {fpclassify(z)};

    if (fp == FP_ZERO || fp == FP_NAN)
    {
        return z;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fp == FP_INFINITE)
    {
        return z < T{0} ? T{-1} : T{1};
    }
    #endif

    return erf_detail::erf_fixed(z, false);
}

template <typename T>
constexpr auto erfc_impl(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // Edge cases
    const auto fp {fpclassify(z)};

    if (fp == FP_NAN)
    {
        return z;
    }
    else if (fp == FP_ZERO)
    {
        return T{1};
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fp == FP_INFINITE)
    {
        return z < T{0} ? T{2} : T{0};
    }
    #endif

    return erf_detail::erf_fixed(z, true);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto erf(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::erf_impl(static_cast<evaluation_type>(z)));
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto erfc(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::erfc_impl(static_cast<evaluation_type>(z)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP
