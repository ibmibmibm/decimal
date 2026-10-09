// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_TANH_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_TANH_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/exp.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace exp_detail {

// For a small |x|, tanh(x) = x * (sinh(x) / x) / cosh(x) by the Taylor series of both in x^2, without the x
template <int words, std::size_t terms>
constexpr auto tanh_small(const fx<words>& frac) noexcept -> fx<words>
{
    fx<words> s[terms] {};
    fx<words> c[terms] {};
    exp_head<1U, 2U>(s);
    exp_head<0U, 2U>(c);
    const auto u {fx_mul(frac, frac)};
    return fx_div(fx_horner(u, s), fx_horner(u, c));
}

// The digits of tanh(x) = m (1 - x^2 R(x^2)) next to a rounding boundary, with x = t * 10^-a for t in [0.1, 1):
// R = sum 2k / (2k + 1)! x^(2k - 2) / cosh(x) and the rest t^2 R 10^-2a take their own precision
template <int words, std::size_t terms>
constexpr auto tanh_lead(const boost::int128::uint128_t& sig, const int k, const fx<words>& frac, const boost::int128::uint128_t& m) noexcept -> boost::int128::uint128_t
{
    fx<words> c[terms] {};
    fx<words> d[terms - 1U] {};
    fx<words> s[terms - 1U] {};
    exp_head<0U, 2U>(c);
    exp_head<3U, 2U>(s);
    for (std::size_t i {}; i < terms - 1U; ++i)
    {
        d[i] = fx_sub(c[i + 1U], s[i]);
    }
    const auto u {fx_mul(frac, frac)};
    int a {};
    const auto t {fx_lead_t<words>(sig, k, a)};
    const auto g {fx_mul(fx_mul(t, t), fx_div(fx_horner(u, d), fx_horner(u, c)))};
    return fx_lead_digits(m, g, 2 * a, true);
}

// tanh(|x|) = (1 - t) / (1 + t) in q for |x| = whole + frac >= 1/16 and t = exp(-2|x|), or false where t is below
// 10^-(19 words): then tanh(|x|) is closer to 1 than the rounding can see
template <int words, bool Precise = false>
BOOST_DECIMAL_FORCE_INLINE constexpr auto tanh_quotient(const std::uint64_t whole, const fx<words>& frac, fx<words>& q) noexcept -> bool
{
    // With 2y = n + h and 0 <= h < 1, t = 10^-(n + 1) * 10^(1 - h)
    const auto y {exp_log10(whole, frac, fx_table::log10_e)};
    auto h {fx_add(y.f, y.f)};
    int n {2 * y.exp10};
    // h >= 1
    if (h.w[words - 1] >= (UINT64_C(1) << 62U))
    {
        h = fx_sub(h, fx_one<words>());
        ++n;
    }
    if (n + 1 >= 19 * words)
    {
        return false;
    }

    // 2t is larger than the error of 1 / (1 + t), thus q < 1
    const auto r {exp10_frac<words, Precise>(fx_sub(fx_one<words>(), h))};
    const auto t {fx_mul(fx_shr(r.m, 3U - r.a), fx_from_decimal<words>(8U, n + 1))};
    q = fx_div(fx_sub(fx_one<words>(), t), fx_add(fx_one<words>(), t));
    return true;
}

// For decimal64_t, tanh(|x|) in one word, rounded once in out where it cannot be next to a rounding boundary. Below
// 1/2, the series take 10 terms, and the quotient has an error below 8 units of 2^-62.
template <typename T>
constexpr auto tanh_try(const boost::int128::uint128_t& sig, const int e, const std::uint64_t whole, const bool neg, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;

    const auto frac {exp_frac<1>(sig, -e, whole)};
    if (whole == 0U && frac.w[0] < (UINT64_C(1) << 61U))
    {
        // For |x| < 10^-8 the digits are next to m
        return e > -24 && exp_try_round(sig * 1000U, tanh_small<1, 10U>(frac), e - 3, neg, 8U, out);
    }

    // 1 - h has an error below 7 units of 2^-62, thus 10^(1 - h) and t a relative error below 24 units and t 4 units
    // more. q > 0.46 then has an error below 12 (1 - q^2) + 9 units, and its 19 digits below 27 (1 - q^2) + 21 units.
    fx<1> q {};
    if (!tanh_quotient(whole, frac, q))
    {
        // t is below 10^-19, thus tanh(|x|) is as close to 1 for 16 digits
        constexpr int words {fx_words(std::numeric_limits<T>::digits10)};
        out = just_below_one<T, fx_digits_t<words>>(fx_scale10(words), neg);
        return true;
    }
    const auto d {fx_scale(pow10(static_cast<sig_type>(19)), q)};
    const auto rest {(UINT64_C(1) << 62U) - fx_mul(q, q).w[0]};
    return fx_try_round(d.low, -19, neg, ((rest >> 4U) * 27U >> 58U) + 22U, out);
}

// For decimal128_t and |x| >= 1/2, tanh(|x|) in two words with the precise series of exp, rounded once in out where it
// cannot be next to a rounding boundary. Below 1/2, two words of the series are not faster than three.
template <typename T>
constexpr auto tanh_precise_try(const boost::int128::uint128_t& sig, const int e, const std::uint64_t whole, const bool neg, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;
    constexpr int words {fx_words(std::numeric_limits<T>::digits10)};
    using digits_type = fx_digits_t<words>;

    const auto frac {exp_frac<2>(sig, -e, whole)};

    // t has an error below 61 t + 2.25 units, and q below 124 t + 14.4, with t <= 1 - q
    fx<2> q {};
    if (!tanh_quotient<2, true>(whole, frac, q))
    {
        out = just_below_one<T, digits_type>(fx_scale10(words), neg);
        return true;
    }
    const auto one {pow10(static_cast<sig_type>(fx_scale10(2) + 1))};
    const auto d {fx_scale(one, q)};
    if (fx_near_boundary(d, static_cast<std::uint64_t>((one - d) >> 118U) + 24U))
    {
        return false;
    }
    out = fx_round<T>(static_cast<digits_type>(d), -fx_scale10(2) - 1, neg);
    return true;
}

// tanh(|x|) for a finite x != 0, rounded once in the current mode with the sign neg
template <typename T>
constexpr auto tanh_fixed(const T x, const bool neg) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};

    const int k {-e};
    const auto wide_whole {exp_whole(sig, k)};
    const auto whole {static_cast<std::uint64_t>(wide_whole)};
    // From |x| >= 22 words, 2|x| log10(e) > 19 words: t is below the bound of the test after exp_log10
    if (wide_whole >= 22U * words)
    {
        return just_below_one<T, digits_type>(scale10, neg);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        T r {};
        if (tanh_try(sig, e, whole, neg, r))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        T r {};
        if (decimal_at_least(sig, e, p, 50U, 0) && tanh_precise_try(sig, e, whole, neg, r))
        {
            return r;
        }
    }
    const auto frac {exp_frac<words>(sig, k, whole)};

    if (exp_is_small(whole, frac))
    {
        const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(scale10 - p))};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // A first try in two words where |x| >= 10^-16, else the digits are next to m: its error is
            // below 128 units of 2^-126, so below d2 * 2^-118 + 2 units of d2
            if (e + p > -16)
            {
                const auto d2 {fx_scale(m, tanh_small<2, 10U>(exp_frac<2>(sig, k, whole)))};
                if (d2 >= pow10(static_cast<boost::int128::uint128_t>(fx_scale10(2) - 1)) && !fx_near_boundary(d2, static_cast<std::uint64_t>(d2 >> 118U) + 2U))
                {
                    return fx_round<T>(static_cast<digits_type>(d2), e - (scale10 - p), neg);
                }
            }
        }

        // tanh(|x|) < |x| = m * 10^e, thus the digits are below m also where the quotient rounds to 1
        const auto digits {fx_scale(m, tanh_small<words, exp_small_terms(words)>(frac))};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            if (fx_near_boundary(digits, 2U))
            {
                return fx_round<T>(static_cast<digits_type>(tanh_lead<words, exp_small_terms(words)>(sig, k, frac, m)), e - (scale10 - p), neg);
            }
        }
        return fx_round<T>(static_cast<digits_type>(digits < m ? digits : m - 1U), e - (scale10 - p), neg);
    }

    fx<words> q {};
    if (!tanh_quotient(whole, frac, q))
    {
        return just_below_one<T, digits_type>(scale10, neg);
    }
    return fx_round<T>(static_cast<digits_type>(fx_scale(pow10(static_cast<boost::int128::uint128_t>(scale10 + 1)), q)), -scale10 - 1, neg);
}

} // namespace exp_detail

template <typename T>
constexpr auto tanh_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpc = fpclassify(x);

    constexpr T zero { 0, 0 };

    auto result = zero;

    if (fpc == FP_ZERO)
    {
        result = x;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpc != FP_NORMAL && fpc != FP_SUBNORMAL)
    {
        if (fpc == FP_INFINITE)
        {
            if (signbit(x))
            {
                result = -T {1, 0};
            }
            else
            {
                result = T {1, 0};
            }
        }
        else if (fpc == FP_NAN)
        {
            result = issignaling(x) ? nan_conversion(x) : x;
        }
    }
    #endif
    else
    {
        result = exp_detail::tanh_fixed(x, signbit(x));
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto tanh(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::tanh_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_TANH_HPP
