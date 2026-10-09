// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_SINH_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_SINH_HPP

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

// The Taylor series of sinh(x) / x or cosh(x) in x^2
template <int words, std::size_t terms>
constexpr auto sinh_cosh_series(const fx<words>& frac, const bool is_cosh) noexcept -> fx<words>
{
    fx<words> c[terms] {};
    if (is_cosh)
    {
        exp_head<0U, 2U>(c);
    }
    else
    {
        exp_head<1U, 2U>(c);
    }
    return fx_horner(fx_mul(frac, frac), c);
}

// The digits of sinh(x) or cosh(x) = m (1 + x^2 R(x^2)) next to a rounding boundary, with x = t * 10^-a for t in
// [0.1, 1): the rest x^2 R = t^2 R 10^-2a takes its own precision
template <int words, std::size_t terms>
constexpr auto sinh_cosh_lead(const boost::int128::uint128_t& sig, const int k, const fx<words>& frac, const boost::int128::uint128_t& m, const bool is_cosh) noexcept -> boost::int128::uint128_t
{
    fx<words> c[terms - 1U] {};
    if (is_cosh)
    {
        exp_head<2U, 2U>(c);
    }
    else
    {
        exp_head<3U, 2U>(c);
    }
    int a {};
    const auto t {fx_lead_t<words>(sig, k, a)};
    const auto g {fx_mul(fx_mul(t, t), fx_horner(fx_mul(frac, frac), c))};
    return fx_lead_digits(m, g, 2 * a, false);
}

// (exp(|x|) +- exp(-|x|)) / 2 = 10^exp10 * 4d for |x| = whole + frac, with exp10 and d in the parts
template <int words, bool Precise = false>
BOOST_DECIMAL_FORCE_INLINE constexpr auto sinh_cosh_parts(const std::uint64_t whole, const fx<words>& frac, const bool is_cosh) noexcept -> exp_parts<words>
{
    // exp(|x|) = 10^exp10 * 10^f and exp(-|x|) = 10^-gap * 10^exp10 * 10^(1 - f) with gap = 2 exp10 + 1.
    // In d, each 10^f is divided by 8 so that it is below 4, and 10^(1 - f) / 8 = 10 / (8 * 2^a * m) = 5/4 * 2^-a / m.
    const auto y {exp_log10(whole, frac, fx_table::log10_e)};
    const auto big {exp10_frac<words, Precise>(y.f)};
    auto d {fx_shr(big.m, 3U - big.a)};

    // Where 10^-gap is below the last bit of d, exp(-|x|) has no effect
    const int gap {2 * y.exp10 + 1};
    if (gap < 19 * words)
    {
        const auto small {fx_mul_shr(fx_recip(big.m), 5U, 2U)};
        const auto low {fx_mul(fx_shr(small, big.a), fx_from_decimal<words>(1U, gap))};
        d = is_cosh ? fx_add(d, low) : fx_sub(d, low);
    }
    return {d, y.exp10};
}

// For decimal64_t, sinh(|x|) or cosh(x) in one word, rounded once in out where it cannot be next to a rounding
// boundary. Below 1, sinh(x) = x S(x^2) with 10 terms of S: its error is below 4 units of 2^-62 for |x| < 1/2, else 20.
template <typename T>
constexpr auto sinh_cosh_try(const boost::int128::uint128_t& sig, const int e, const std::uint64_t whole, const bool is_cosh, const bool neg, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;

    const auto frac {exp_frac<1>(sig, -e, whole)};
    if (whole == 0U && (!is_cosh || exp_is_small(whole, frac)))
    {
        // For |x| < 10^-8 the digits are next to m. cosh takes the 5 terms below, with an error below 2^-60.
        const bool half {frac.w[0] < (UINT64_C(1) << 61U)};
        return e > -24 && (is_cosh ? exp_try_round(pow10(static_cast<sig_type>(18)), sinh_cosh_series<1, 5U>(frac, true), -18, neg, 4U, out)
                                   : exp_try_round(sig * 1000U, sinh_cosh_series<1, 10U>(frac, false), e - 3, neg, half ? 4U : 20U, out));
    }

    // 10^f has a relative error below 2^-58, the shift adds 1 unit of 2^-62 and exp(-|x|) below 5 for gap = 1, else 3.
    // For sinh, |x| >= 1 here, and exp(|x|) is below 1.16 sinh(|x|).
    const auto r {sinh_cosh_parts(whole, frac, is_cosh)};
    const std::uint64_t err {r.exp10 == 0 ? 6U : r.exp10 < 9 ? 4U : 1U};

    // 19 digits of 4d > 1/4: 4d * 10^18 for 4d >= 1, else 4d * 10^19
    const bool below {r.f.w[0] < (UINT64_C(1) << 60U)};
    const sig_type m {pow10(static_cast<sig_type>(below ? 19 : 18)) * 4U};
    const auto d {fx_scale(m, r.f)};
    const auto w {(d >> 58U) + (is_cosh ? sig_type {} : (d >> 61U) + (d >> 62U) + 2U) + ((m * err) >> 62U) + 3U};
    return fx_try_round(d.low, r.exp10 - (below ? 19 : 18), neg, w.low, out);
}

// For decimal128_t and |x| >= 1/16, sinh(|x|) or cosh(x) in two words with the precise series of exp, rounded once in
// out where it cannot be next to a rounding boundary. Below 1, sinh(x) = x S(x^2) with 17 terms of S.
template <typename T>
constexpr auto sinh_cosh_precise_try(const boost::int128::uint128_t& sig, const int e, const std::uint64_t whole, const bool is_cosh, const bool neg, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;
    using digits_type = fx_digits_t<fx_words(std::numeric_limits<T>::digits10)>;
    constexpr int scale10 {fx_scale10(2)};
    constexpr int up {scale10 - std::numeric_limits<T>::digits10};

    const auto frac {exp_frac<2>(sig, -e, whole)};
    if (exp_is_small(whole, frac))
    {
        return false;
    }
    if (whole == 0U && !is_cosh)
    {
        // S has an error below 34 units of 2^-126, so the digits below d * 2^-120 + 2 units
        const auto d {fx_scale(sig * pow10(static_cast<sig_type>(up)), sinh_cosh_series<2, 17U>(frac, false))};
        if (fx_near_boundary(d, static_cast<std::uint64_t>(d >> 120U) + 2U))
        {
            return false;
        }
        out = fx_round<T>(static_cast<digits_type>(d), e - up, neg);
        return true;
    }

    // 10^f has a relative error below 46 units of 2^-126, and 4d below 110 for cosh and for sinh with |x| >= 1
    const auto r {sinh_cosh_parts<2, true>(whole, frac, is_cosh)};
    const auto d {fx_scale(pow10(static_cast<sig_type>(scale10)) * 4U, r.f)};
    if (fx_near_boundary(d, static_cast<std::uint64_t>(d >> 119U) + 2U))
    {
        return false;
    }
    out = fx_round<T>(static_cast<digits_type>(d), r.exp10 - scale10, neg);
    return true;
}

// sinh(|x|) for is_cosh false and cosh(x) for is_cosh true, for a finite x != 0, rounded once in the
// current mode with the sign neg
template <typename T>
constexpr auto sinh_cosh_fixed(const T x, const bool is_cosh, const bool neg) noexcept -> T
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
    if (wide_whole >= exp_bound<T>())
    {
        return T {1U, 2 * std::numeric_limits<T>::max_exponent10, neg};
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        T r {};
        if (sinh_cosh_try(sig, e, whole, is_cosh, neg, r))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        T r {};
        if (sinh_cosh_precise_try(sig, e, whole, is_cosh, neg, r))
        {
            return r;
        }
    }
    const auto frac {exp_frac<words>(sig, k, whole)};

    // For |x| < 1/16, the Taylor series sinh(x) = x * (1 + x^2/3! + ...) and cosh(x) = 1 + x^2/2! + ...
    // The error is below 2^-60, 2^-100 and 2^-160 for one, two and three words.
    if (exp_is_small(whole, frac))
    {
        const auto m {is_cosh ? pow10(static_cast<boost::int128::uint128_t>(scale10)) : sig * pow10(static_cast<boost::int128::uint128_t>(scale10 - p))};
        const int exp10 {is_cosh ? -scale10 : e - (scale10 - p)};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // A first try in two words where |x| >= 10^-16, else the digits are next to m: its error is
            // below 64 units of 2^-126, so below d2 * 2^-119 + 2 units of d2
            if (e + p > -16)
            {
                const auto d2 {fx_scale(m, sinh_cosh_series<2, 10U>(exp_frac<2>(sig, k, whole), is_cosh))};
                if (!fx_near_boundary(d2, static_cast<std::uint64_t>(d2 >> 119U) + 2U))
                {
                    return fx_round<T>(static_cast<digits_type>(d2), exp10, neg);
                }
            }
        }
        const auto d {fx_scale(m, sinh_cosh_series<words, exp_small_terms(words)>(frac, is_cosh))};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            if (fx_near_boundary(d, 2U))
            {
                return fx_round<T>(static_cast<digits_type>(sinh_cosh_lead<words, exp_small_terms(words)>(sig, k, frac, m, is_cosh)), exp10, neg);
            }
        }
        return fx_round<T>(static_cast<digits_type>(d), exp10, neg);
    }

    const auto r {sinh_cosh_parts(whole, frac, is_cosh)};
    return fx_round<T>(static_cast<digits_type>(fx_scale(pow10(static_cast<boost::int128::uint128_t>(scale10)) * 4U, r.f)), r.exp10 - scale10, neg);
}

} // namespace exp_detail

template <typename T>
constexpr auto sinh_impl(const T x) noexcept
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
            result = x;
        }
        else if (fpc == FP_NAN)
        {
            result = issignaling(x) ? nan_conversion(x) : x;
        }
    }
    #endif
    else
    {
        result = exp_detail::sinh_cosh_fixed(x, false, signbit(x));
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto sinh(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::sinh_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_SINH_HPP
