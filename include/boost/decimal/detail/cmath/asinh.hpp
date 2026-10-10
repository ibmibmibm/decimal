// Copyright 2023 Matt Borland
// Copyright 2023 - 2026 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ASINH_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ASINH_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/concepts.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace log_detail {

// asinh(x) / x = 1 - t/6 + 3t^2/40 - ... for t = x^2 < 2^-16, to below 2^-69, 2^-133 and 2^-197 in one, two
// and three words. The terms from t^4 and t^8 on are below 2^-64 and 2^-128.
template <int N>
constexpr auto asinh_series(const fx<N>& t) noexcept -> fx<N>
{
    return fx_series<true, N == 1 ? 4U : N == 2 ? 8U : 12U, 4U, 8U>(t, fx_table::asinh_coeff);
}

// The s with sig * 2^s in [0.5, 1) * 10^p, for sig of p digits
constexpr auto asinh_shift(const boost::int128::uint128_t& sig, const int p) noexcept -> unsigned
{
    const auto q {pow10(static_cast<boost::int128::uint128_t>(p - 1))};
    return static_cast<unsigned>(sig < q * 5U) + static_cast<unsigned>(sig * 2U < q * 5U) + static_cast<unsigned>(sig * 4U < q * 5U);
}

// asinh(x) = ln(|x| + sqrt(x^2 + 1)) for |x| >= 1/256. Below 1, it is ln(u + sqrt(u^2 + 1)) with u = |x|. Else
// |x| = u * 2^-s * 10^n10 with u in [0.5, 1), and it is ln((u + sqrt(u^2 + 4^s * 10^(-2 n10))) * 2^-s * 10^n10).
template <int N, int Steps = 2>
constexpr auto asinh_ln(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> fx_wide<N>
{
    if (e + p <= 0)
    {
        const auto u {fx_from_decimal<N>(sig, -e)};
        return log_wide_of(fx_add(u, fx_sqrt<N, Steps>(fx_add(fx_one<N>(), fx_mul(u, u)))), 0, 0);
    }
    const unsigned s {asinh_shift(sig, p)};
    const auto u {fx_from_decimal<N>(sig << s, p)};
    const auto c {fx_from_decimal<N>(1U << (2U * s), 2 * (e + p))};
    return log_wide_of(fx_add(u, fx_sqrt<N, Steps>(fx_add(fx_mul(u, u), c))), -static_cast<int>(s), e + p);
}

// asinh of a finite x != 0, rounded once in the current mode
template <typename T>
constexpr auto asinh_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // For |x| < 1/256, asinh(x) = x * asinh(x) / x < |x|, also where the series rounds to 1
    if (decimal_tiny(sig, e))
    {
        const auto u {fx_from_decimal<words>(sig, -e)};
        return fx_round_below_x<T>(sig, e, asinh_series(fx_mul(u, u)), neg);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer for |x| >= 0.2 with an error below 256 units of its last bit, and in one word
        // for |x| >= 0.11 with an error in (-10.8, 0.6) units of 2^-62 below 1 and in (-16, 4.4) from 1 on
        T r {};
        if (decimal_at_least(sig, e, p, words == 3 ? 20U : 11U, 0)
            && fx_wide_try_round(asinh_ln<words == 3 ? 2 : 1>(sig, e, p), neg, r, words == 3 ? 8U : 4U))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
    {
        // A first try with one Newton step fewer in the sqrt, within 2^28 units of 2^-62 of the full value for every
        // decimal32_t x, for |x| >= 0.11, where the result is above 0.1
        T r {};
        if (decimal_at_least(sig, e, p, 11U, 0) && fx_wide_try_round(asinh_ln<1, 1>(sig, e, p), neg, r, 29U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(asinh_ln<words>(sig, e, p), neg);
}

} // namespace log_detail

template <typename T>
constexpr auto asinh_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    T result { };

    const auto fpc {fpclassify(x)};

    if (fpc != FP_NORMAL && fpc != FP_SUBNORMAL)
    {
        result = issignaling(x) ? nan_conversion(x) : x;
    }
    else
    {
        result = log_detail::asinh_fixed(x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto asinh(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::asinh_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_ASINH_HPP
