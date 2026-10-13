//  Copyright 2002-2011, 2024 Christopher Kormanyos
//  Copyright 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ELLINT_2_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ELLINT_2_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/impl/ellint_impl.hpp>
#include <boost/decimal/detail/cmath/fabs.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/check_non_finite.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto ellint_2_impl(const T k, const T phi) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(k) || isnan(phi))
    {
        return detail::check_non_finite(k, phi);
    }
    #endif

    const auto x {ellint_detail::ellint_args_of(k, phi)};
    if (x.band == 5U && ellint_detail::ellint_above_one(x.ks(), x.ke, precision_v<T>))
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    // The integral grows without bound in phi
    if (isinf(phi))
    {
        return phi;
    }
    #endif

    if (x.ks() == 0U || x.ps() == 0U)
    {
        return phi;
    }
    return ellint_detail::ellint_fixed<T>(x, ellint_detail::ellint_kind::second, signbit(phi));
}

template <typename T>
constexpr auto comp_ellint_2_impl(const T k) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(k))
    {
        return detail::check_non_finite(k);
    }
    #endif

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(k, &e))};
    const auto band {ellint_detail::ellint_band(sig, e, precision_v<T>)};
    if (band == 5U && ellint_detail::ellint_above_one(sig, e, precision_v<T>))
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    // E(1) = 1
    if (band == 5U)
    {
        return T {1};
    }
    return ellint_detail::comp_ellint_fixed<T>(sig, e, band, ellint_detail::ellint_kind::second);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto ellint_2(const T k, const T phi) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::ellint_2_impl(static_cast<evaluation_type>(k), static_cast<evaluation_type>(phi)));
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto comp_ellint_2(const T k) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::comp_ellint_2_impl(static_cast<evaluation_type>(k)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ELLINT_2_HPP
