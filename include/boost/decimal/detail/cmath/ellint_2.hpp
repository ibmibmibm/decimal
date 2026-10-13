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
    if (isnan(k) || isnan(phi) || fabs(k) > T {1})
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    // The integral grows without bound in phi
    if (isinf(phi))
    {
        return phi;
    }
    #else
    if (fabs(k) > T {1})
    {
        return T {0};
    }
    #endif

    if (k == T {0} || phi == T {0})
    {
        return phi;
    }
    return ellint_detail::ellint_fixed(fabs(k), fabs(phi), true, signbit(phi));
}

template <typename T>
constexpr auto comp_ellint_2_impl(const T k) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(k) || fabs(k) > T {1})
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    #else
    if (fabs(k) > T {1})
    {
        return T {0};
    }
    #endif

    // E(1) = 1
    if (fabs(k) == T {1})
    {
        return T {1};
    }
    return ellint_detail::comp_ellint_fixed(fabs(k), true);
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
