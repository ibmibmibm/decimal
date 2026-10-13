//  (C) Copyright John Maddock 2006.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_DECIMAL_DETAIL_CMATH_LEGENDRE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_LEGENDRE_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/impl/polynomial_impl.hpp>
#include <boost/decimal/detail/promotion.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The Legendre polynomial P_l(x)
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
constexpr auto legendre_impl(const unsigned l, const T x) noexcept
{
    if (x < -1 || x > 1 || l > 128)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T{0};
        #endif
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (isnan(x))
    {
        return x;
    }
    #endif

    // P_0 = 1, P_l(1) = 1 and P_l(-1) = (-1)^l, and P_1(x) = x with +0 for a zero x as from the recurrence
    if (l == 0U || x == 1)
    {
        return T {1};
    }
    if (x == -1)
    {
        return (l & 1U) != 0U ? T {-1} : T {1};
    }
    if (l == 1U)
    {
        return x == 0 ? T {0} : x;
    }
    return poly_detail::poly_fixed(poly_detail::poly_kind::legendre, l, 0U, x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto legendre(const unsigned n, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::legendre_impl(n, static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_LEGENDRE_HPP
