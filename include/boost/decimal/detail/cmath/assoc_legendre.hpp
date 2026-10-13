//  (C) Copyright John Maddock 2006.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LEGENDRE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LEGENDRE_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/cmath/legendre.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
constexpr auto assoc_legendre_impl(const unsigned l, const unsigned m, const T x) noexcept
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
        return issignaling(x) ? nan_conversion(x) : x;
    }
    #endif

    if (m == 0U)
    {
        return legendre_impl(l, x);
    }
    if (m > l || x == 1 || x == -1)
    {
        return T {0};
    }
    return poly_detail::poly_fixed<T, true>(poly_detail::poly_kind::legendre, l - m, m, x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto assoc_legendre(const unsigned n, const unsigned m, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::assoc_legendre_impl(n, m, static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LEGENDRE_HPP
