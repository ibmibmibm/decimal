// Copyright 2006 John Maddock
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_LAGUERRE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_LAGUERRE_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/impl/polynomial_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The Laguerre polynomial of degree n, and of order m for assoc_laguerre
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
constexpr auto laguerre_impl(const unsigned n, const T x, const unsigned m = 0U)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(x))
    {
        return x;
    }
    if (isinf(x))
    {
        // (-x)^n / n! leads
        const T a {signbit(x) ? -x : x};
        return n == 0U ? T {1} : (signbit(x) || (n & 1U) == 0U ? a : -a);
    }
    #endif

    if (n == 0U || (m == 0U && x == 0))
    {
        return T {1};
    }
    // L_1^m(x) = m + 1 - x takes one rounding where m + 1 has 7 digits at most, and a zero result is +0 as from the recurrence
    if (n == 1U && m < 9999999U)
    {
        const T r {T {m + 1U} - x};
        return r == 0 ? T {0} : r;
    }
    return poly_detail::poly_fixed(poly_detail::poly_kind::laguerre, n, m, x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto laguerre(const unsigned n, const T x)
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::laguerre_impl(n, static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_LAGUERRE_HPP
