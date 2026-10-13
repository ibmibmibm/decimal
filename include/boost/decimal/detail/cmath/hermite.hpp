// Copyright 2006 John Maddock
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_HERMITE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_HERMITE_HPP

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

template <typename T>
constexpr auto hermite_impl(const unsigned n, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(x))
    {
        return x;
    }
    if (isinf(x))
    {
        // (2x)^n leads
        return n == 0U ? T {1} : (signbit(x) && (n & 1U) == 0U ? -x : x);
    }
    #endif

    if (n == 0U)
    {
        return T {1};
    }
    // H_1(x) = 2x takes one rounding, and a zero x gives +0 as the recurrence does
    if (n == 1U)
    {
        return x == 0 ? T {0} : x + x;
    }
    return poly_detail::poly_fixed(poly_detail::poly_kind::hermite, n, 0U, x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto hermite(const unsigned n, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::hermite_impl(n, static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_HERMITE_HPP
