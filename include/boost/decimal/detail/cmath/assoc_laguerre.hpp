// Copyright 2006 John Maddock
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LAGUERRE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LAGUERRE_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/laguerre.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto assoc_laguerre_impl(const unsigned n, const unsigned m, const T x)
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    return laguerre_impl(n, x, m);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto assoc_laguerre(const unsigned n, const unsigned m, const T x)
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::assoc_laguerre_impl(n, m, static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ASSOC_LAGUERRE_HPP
