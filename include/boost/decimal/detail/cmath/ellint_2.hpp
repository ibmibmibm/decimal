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

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto ellint_2(const T k, const T phi) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::ellint_detail::ellint_checked(static_cast<evaluation_type>(k), static_cast<evaluation_type>(phi), detail::ellint_detail::ellint_kind::second));
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto comp_ellint_2(const T k) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::ellint_detail::comp_ellint_checked(static_cast<evaluation_type>(k), detail::ellint_detail::ellint_kind::second));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ELLINT_2_HPP
