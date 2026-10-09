// Copyright 2024 Matt Borland
// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP

#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/cmath/atan.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <cstdint>
#include <cmath>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto atan2_impl(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpcx {fpclassify(x)};
    const auto fpcy {fpclassify(y)};
    const auto signx {signbit(x)}; // True if neg
    const auto signy {signbit(y)};
    const auto isfinitex {fpcx != FP_INFINITE && fpcx != FP_NAN};
    const auto isfinitey {fpcy != FP_INFINITE && fpcy != FP_NAN};

    T result { };

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (fpcx == FP_NAN)
    {
        result = x;
    }
    else if (fpcy == FP_NAN)
    {
        result = y;
    }
    else
    #endif
    if (fpcy == FP_ZERO && signx)
    {
        result = atan_detail::pi_quarters<T>(4, signy);
    }
    else if (fpcy == FP_ZERO && !signx)
    {
        result = y;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpcy == FP_INFINITE && isfinitex)
    {
        result = atan_detail::pi_quarters<T>(2, signy);
    }
    else if (fpcy == FP_INFINITE && fpcx == FP_INFINITE && signx)
    {
        result = atan_detail::pi_quarters<T>(3, signy);
    }
    else if (fpcy == FP_INFINITE && fpcx == FP_INFINITE && !signx)
    {
        result = atan_detail::pi_quarters<T>(1, signy);
    }
    #endif
    else if (fpcx == FP_ZERO)
    {
        result = atan_detail::pi_quarters<T>(2, signy);
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpcx == FP_INFINITE && signx && isfinitey)
    {
        result = atan_detail::pi_quarters<T>(4, signy);
    }
    else if (fpcx == FP_INFINITE && !signx && isfinitey)
    {
        constexpr T zero { 0, 0 };

        result = signy ? -zero : zero;
    }
    #endif
    else
    {
        result = atan_detail::atan2_fixed(y, x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto atan2(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::atan2_impl(static_cast<evaluation_type>(y), static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP
