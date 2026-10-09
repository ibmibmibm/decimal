// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_COS_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_COS_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/cmath/impl/trig_reduce.hpp>
#include <boost/decimal/detail/cmath/impl/sin_impl.hpp>
#include <boost/decimal/detail/cmath/impl/cos_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto cos_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(x))
    {
        return issignaling(x) ? nan_conversion(x) : x;
    }
    if (isinf(x))
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    #endif

    BOOST_DECIMAL_IF_CONSTEXPR (trig::trig_traits<T>::words > 1)
    {
        // A first try in one word less
        constexpr int words {trig::trig_traits<T>::words > 1 ? trig::trig_traits<T>::words - 1 : 1};
        trig::trig_arg<words> r1 {};
        T res {};
        if (trig::trig_prepare_try(x, r1) && trig::sin_cos_try(r1, (r1.n & 1U) == 0U, trig::cos_neg(r1), res))
        {
            return res;
        }
    }

    // x = n*pi/2 + r: cos(x) is cos(r), -sin(r), -cos(r) or sin(r) for n = 0 to 3.
    const auto r {trig::trig_prepare(x)};
    if (r.zero)
    {
        return T {1};
    }

    return trig::sin_cos_of<T>(r, (r.n & 1U) == 0U, trig::cos_neg(r));
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto cos(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::cos_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_COS_HPP
