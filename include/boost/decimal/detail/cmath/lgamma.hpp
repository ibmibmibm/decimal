// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_LGAMMA_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_LGAMMA_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/impl/gamma_impl.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto lgamma_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    T result { };

    const auto fpc = fpclassify(x);

    #ifndef BOOST_DECIMAL_FAST_MATH
    if ((fpc == FP_ZERO) || (fpc == FP_INFINITE))
    {
        result = std::numeric_limits<T>::infinity();
    }
    else if (fpc == FP_NAN)
    {
        result = x;
    }
    #else
    if (fpc == FP_ZERO)
    {
        result = std::numeric_limits<T>::max();
    }
    #endif
    else
    {
        int e {};
        const auto sig {frexp10(x, &e)};

        const bool neg_x {signbit(x)};

        if (neg_x && gamma_detail::gamma_is_integer(sig, e))
        {
            // Pure negative integer argument.
            #ifndef BOOST_DECIMAL_FAST_MATH
            result = std::numeric_limits<T>::infinity();
            #else
            result = T{0};
            #endif
        }
        else
        {
            result = gamma_detail::lgamma_fixed<T>(sig, e, neg_x);
        }
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto lgamma(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::lgamma_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_LGAMMA_HPP
