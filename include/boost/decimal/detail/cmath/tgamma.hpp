// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_TGAMMA_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_TGAMMA_HPP

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
constexpr auto tgamma_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    T result { };

    const bool is_neg = signbit(x);

    const auto fpc = fpclassify(x);

    if (fpc == FP_ZERO)
    {
        result = (is_neg ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity());
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpc == FP_INFINITE)
    {
        result = (is_neg ? std::numeric_limits<T>::quiet_NaN() : std::numeric_limits<T>::infinity());
    }
    else if (fpc == FP_NAN)
    {
        result = x;
    }
    #endif
    else
    {
        int e {};
        const auto sig {frexp10(x, &e)};

        if (is_neg && gamma_detail::gamma_is_integer(sig, e))
        {
            // Pure negative integer argument.
            #ifndef BOOST_DECIMAL_FAST_MATH
            result = std::numeric_limits<T>::quiet_NaN();
            #else
            result = T{0};
            #endif
        }
        else
        {
            result = gamma_detail::tgamma_fixed<T>(sig, e, is_neg);
        }
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto tgamma(const T x) noexcept                           // LCOV_EXCL_LINE
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)  // LCOV_EXCL_LINE
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::tgamma_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_TGAMMA_HPP
