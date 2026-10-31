// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_LOGB_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_LOGB_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/ilogb.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cmath>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto logb(const T num) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(num))
    {
        return isnan(num) ? num : std::numeric_limits<T>::infinity();
    }
    #endif

    const auto sig {num.full_significand()};
    if (sig == 0U)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return -std::numeric_limits<T>::infinity();
        #else
        return T{0};
        #endif
    }

    return static_cast<T>(detail::leading_digit_exponent(static_cast<int>(num.biased_exponent()), sig));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_LOGB_HPP
