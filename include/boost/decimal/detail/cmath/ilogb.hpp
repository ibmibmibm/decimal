// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ILOGB_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ILOGB_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cmath>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The exponent of the leading digit of a finite nonzero value, from its exponent and its nonzero significand.
// A subnormal value needs no other path.
template <typename Sig>
constexpr auto leading_digit_exponent(const int exp, const Sig sig) noexcept -> int
{
    return exp + num_digits(sig) - 1;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto ilogb(const T d) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_decimal_floating_point_v, T, int)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(d))
    {
        return isnan(d) ? FP_ILOGBNAN : INT_MAX;
    }
    #endif

    const auto sig {d.full_significand()};
    return sig == 0U ? FP_ILOGB0 : detail::leading_digit_exponent(static_cast<int>(d.biased_exponent()), sig);
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_ILOGB_HPP
