// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT32_HPP
#define BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT32_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/fast_float/compute_float64.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <cstdint>
#include <cmath>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace fast_float {

// i * 10^power, negated if negative
BOOST_DECIMAL_CXX20_CONSTEXPR auto compute_float32(const std::int64_t power, const std::uint64_t i, const bool negative, bool& success) noexcept -> float
{
    success = true;
    return to_binary<float>(compute_float<float>(power, i), negative);
}

} // namespace fast_float
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT32_HPP
