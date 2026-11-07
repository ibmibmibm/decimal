// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <cstdint>

auto second_significand_64() -> const std::uint64_t*
{
    return boost::decimal::detail::fast_float::significand_table::significand_64;
}

auto second_significand_128() -> const std::uint64_t*
{
    return boost::decimal::detail::fast_float::significand_table::significand_128;
}

auto second_coarse_power_high() -> const std::uint64_t*
{
    return boost::decimal::detail::fast_float::coarse_power_table::high;
}

auto second_coarse_power_low() -> const std::uint64_t*
{
    return boost::decimal::detail::fast_float::coarse_power_table::low;
}

auto second_coarse_power_exponent() -> const std::int32_t*
{
    return boost::decimal::detail::fast_float::coarse_power_table::exponent;
}

auto second_clinger_powers() -> const double*
{
    return boost::decimal::detail::fast_float::clinger_table::powers_of_ten;
}

auto second_long_double_clinger_powers() -> const long double*
{
    return boost::decimal::detail::fast_float::long_double_clinger_table::powers_of_ten;
}
