// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

auto second_significand_64() -> const std::uint64_t*;
auto second_significand_128() -> const std::uint64_t*;
auto second_coarse_power_high() -> const std::uint64_t*;
auto second_coarse_power_low() -> const std::uint64_t*;
auto second_coarse_power_exponent() -> const std::int32_t*;
auto second_clinger_powers() -> const double*;
auto second_long_double_clinger_powers() -> const long double*;
auto second_ryu_inv_pow5() -> const boost::int128::uint128_t*;

// Each source file must see the same table, also in C++14
int main()
{
    using boost::decimal::detail::fast_float::significand_table;
    using boost::decimal::detail::fast_float::coarse_power_table;
    using boost::decimal::detail::fast_float::clinger_table;
    using boost::decimal::detail::fast_float::long_double_clinger_table;

    BOOST_TEST(significand_table::significand_64 == second_significand_64());
    BOOST_TEST(significand_table::significand_128 == second_significand_128());
    BOOST_TEST(coarse_power_table::high == second_coarse_power_high());
    BOOST_TEST(coarse_power_table::low == second_coarse_power_low());
    BOOST_TEST(coarse_power_table::exponent == second_coarse_power_exponent());
    BOOST_TEST(clinger_table::powers_of_ten == second_clinger_powers());
    BOOST_TEST(long_double_clinger_table::powers_of_ten == second_long_double_clinger_powers());
    BOOST_TEST(boost::decimal::detail::ryu::inv_pow5_table::values == second_ryu_inv_pow5());
    BOOST_TEST_EQ(static_cast<boost::decimal::decimal32_t>(0.5), boost::decimal::decimal32_t(5, -1));

    return boost::report_errors();
}
