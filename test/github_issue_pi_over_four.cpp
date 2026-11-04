// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

int main()
{
    // pi / 4 to 50 digits, which each literal rounds to the precision of its type
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal32_t>, 0.78539816339744830961566084581987572104929234984378_DF);
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal64_t>, 0.78539816339744830961566084581987572104929234984378_DD);
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal128_t>, 0.78539816339744830961566084581987572104929234984378_DL);
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal_fast32_t>, 0.78539816339744830961566084581987572104929234984378_DFF);
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal_fast64_t>, 0.78539816339744830961566084581987572104929234984378_DDF);
    BOOST_TEST_EQ(numbers::pi_over_four_v<decimal_fast128_t>, 0.78539816339744830961566084581987572104929234984378_DLF);

    return boost::report_errors();
}
