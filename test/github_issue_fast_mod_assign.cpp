// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal::literals;

int main()
{
    auto fast64 {7.5_DDF};
    auto& fast64_result {fast64 %= 2_DDF};
    BOOST_TEST_EQ(fast64, 1.5_DDF);
    BOOST_TEST_EQ(&fast64_result, &fast64);

    auto fast128 {7.5_DLF};
    auto& fast128_result {fast128 %= 2_DLF};
    BOOST_TEST_EQ(fast128, 1.5_DLF);
    BOOST_TEST_EQ(&fast128_result, &fast128);

    auto self {7.5_DLF};
    self %= self;
    BOOST_TEST_EQ(self, 0_DLF);

    return boost::report_errors();
}
