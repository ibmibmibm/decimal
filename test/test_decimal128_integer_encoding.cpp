// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The integer constructor of decimal128_t gives the encoding of exponent 0

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <limits>

using namespace boost::decimal;

template <typename Int>
void test(const Int value)
{
    const decimal128_t a {value};
    const decimal128_t b {value, 0};
    BOOST_TEST(to_bid(a) == to_bid(b));
}

int main()
{
    test(0);
    test(-7);
    test((std::numeric_limits<std::int64_t>::min)());
    test((std::numeric_limits<std::uint64_t>::max)());

    // A value wider than 64 bits takes the general path of the constructor
    test(boost::int128::uint128_t {1U, 0U});

    return boost::report_errors();
}
