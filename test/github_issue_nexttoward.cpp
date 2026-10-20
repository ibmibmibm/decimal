// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// A direction that rounds to x in the decimal type is still above or below x
void test_direction_rounds_to_x()
{
    BOOST_TEST_EQ(nexttoward("1"_DF, 1.0L + 1e-12L), "1.000001"_DF);
    BOOST_TEST_EQ(nexttoward("1"_DF, 1.0L - 1e-12L), "0.9999999"_DF);
    BOOST_TEST_EQ(nexttoward("1"_DD, 1.0L + 4e-16L), "1.000000000000001"_DD);
    BOOST_TEST_EQ(nexttoward("1"_DFF, 1.0L - 1e-12L), "0.9999999"_DFF);

    // The binary value of 0.1L is above 0.1
    BOOST_TEST_EQ(nexttoward("0.1"_DF, 0.1L), "0.1000001"_DF);
    BOOST_TEST_EQ(nexttoward("-0.1"_DF, -0.1L), "-0.1000001"_DF);
}

// The shortest decimal string of a long double can be on the other side of x
void test_shortest_string_other_side()
{
    #if BOOST_DECIMAL_LDBL_BITS == 80
    // y = 1 - 2^-64 is 0.99999999999999999994578989137572477...
    const long double y {1.0L - 1.0L / 18446744073709551616.0L};
    BOOST_TEST_EQ(nexttoward("0.9999999999999999999457898913757248"_DL, y), "0.9999999999999999999457898913757247"_DL);
    BOOST_TEST_EQ(nexttoward("0.9999999999999999999457898913757247"_DL, y), "0.9999999999999999999457898913757248"_DL);
    #endif
}

// An equal direction gives x, with the sign of direction for a zero
void test_equal()
{
    BOOST_TEST_EQ(nexttoward("0.5"_DF, 0.5L), "0.5"_DF);

    const auto neg_zero {nexttoward("0"_DF, -0.0L)};
    BOOST_TEST_EQ(neg_zero, "0"_DF);
    BOOST_TEST(signbit(neg_zero));

    BOOST_TEST_EQ(nexttoward("0"_DD, 1e-300L), std::numeric_limits<decimal64_t>::denorm_min());
    BOOST_TEST_EQ(nexttoward(std::numeric_limits<decimal32_t>::denorm_min(), 0.0L), "0"_DF);
}

void test_non_finite()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto inf {std::numeric_limits<long double>::infinity()};

    BOOST_TEST(isinf(nexttoward((std::numeric_limits<decimal32_t>::max)(), inf)));
    BOOST_TEST_EQ(nexttoward(std::numeric_limits<decimal32_t>::infinity(), 1.0L), (std::numeric_limits<decimal32_t>::max)());
    BOOST_TEST(isinf(nexttoward(std::numeric_limits<decimal32_t>::infinity(), inf)));
    BOOST_TEST_EQ(nexttoward("1"_DF, -inf), "0.9999999"_DF);
    BOOST_TEST(isnan(nexttoward("1"_DF, std::numeric_limits<long double>::quiet_NaN())));

    // A decimal128_t value can be beyond the long double range
    BOOST_TEST_EQ(nexttoward("1e5000"_DL, inf), "1.000000000000000000000000000000001e5000"_DL);
    BOOST_TEST_EQ(nexttoward("1e5000"_DL, 1.0L), "9.999999999999999999999999999999999e4999"_DL);
    #endif
}

int main()
{
    test_direction_rounds_to_x();
    test_shortest_string_other_side();
    test_equal();
    test_non_finite();

    return boost::report_errors();
}
