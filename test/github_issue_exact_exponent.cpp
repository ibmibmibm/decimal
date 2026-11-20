// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// A sum with a zero operand has the smaller exponent, or the nearest one that the
// digits of the other operand allow (IEEE 754 5.4.1)
template <typename T>
void test_zero_operand()
{
    constexpr int digits {std::numeric_limits<T>::digits10};

    BOOST_TEST_EQ(quantexp(T{0, -2} + T{0}), -2);
    BOOST_TEST_EQ(quantexp(T{0, -2} + 0), -2);
    BOOST_TEST_EQ(quantexp(T{0, -2} - T{0, 3}), -2);
    BOOST_TEST_EQ(quantexp(T{0, 3} + T{0, 5}), 3);

    BOOST_TEST_EQ(quantexp(T{1, 5} + T{0}), 0);
    BOOST_TEST_EQ(quantexp(T{0} + T{1, 5}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 5} - T{0}), 0);
    BOOST_TEST_EQ(quantexp(T{0} - T{1, 5}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 5} + 0), 0);
    BOOST_TEST_EQ(quantexp(0 - T{1, 5}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 3} + T{0, -2}), -2);

    // The zero has the larger exponent, or the other operand has too few places
    BOOST_TEST_EQ(quantexp(T{15, -1} + T{0, 3}), -1);
    BOOST_TEST_EQ(quantexp(T{0, 3} - T{15, -1}), -1);
    BOOST_TEST_EQ(quantexp(T{1, 40} + T{0}), 41 - digits);
    BOOST_TEST_EQ(quantexp(T{0} - T{12, 40}), 42 - digits);

    BOOST_TEST_EQ(T(1, 5) + T(0), T(1, 5));
    BOOST_TEST_EQ(T(0) - T(1, 5), T(-1, 5));
    BOOST_TEST_EQ(T(0, 3) - T(15, -1), T(-15, -1));

    BOOST_TEST(!signbit(T{0} - T{0, 2}));
    BOOST_TEST(!signbit(T{0, 2} + -T{0}));
    BOOST_TEST(signbit(-T{0, 2} + -T{0}));
    BOOST_TEST(signbit(-T{0} - T{0, 2}));
}

// An exact sum of short operands has the smaller exponent at each exponent gap
template <typename T>
void test_short_sum()
{
    BOOST_TEST_EQ(quantexp(T{1, 4} + 1), 0);
    BOOST_TEST_EQ(quantexp(1 + T{1, 4}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 4} + T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 4} - T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{1} - T{1, 4}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 5} + T{25, -1}), -1);
    BOOST_TEST_EQ(quantexp(T{1, 6} + T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{12, 2} - T{5, -3}), -3);

    BOOST_TEST_EQ(T(1, 4) + T(1), T(10001));
    BOOST_TEST_EQ(T(1) - T(1, 4), T(-9999));
    BOOST_TEST_EQ(T(1, 5) + T(25, -1), T(1000025, -1));
    BOOST_TEST_EQ(T(12, 2) - T(5, -3), T(1199995, -3));
}

int main()
{
    test_zero_operand<decimal32_t>();
    test_zero_operand<decimal64_t>();
    test_zero_operand<decimal128_t>();

    test_short_sum<decimal32_t>();
    test_short_sum<decimal64_t>();
    test_short_sum<decimal128_t>();

    return boost::report_errors();
}
