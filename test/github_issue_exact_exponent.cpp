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

// A product of at most p digits is exact and has the sum of the exponents
template <typename T>
void test_exact_product()
{
    BOOST_TEST_EQ(quantexp(T{1234567} * T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{1234567, 3} * T{1, -5}), -2);
    BOOST_TEST_EQ(quantexp(T{1234567, 3} * T{0, -5}), -2);
    BOOST_TEST_EQ(quantexp(T{0, 4} * T{0, -3}), 1);
}

void test_exact_product_128()
{
    const decimal128_t x {123456789012345678LL};

    BOOST_TEST_EQ(quantexp(x * 1), 0);
    BOOST_TEST_EQ(quantexp(1 * x), 0);
    BOOST_TEST_EQ(quantexp(x * decimal128_t{1}), 0);
    BOOST_TEST_EQ(quantexp(x * decimal128_t{1000000000000000LL, -2}), -2);
    BOOST_TEST_EQ(quantexp(decimal128_t{1, 20} * decimal128_t{0, -3}), 17);
    BOOST_TEST_EQ(quantexp(decimal128_t{0, -3} * x), -3);

    BOOST_TEST_EQ(x * 1, x);
    BOOST_TEST_EQ(x * decimal128_t(1000000000000000LL, -2), decimal128_t(123456789012345678LL, 13));

    // 34 digits fit, and 35 digits round to 34
    const decimal128_t nines {9999999999999999LL};
    const decimal128_t wide {999999999999999999LL};
    BOOST_TEST_EQ(quantexp(nines * wide), 0);
    BOOST_TEST_EQ(quantexp(wide * (nines * 10 + 9)), 1);
    BOOST_TEST_EQ(nines * wide, 9999999999999998990000000000000001_DL);
}

// An exact quotient has the difference of the exponents, or the nearest exponent
// that holds its digits. An inexact quotient has all the digits.
template <typename T>
void test_exact_quotient()
{
    constexpr int digits {std::numeric_limits<T>::digits10};

    BOOST_TEST_EQ(quantexp(T{1} / T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{10} / T{2}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 3} / T{1, -2}), 5);
    BOOST_TEST_EQ(quantexp(T{15, -1} / T{5}), -1);
    BOOST_TEST_EQ(quantexp(T{1} / T{4}), -2);
    BOOST_TEST_EQ(quantexp(T{1, 2} / T{8, 5}), -6);
    BOOST_TEST_EQ(quantexp(T{6} / 2), 0);
    BOOST_TEST_EQ(quantexp(6 / T{2}), 0);
    BOOST_TEST_EQ(quantexp(T{15, -1} / -5), -1);
    BOOST_TEST_EQ(quantexp(1 / T{4, 1}), -3);
    BOOST_TEST_EQ(quantexp(T{1} / T{3}), -digits);
    BOOST_TEST_EQ(quantexp(T{2} / 3), -digits);

    BOOST_TEST_EQ(quantexp(T{0, 3} / T{1, -2}), 5);
    BOOST_TEST_EQ(quantexp(T{0, -2} / 7), -2);
    BOOST_TEST_EQ(quantexp(0 / T{25, -1}), 1);
    BOOST_TEST(signbit(T{0, 3} / T{-1, -2}));

    BOOST_TEST_EQ(T(1) / T(4), T(25, -2));
    BOOST_TEST_EQ(T(1, 3) / T(1, -2), T(1, 5));
    BOOST_TEST_EQ(T(15, -1) / -5, T(-3, -1));
}

int main()
{
    test_zero_operand<decimal32_t>();
    test_zero_operand<decimal64_t>();
    test_zero_operand<decimal128_t>();

    test_short_sum<decimal32_t>();
    test_short_sum<decimal64_t>();
    test_short_sum<decimal128_t>();

    test_exact_product<decimal32_t>();
    test_exact_product<decimal64_t>();
    test_exact_product<decimal128_t>();
    test_exact_product_128();

    test_exact_quotient<decimal32_t>();
    test_exact_quotient<decimal64_t>();
    test_exact_quotient<decimal128_t>();

    return boost::report_errors();
}
