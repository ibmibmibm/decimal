// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <limits>

using namespace boost::decimal;

// Zero times infinity is invalid (IEEE 754 7.2), as when the zero is a decimal
template <typename T>
void test_infinity_times_zero()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const T inf {std::numeric_limits<T>::infinity()};

    BOOST_TEST(isnan(inf * 0));
    BOOST_TEST(isnan(0 * inf));
    BOOST_TEST(isnan(-inf * 0));
    BOOST_TEST(isnan(0U * -inf));
    #endif
}

// The product has the preferred exponent, the sum of the two exponents (IEEE 754 5.4.1)
template <typename T>
void test_preferred_exponent()
{
    BOOST_TEST_EQ(quantexp(T{3} * 5), 0);
    BOOST_TEST_EQ(quantexp(5 * T{3}), 0);
    BOOST_TEST_EQ(quantexp(T{15, -1} * 4), -1);
    BOOST_TEST_EQ(quantexp(T{15, -1} * -4LL), -1);
}

int main()
{
    test_infinity_times_zero<decimal32_t>();
    test_infinity_times_zero<decimal64_t>();
    test_infinity_times_zero<decimal128_t>();
    test_infinity_times_zero<decimal_fast32_t>();
    test_infinity_times_zero<decimal_fast64_t>();
    test_infinity_times_zero<decimal_fast128_t>();

    test_preferred_exponent<decimal32_t>();
    test_preferred_exponent<decimal64_t>();
    test_preferred_exponent<decimal128_t>();

    return boost::report_errors();
}
