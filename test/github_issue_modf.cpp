// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename T>
bool is_neg_zero(const T x)
{
    return x == 0 && signbit(x);
}

// The fractional part has the sign of x, also when it is zero
template <typename T>
void test_sign()
{
    T ip {};

    BOOST_TEST(is_neg_zero(modf(T(3U, 0, true), &ip)));
    BOOST_TEST_EQ(ip, T(-3));
    BOOST_TEST(is_neg_zero(modf(T(30U, -1, true), &ip)));
    BOOST_TEST_EQ(ip, T(-3));
    BOOST_TEST(is_neg_zero(modf(T(5U, 3, true), &ip)));
    BOOST_TEST_EQ(ip, T(-5000));
    BOOST_TEST(is_neg_zero(modf(T(0U, 0, true), &ip)));
    BOOST_TEST(is_neg_zero(ip));
    BOOST_TEST_EQ(modf(T(375U, -2, true), &ip), T(75U, -2, true));
    BOOST_TEST_EQ(ip, T(-3));
    BOOST_TEST_EQ(modf(T(25U, -2, true), &ip), T(25U, -2, true));
    BOOST_TEST(is_neg_zero(ip));
}

// The integer part has the exponent max(q(x), 0), and the fractional part keeps q(x)
void test_quantum()
{
    decimal64_t ip {};

    BOOST_TEST_EQ(quantexp(modf("0.25"_DD, &ip)), -2);
    BOOST_TEST_EQ(quantexp(ip), 0);
    BOOST_TEST_EQ(quantexp(modf("3.750"_DD, &ip)), -3);
    BOOST_TEST_EQ(quantexp(ip), 0);
    BOOST_TEST_EQ(quantexp(modf("5e3"_DD, &ip)), 3);
    BOOST_TEST_EQ(quantexp(ip), 3);
    BOOST_TEST_EQ(quantexp(modf("0e-5"_DD, &ip)), -5);
    BOOST_TEST_EQ(quantexp(ip), 0);

    decimal32_t ip32 {};
    BOOST_TEST_EQ(quantexp(modf("-0.0025"_DF, &ip32)), -4);
    BOOST_TEST_EQ(quantexp(ip32), 0);

    decimal128_t ip128 {};
    BOOST_TEST_EQ(modf("1234567890123456789012345678901234e-10"_DL, &ip128), "0.5678901234"_DL);
    BOOST_TEST_EQ(ip128, "123456789012345678901234"_DL);
    BOOST_TEST_EQ(quantexp(modf("1.5e-6000"_DL, &ip128)), -6001);
    BOOST_TEST_EQ(quantexp(ip128), 0);
}

// from_bits can make a significand with precision + 1 digits, here 10000001e-7 and 10000000000000001e-16
void test_non_canonical()
{
    decimal32_t ip32 {};
    BOOST_TEST_EQ(modf(from_bits(UINT32_C(0x6BD89681)), &ip32), "1e-7"_DF);
    BOOST_TEST_EQ(ip32, 1);

    decimal64_t ip64 {};
    BOOST_TEST_EQ(modf(from_bits(UINT64_C(0x6BF386F26FC10001)), &ip64), "1e-16"_DD);
    BOOST_TEST_EQ(ip64, 1);
}

int main()
{
    test_sign<decimal32_t>();
    test_sign<decimal64_t>();
    test_sign<decimal128_t>();
    test_sign<decimal_fast32_t>();
    test_sign<decimal_fast64_t>();
    test_sign<decimal_fast128_t>();

    test_quantum();
    test_non_canonical();

    return boost::report_errors();
}
