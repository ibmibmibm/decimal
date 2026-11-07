// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// Zero divided by zero is NaN (IEEE 754 7.2), also when one zero is an integer
template <typename T, typename Integer>
void test_zero(const T zero, const Integer n)
{
    BOOST_TEST(isnan(zero / n));
    BOOST_TEST(isnan(n / zero));
    BOOST_TEST(!issignaling(zero / n));
    BOOST_TEST(!issignaling(n / zero));

    T x {zero};
    x /= n;
    BOOST_TEST(isnan(x));
}

template <typename T>
void test()
{
    const T zero {0};

    #ifndef BOOST_DECIMAL_FAST_MATH
    test_zero(zero, 0);
    test_zero(-zero, 0);
    test_zero(zero, 0U);
    test_zero(-zero, 0LL);
    test_zero(zero, 0ULL);

    // A non zero integer divided by a zero, and the reverse, do not change
    BOOST_TEST(isinf(5 / zero));
    BOOST_TEST(signbit(-5 / zero));
    BOOST_TEST(isinf(T{1} / 0));

    // An infinity divided by an integer has the sign of the quotient
    const T inf {std::numeric_limits<T>::infinity()};
    BOOST_TEST_EQ(-inf / 2, -inf);
    BOOST_TEST_EQ(-inf / -2, inf);
    BOOST_TEST_EQ(inf / -2, -inf);
    BOOST_TEST_EQ(-inf / 2U, -inf);
    #endif

    BOOST_TEST_EQ(zero / 5, zero);
    BOOST_TEST(signbit(zero / -5));
}

int main()
{
    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    return boost::report_errors();
}
