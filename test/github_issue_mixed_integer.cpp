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

// An integer operand gives the same result as the operand converted to T first
template <typename T, typename Integer>
void test_operand(const T x, const Integer n)
{
    const T dn {n};

    BOOST_TEST_EQ(x + n, x + dn);
    BOOST_TEST_EQ(n + x, dn + x);
    BOOST_TEST_EQ(x - n, x - dn);
    BOOST_TEST_EQ(n - x, dn - x);
    BOOST_TEST_EQ(x * n, x * dn);
    BOOST_TEST_EQ(n * x, dn * x);
    BOOST_TEST_EQ(x / n, x / dn);
    BOOST_TEST_EQ(n / x, dn / x);
}

// big is 10^(p + 1), so each operand has more digits than T and the conversion rounds it
template <typename T, typename Integer>
void test_wide(const Integer big)
{
    const T x {3, -1};

    test_operand(x, big + 5);
    test_operand(x, -(big + 5));
    test_operand(x, big * 9 + 15);
    test_operand(x, -(big / 10 + 5));
}

template <typename T, typename Integer>
void test_all_modes(const Integer big)
{
    test_wide<T>(big);

    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    for (const auto mode : {rounding_mode::fe_dec_downward, rounding_mode::fe_dec_upward,
                            rounding_mode::fe_dec_toward_zero, rounding_mode::fe_dec_to_nearest_from_zero})
    {
        fesetround(mode);
        test_wide<T>(big);
    }
    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

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

    test_all_modes<decimal32_t>(INT32_C(100000000));
    test_all_modes<decimal_fast32_t>(INT32_C(100000000));
    test_all_modes<decimal64_t>(INT64_C(100000000000000000));
    test_all_modes<decimal_fast64_t>(INT64_C(100000000000000000));

    const auto big128 {boost::int128::int128_t{INT64_C(10000000000000000)} * INT64_C(1000000000000000000) * 10};
    test_all_modes<decimal128_t>(big128);
    test_all_modes<decimal_fast128_t>(big128);

    return boost::report_errors();
}
