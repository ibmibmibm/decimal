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

    // The middle words of the product carry, and the product has 39 digits
    BOOST_TEST_EQ(quantexp(36893488147419103231_DL * 9223372036854775809_DL), 5);
    BOOST_TEST_EQ(36893488147419103231_DL * 9223372036854775809_DL, 3402823669209384634910447235423325e5_DL);
}

// An exact quotient has the difference of the exponents, or the nearest exponent
// that holds its digits. An inexact quotient has all the digits.
template <typename T>
void test_exact_quotient()
{
    constexpr int digits {std::numeric_limits<T>::digits10};

    BOOST_TEST_EQ(quantexp(T{1} / T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{10} / T{2}), 0);
    BOOST_TEST_EQ(quantexp(T{100} / T{1}), 0);
    BOOST_TEST_EQ(quantexp(100 / T{1}), 0);
    BOOST_TEST_EQ(quantexp(T{1, 3} / T{1, -2}), 5);
    BOOST_TEST_EQ(quantexp(T{15, -1} / T{5}), -1);
    BOOST_TEST_EQ(quantexp(T{1} / T{2}), -1);
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

    // A preferred exponent out of range gives the nearest exponent in range
    const int max_exp {quantexp((std::numeric_limits<T>::max)())};
    const int min_exp {quantexp(std::numeric_limits<T>::denorm_min())};
    BOOST_TEST_EQ(quantexp(T{1, max_exp - 5} / T{1, -10}), max_exp);
    BOOST_TEST_EQ((T{1, max_exp - 5} / T{1, -10}), T(1, max_exp + 5));
    BOOST_TEST_EQ(quantexp(T{1000000, min_exp} / T{1, 2}), min_exp);
    BOOST_TEST_EQ((T{1000000, min_exp} / T{1, 2}), T(1, min_exp + 4));
}

#if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)

// Each entry of the table times 5^k is 1 modulo 2^128
void test_inverse_table()
{
    for (unsigned k {0}; k < 34U; ++k)
    {
        const auto five_to_k {boost::decimal::detail::pow10(boost::int128::uint128_t{k}) >> k};
        BOOST_TEST(boost::decimal::detail::impl::inverse_pow5_u128[k] * five_to_k == 1U);
    }
}

#endif

// The rules use an integer operand after its conversion, which rounds a wide integer
void test_wide_integer_operand()
{
    BOOST_TEST_EQ(quantexp(decimal32_t{1234567, 1} - 12345670), 1);
    BOOST_TEST_EQ(quantexp(decimal32_t{0, 5} * 12345670), 6);
    BOOST_TEST_EQ(quantexp(decimal64_t{1234567890123456, 1} - 12345678901234560), 1);
    BOOST_TEST_EQ(quantexp(decimal64_t{0, 5} * 12345678901234560), 6);
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
    test_wide_integer_operand();
    #if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)
    test_inverse_table();
    #endif

    return boost::report_errors();
}
