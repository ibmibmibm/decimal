// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// https://github.com/boostorg/decimal/issues/NNNN
//
// The u256 division and remainder counted the words of a zero dividend past the start of
// the array, thus fmod of a zero decimal_fast128_t read outside the stack frame.

#include <boost/decimal.hpp>
#include <boost/decimal/detail/u256.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The divisor has more than one word, thus the division does not take the 64-bit path.
void test_zero_dividend()
{
    const detail::u256 zero {UINT64_C(0)};
    const detail::u256 divisor {detail::pow10_256(30)};

    BOOST_TEST(zero / divisor == zero);
    BOOST_TEST(zero % divisor == zero);

    const auto result {detail::impl::div_mod(zero, divisor)};
    BOOST_TEST(result.quotient == zero);
    BOOST_TEST(result.remainder == zero);
}

// The fast type keeps a significand with 34 digits, thus 3 is above 2^64 in it.
void test_fmod_of_zero()
{
    BOOST_TEST_EQ(fmod(0_DLF, 3_DLF), 0_DLF);
    BOOST_TEST_EQ(fmod(-0_DLF, 3_DLF), -0_DLF);
    BOOST_TEST(signbit(fmod(-0_DLF, 3_DLF)));
    BOOST_TEST_EQ(fmod(0_DL, 1.234567890123456789012_DL), 0_DL);
    BOOST_TEST_EQ(remainder(0_DLF, 3_DLF), 0_DLF);
}

int main()
{
    test_zero_dividend();
    test_fmod_of_zero();

    return boost::report_errors();
}
