// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// tiny is too far below big to align with it, so only the rounding mode decides the result.
// above and below are the neighbors of big.
template <typename T>
void test_far_operands(const T big, const T tiny, const T above, const T below)
{
    for (const auto mode : {rounding_mode::fe_dec_to_nearest, rounding_mode::fe_dec_to_nearest_from_zero})
    {
        fesetround(mode);
        BOOST_TEST_EQ(big + tiny, big);
        BOOST_TEST_EQ(big - tiny, big);
        BOOST_TEST_EQ(tiny + big, big);
        BOOST_TEST_EQ(tiny - big, -big);
        BOOST_TEST_EQ(fma(big, T{1}, tiny), big);
        BOOST_TEST_EQ(fma(big, T{1}, -tiny), big);
    }

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(fma(big, T{1}, tiny), above);
    BOOST_TEST_EQ(fma(big, T{1}, -tiny), big);
    BOOST_TEST_EQ(fma(-big, T{1}, tiny), -below);
    BOOST_TEST_EQ(fma(tiny, T{1}, big), above);
    BOOST_TEST_EQ(big + tiny, above);
    BOOST_TEST_EQ(big - tiny, big);
    BOOST_TEST_EQ(-big + tiny, -below);
    BOOST_TEST_EQ(-big - tiny, -big);
    BOOST_TEST_EQ(tiny + big, above);
    BOOST_TEST_EQ(tiny - big, -below);

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(big + tiny, big);
    BOOST_TEST_EQ(big - tiny, below);
    BOOST_TEST_EQ(-big + tiny, -big);
    BOOST_TEST_EQ(-big - tiny, -above);
    BOOST_TEST_EQ(tiny + big, big);
    BOOST_TEST_EQ(tiny - big, -big);
    BOOST_TEST_EQ(fma(big, T{1}, tiny), big);
    BOOST_TEST_EQ(fma(big, T{1}, -tiny), below);
    BOOST_TEST_EQ(fma(-big, T{1}, -tiny), -above);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(fma(big, T{1}, -tiny), below);
    BOOST_TEST_EQ(fma(-big, T{1}, tiny), -below);
    BOOST_TEST_EQ(big + tiny, big);
    BOOST_TEST_EQ(big - tiny, below);
    BOOST_TEST_EQ(-big + tiny, -below);
    BOOST_TEST_EQ(-big - tiny, -big);

    fesetround(rounding_mode::fe_dec_default);
}

// A subnormal result can not use the kernels for a normal result
void test_subnormal_decimal128()
{
    BOOST_TEST_EQ(1e-6170_DL + 1e-6176_DL, 1000001e-6176_DL);
    BOOST_TEST_EQ(1e-6176_DL + 1e-6170_DL, 1000001e-6176_DL);
    BOOST_TEST_EQ(1e-6170_DL - 1e-6176_DL, 999999e-6176_DL);
    BOOST_TEST_EQ(1e-6176_DL - 1e-6170_DL, -999999e-6176_DL);
    BOOST_TEST_EQ(-1e-6170_DL + 1e-6176_DL, -999999e-6176_DL);
    BOOST_TEST_EQ(-1e-6170_DL - 1e-6176_DL, -1000001e-6176_DL);

    BOOST_TEST(!signbit(1e-6170_DL - 1000000e-6176_DL));
    BOOST_TEST_EQ(1e-6170_DL - 1000000e-6176_DL, 0_DL);

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST(signbit(1e-6170_DL - 1000000e-6176_DL));
    fesetround(rounding_mode::fe_dec_default);
}

int main()
{
    test_far_operands(1e30_DF, 1e-30_DF, 1.000001e30_DF, 9.999999e29_DF);
    test_far_operands(1.5e30_DF, 1e-30_DF, 1.500001e30_DF, 1.499999e30_DF);
    test_far_operands(1e30_DFF, 1e-30_DFF, 1.000001e30_DFF, 9.999999e29_DFF);

    test_far_operands(1e100_DD, 1e-100_DD, 1.000000000000001e100_DD, 9.999999999999999e99_DD);
    test_far_operands(1.5e100_DD, 1e-100_DD, 1.500000000000001e100_DD, 1.499999999999999e100_DD);
    test_far_operands(1e100_DDF, 1e-100_DDF, 1.000000000000001e100_DDF, 9.999999999999999e99_DDF);

    test_far_operands(1e1000_DL, 1e-1000_DL, 1.000000000000000000000000000000001e1000_DL, 9.999999999999999999999999999999999e999_DL);
    test_far_operands(1.5e1000_DL, 1e-1000_DL, 1.500000000000000000000000000000001e1000_DL, 1.499999999999999999999999999999999e1000_DL);
    test_far_operands(1e1000_DLF, 1e-1000_DLF, 1.000000000000000000000000000000001e1000_DLF, 9.999999999999999999999999999999999e999_DLF);

    test_subnormal_decimal128();

    return boost::report_errors();
}
