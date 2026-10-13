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

// The expected values are the integrals rounded once, from MPFR
void test_nearest()
{
    // phi below pi/4, between pi/4 and pi/2, and above, where F(phi) = 2K/pi phi + (F(psi) - 2K/pi psi)
    BOOST_TEST_EQ(ellint_1(0.5_DF, 0.5_DF), 0.5050887_DF);
    BOOST_TEST_EQ(ellint_2(0.8_DL, 0.7_DL), "0.6653218748449987907012271726281617"_DL);
    BOOST_TEST_EQ(ellint_1(0.9_DD, 1.2_DD), 1.495205431757420_DD);
    BOOST_TEST_EQ(ellint_2(0.3_DF, 1.5_DF), 1.467292_DF);
    BOOST_TEST_EQ(ellint_1(0.7_DL, 5_DL), "5.934746266872147354936963225687505"_DL);
    BOOST_TEST_EQ(ellint_2(0.4_DD, 10_DD), 9.606086712414683_DD);
    BOOST_TEST_EQ(ellint_1(0.6_DF, -3_DF), -3.359745_DF);
    BOOST_TEST_EQ(ellint_2(-0.6_DF, 3_DF), 2.694744_DF);
    BOOST_TEST_EQ(ellint_1(0.5_DL, 1E+40_DL), "1.073182007149364375052841707970350E+40"_DL);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 1E+80_DD), 1.073182007149364E+80_DD);
    BOOST_TEST_EQ(ellint_2(0.5_DF, 1E+50_DF), 9.342155E+49_DF);
    BOOST_TEST_EQ(ellint_1(0.9_DDF, 1.2_DDF), 1.495205431757420_DDF);
    BOOST_TEST_EQ(ellint_2(0.3_DFF, 1.5_DFF), 1.467292_DFF);

    // k = 1: F(phi) = ln((1 + sin phi) / cos phi) below pi/2 and infinite beyond, and E(phi) = 2m + sin psi
    BOOST_TEST_EQ(ellint_1(1_DD, 0.5_DD), 0.5222381032784403_DD);
    BOOST_TEST_EQ(ellint_1(1_DL, 1.5_DL), "3.340677542798311003320812669037689"_DL);
    BOOST_TEST_EQ(ellint_1(1_DF, 1.570796_DF), 15.62708_DF);
    BOOST_TEST_EQ(ellint_1(1_DL, 1E-20_DL), 1E-20_DL);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isinf(ellint_1(1_DD, 2_DD)));
    #endif
    BOOST_TEST_EQ(ellint_2(1_DD, 0.5_DD), 0.4794255386042030_DD);
    BOOST_TEST_EQ(ellint_2(1_DL, 7_DL), "4.656986598718789090396999091593635"_DL);

    // k next to 1, and values where the polynomials of sin and cos are too close to a rounding boundary
    BOOST_TEST_EQ(ellint_1(0.9999999_DF, 1.5_DF), 3.340668_DF);
    BOOST_TEST_EQ(ellint_1(0.2006_DL, 7.519_DL), "7.593163036177046711748955109304458"_DL);
    BOOST_TEST_EQ(ellint_2(0.3777_DL, 5.751_DL), "5.523939034619259364629137849338906"_DL);

    // The complete integrals, also where the truncation leaves b above a in the AGM
    BOOST_TEST_EQ(comp_ellint_1(0.5_DL), "1.685750354812596042871203657799077"_DL);
    BOOST_TEST_EQ(comp_ellint_2(0.9_DD), 1.171697052781614_DD);
    BOOST_TEST_EQ(comp_ellint_1(0.9999999999999999_DD), 19.46040151479228_DD);
    BOOST_TEST_EQ(comp_ellint_2("0.99"_DL), "1.028475809028804000983887138518022"_DL);
    BOOST_TEST_EQ(comp_ellint_1(0_DF), 1.570796_DF);
    BOOST_TEST_EQ(comp_ellint_1(0.2158582_DF), 1.589590_DF);
    BOOST_TEST_EQ(comp_ellint_2(1_DD), 1_DD);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(comp_ellint_1(-1_DL), std::numeric_limits<decimal128_t>::infinity());
    #endif

    // Next to a rounding boundary, where the first try in one word cannot round
    BOOST_TEST_EQ(comp_ellint_1(0.5788157013434914_DD), 1.734964296895059_DD);
    BOOST_TEST_EQ(comp_ellint_2(0.8877832539736914_DD), 1.186404531467351_DD);

    // Zero, subnormal and too large arguments
    BOOST_TEST(signbit(ellint_1(0.5_DD, -0_DD)));
    BOOST_TEST_EQ(ellint_2(0_DF, 3_DF), 3_DF);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 1E-398_DD), 1E-398_DD);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isinf(ellint_1(0.9_DF, std::numeric_limits<decimal32_t>::max())));
    BOOST_TEST(isnan(ellint_1(1.5_DF, 1_DF)));
    BOOST_TEST_EQ(ellint_1(0.5_DD, -std::numeric_limits<decimal64_t>::infinity()), -std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(ellint_2(1_DF, std::numeric_limits<decimal32_t>::infinity()), std::numeric_limits<decimal32_t>::infinity());
    #endif
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(comp_ellint_1(0_DD), 1.570796326794897_DD);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 1E-20_DD), 1.000000000000001E-20_DD);
    BOOST_TEST_EQ(ellint_2(0.75_DD, 100_DD), 83.86482140325315_DD);

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(comp_ellint_1(0_DD), 1.570796326794896_DD);
    BOOST_TEST_EQ(ellint_2(0.5_DD, 1E-20_DD), 9.999999999999999E-21_DD);
    BOOST_TEST_EQ(ellint_2(0.999_DL, 1.57_DL), "1.003958804174344573363203047810723"_DL);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(ellint_1(0.3_DF, 2.5_DF), 2.570801_DF);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    test_nearest();
    test_directed();

    return boost::report_errors();
}
