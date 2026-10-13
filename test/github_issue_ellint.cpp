// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include "quiet_nan.hpp"
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The expected values are the integrals rounded once, from MPFR
void test_nearest()
{
    // Arguments made to put the result next to a rounding boundary, closer than the last stage told apart
    BOOST_TEST_EQ(ellint_2("9991084742670618412903739033515207e-34"_DL, "7990428525308400572612436592297663e-41"_DL), "7990428525308392085025258743391852e-41"_DL);
    BOOST_TEST_EQ(ellint_1("5781511315714688498080653313135854e-39"_DL, "3146485975252437351042033248266420e-35"_DL), "3146485975252454701978381849466136e-35"_DL);
    #if BOOST_DECIMAL_DEC_EVAL_METHOD == 0 // A wider evaluation type rounds twice
    BOOST_TEST_EQ(ellint_1(3093652014346328e-17_DD, 6487570877300547e-20_DD), 6487570877304903e-20_DD);
    #endif
    BOOST_TEST_EQ(ellint_2(8239611643240809e-16_DD, 1276191344458994e-20_DD), 1276191344435476e-20_DD);
    BOOST_TEST_EQ(ellint_1(2230806e-7_DF, 3921025e-9_DF), 3921026e-9_DF);
    BOOST_TEST_EQ(ellint_2(2227397e-7_DF, 3925025e-9_DF), 3925025e-9_DF);
    // phi below pi/4, between pi/4 and pi/2, and above, where F(phi) = 2K/pi phi + (F(psi) - 2K/pi psi)
    // Next to a rounding boundary: the sin and cos polynomials of the type, and phi^5 past k^2 phi^3 / 6
    BOOST_TEST_EQ(ellint_2(0.3_DD, 1e-7_DD), 9.999999999999999e-8_DD);
    BOOST_TEST_EQ(ellint_1(0.8_DD, 7.5e-8_DD), 7.500000000000005e-8_DD);
    BOOST_TEST_EQ(ellint_1(0.9757347_DF, 9.814307e-4_DF), 9.814309e-4_DF);
    BOOST_TEST_EQ(ellint_1(1_DL, 1.5e-15_DL), "1.500000000000000000000000000000563e-15"_DL);
    BOOST_TEST_EQ(ellint_1(0.75_DL, 6e-16_DL), "6.000000000000000000000000000000203e-16"_DL);
    // The terms to phi^5 add up to a tie, and the term in phi^7 takes the value below it
    BOOST_TEST_EQ(ellint_2(0.3_DL, 1e-7_DL), "9.999999999999998500000000000002797e-8"_DL);
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
    // decimal32_t next to a rounding boundary, where the AGM that stops early does not decide
    BOOST_TEST_EQ(ellint_1(0.9639035_DF, 0.7262805_DF), 0.7933881_DF);
    BOOST_TEST_EQ(ellint_1(0.5459291_DF, 0.7367823_DF), 0.7558178_DF);
    BOOST_TEST_EQ(ellint_2(0.965988_DF, 0.668615_DF), 0.6233886_DF);
    BOOST_TEST_EQ(ellint_2(0.5010182_DF, 9.443887_DF), 8.821212_DF);
    // F next to phi, where the first try in one word for decimal64_t does not round
    BOOST_TEST_EQ(ellint_1(1E-10_DD, 0.7_DD), 0.7_DD);
    // E next to a rounding boundary, where the first try in two words for decimal128_t cannot round
    BOOST_TEST_EQ(ellint_2(0.11_DL, 0.58_DL), "0.5796319597982352716674100393140373"_DL);

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
    BOOST_TEST_EQ(comp_ellint_1(0.2952638675311015_DD), 1.606819968759654_DD);
    BOOST_TEST_EQ(comp_ellint_2(0.751633862523325_DD), 1.317178696265717_DD);

    // Zero, subnormal and too large arguments
    BOOST_TEST(signbit(ellint_1(0.5_DD, -0_DD)));
    BOOST_TEST_EQ(ellint_2(0_DF, 3_DF), 3_DF);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 1E-398_DD), 1E-398_DD);
    // A tiny phi needs no reduction and no 1 / phi, which has no power of 10 in one word, also in a constant expression
    #ifdef __GNUC__
    static_assert(ellint_1(0.5_DD, 1E-30_DD) == 1E-30_DD, "");
    static_assert(ellint_2(0.5_DD, 1E-30_DD) == 1E-30_DD, "");
    #endif
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isinf(ellint_1(0.9_DF, std::numeric_limits<decimal32_t>::max())));
    BOOST_TEST(isnan(ellint_1(1.5_DF, 1_DF)));
    BOOST_TEST_EQ(ellint_1(0.5_DD, -std::numeric_limits<decimal64_t>::infinity()), -std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(ellint_2(1_DF, std::numeric_limits<decimal32_t>::infinity()), std::numeric_limits<decimal32_t>::infinity());
    #endif
    // A result next to a rounding boundary that is a constant plus a small term, where the values of the argument next to
    // each other give results that three words do not tell apart
    BOOST_TEST_EQ(ellint_1("-3006656887792573886287950455042040e-46"_DL, "1e0"_DL), "1000000000000000000000000012324920e-33"_DL);
    BOOST_TEST_EQ(ellint_2("-3435607827448873799365644526120561e-48"_DL, "1e0"_DL), "9999999999999999999999999999983908e-34"_DL);
    BOOST_TEST_EQ(ellint_1("1424559085601768668508201297357079e-45"_DL, "3e-2"_DL), "3000000000000000000000000000913052e-35"_DL);
    BOOST_TEST_EQ(ellint_2("5151118648942231785607464430572648e-45"_DL, "12345e0"_DL), "1234499999999999999999991810863587e-29"_DL);
    BOOST_TEST_EQ(ellint_1("-1727692722859575369195825634630249e-34"_DL, "3e-7"_DL), "3000000000000001343214965079847523e-40"_DL);
    BOOST_TEST_EQ(ellint_2("1173227696183664486913706064303173e-34"_DL, "3e-7"_DL), "2999999999999999380591547808418011e-40"_DL);
    BOOST_TEST_EQ(ellint_1("7131589166237953094544880530782265e-39"_DL, "5e-7"_DL), "5000000000000000000000010595742508e-40"_DL);
    BOOST_TEST_EQ(comp_ellint_1("1214268270659259244023318575718854e-50"_DL), "1570796326794896619231321691639752e-33"_DL);
    BOOST_TEST_EQ(comp_ellint_2("1121826666826860737516374540562821e-49"_DL), "1570796326794896619231321691639747e-33"_DL);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_upward);
    // Arguments made to put the result next to a rounding boundary, closer than the last stage told apart
    BOOST_TEST_EQ(ellint_1("3633787057067939556879135838586963e-34"_DL, "1590613736138992618535875461667164e-42"_DL), "1590613736138992618624440411872762e-42"_DL);
    BOOST_TEST_EQ(ellint_2("6317249044032936991496778567695123e-34"_DL, "1440887931162022915196221336389572e-42"_DL), "1440887931162022914997247782247813e-42"_DL);
    BOOST_TEST_EQ(ellint_1(2008574869854960e-17_DD, 4987088600075363e-20_DD), 4987088600076198e-20_DD);
    BOOST_TEST_EQ(ellint_1("9344498040759280280859708036279182e-47"_DL, "1e0"_DL), "1000000000000000000000000001190498e-33"_DL);
    BOOST_TEST_EQ(ellint_2("-1351990642936400695598115948509336e-47"_DL, "1e0"_DL), "9999999999999999999999999999750792e-34"_DL);
    BOOST_TEST_EQ(ellint_1(0.1138621322655823_DD, 0.8071231987_DD), 0.8081233363326548_DD);
    BOOST_TEST_EQ(ellint_1(0.9_DL, 3e-12_DL), "3.000000000000000000000003645000001e-12"_DL);
    // Next to pi/2 and 3 pi/2, E(1, phi) is an odd integer plus or minus 1 - cos(r), below three words
    BOOST_TEST_EQ(ellint_2(1_DL, "1.570796326794896619231321691639754"_DL), "1.000000000000000000000000000000001"_DL);
    BOOST_TEST_EQ(ellint_2(1_DL, "4.712388980384689857693965074919255"_DL), "3.000000000000000000000000000000001"_DL);
    BOOST_TEST_EQ(comp_ellint_1(0_DD), 1.570796326794897_DD);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 1E-20_DD), 1.000000000000001E-20_DD);
    BOOST_TEST_EQ(ellint_2(0.75_DD, 100_DD), 83.86482140325315_DD);
    // The terms to phi^5 add up to a value of the type, and the term in phi^7 takes F above it and E below it
    BOOST_TEST_EQ(ellint_1(1_DL, 6e-8_DL), "6.000000000000003600000000000003241e-8"_DL);
    BOOST_TEST_EQ(ellint_1(0.5_DD, 3E-398_DD), 4E-398_DD);

    fesetround(rounding_mode::fe_dec_downward);
    // Arguments made to put the result next to a rounding boundary, closer than the last stage told apart
    BOOST_TEST_EQ(ellint_1("3633787057067939556879135838586963e-34"_DL, "1590613736138992618535875461667164e-42"_DL), "1590613736138992618624440411872761e-42"_DL);
    BOOST_TEST_EQ(ellint_2("6317249044032936991496778567695123e-34"_DL, "1440887931162022915196221336389572e-42"_DL), "1440887931162022914997247782247812e-42"_DL);
    BOOST_TEST_EQ(ellint_1(2008574869854960e-17_DD, 4987088600075363e-20_DD), 4987088600076197e-20_DD);
    BOOST_TEST_EQ(ellint_2(1_DL, 6e-8_DL), "5.999999999999996400000000000000647e-8"_DL);
    BOOST_TEST_EQ(ellint_1(0.3_DL, 2e-7_DL), "2.000000000000001199999999999992344e-7"_DL);
    BOOST_TEST_EQ(ellint_1(0.9_DL, 3e-12_DL), "3.000000000000000000000003645e-12"_DL);
    BOOST_TEST_EQ(ellint_2(1_DL, "1.570796326794896619231321691639754"_DL), 1_DL);
    BOOST_TEST_EQ(ellint_2(1_DL, "4.712388980384689857693965074919250"_DL), "2.999999999999999999999999999999999"_DL);
    BOOST_TEST_EQ(comp_ellint_1(0_DD), 1.570796326794896_DD);
    BOOST_TEST_EQ(ellint_2(0.5_DD, 1E-20_DD), 9.999999999999999E-21_DD);
    BOOST_TEST_EQ(ellint_2(0.5_DD, 3E-398_DD), 2E-398_DD);
    BOOST_TEST_EQ(ellint_2(0.999_DL, 1.57_DL), "1.003958804174344573363203047810723"_DL);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(ellint_1(0.3_DF, 2.5_DF), 2.570801_DF);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

// A signaling NaN operand gives a quiet NaN with its payload (IEEE 754 6.2)
template <typename T>
void test_snan()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto s : {snan<T>("7"), -snan<T>("7")})
    {
        BOOST_TEST(is_quiet_of(ellint_1(s, T {1}), s));
        BOOST_TEST(is_quiet_of(ellint_2(T {1, -1}, s), s));
        BOOST_TEST(is_quiet_of(comp_ellint_1(s), s));
        BOOST_TEST(is_quiet_of(comp_ellint_2(s), s));
    }
    #endif
}

int main()
{
    test_nearest();
    test_directed();
    test_snan<decimal32_t>();
    test_snan<decimal64_t>();
    test_snan<decimal128_t>();
    test_snan<decimal_fast32_t>();
    test_snan<decimal_fast64_t>();
    test_snan<decimal_fast128_t>();

    return boost::report_errors();
}
