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

// The expected values are atan(x) and atan2(y, x) rounded once, from Python's decimal module
void test_nearest()
{
    BOOST_TEST_EQ(atan("-0.5892067299081487533440595915372664"_DL), "-0.5324454702390584674534036153852551"_DL);
    // Where two words leave the result to three, and one word to two
    BOOST_TEST_EQ(atan(1.07_DL), "0.8192017070909034209463922188695074"_DL);
    BOOST_TEST_EQ(atan2(1_DL, 11_DL), "0.09065988720074511349838667530846562"_DL);
    BOOST_TEST_EQ(atan(2_DD), 1.107148717794091_DD);
    BOOST_TEST_EQ(atan2(2_DD, 1_DD), 1.107148717794091_DD);
    BOOST_TEST_EQ(atan2(5.084777E-8_DF, -1.221514E-12_DF), 1.570820_DF);
    BOOST_TEST_EQ(atan2(-9.889701E-11_DF, -5.737004E-12_DF), -1.628741_DF);
    BOOST_TEST_EQ(atan2(-1234.564356728940_DD, -1189779.630989923_DD), -3.140555012768344_DD);
    BOOST_TEST_EQ(atan2(-6.862956853956393E-10_DD, -0.000005364629725944915_DD), -3.141464723850611_DD);
    BOOST_TEST_EQ(atan2("7.008931184813850285755811658047442E-8"_DL, "-0.04225572226405051731999031722715888"_DL), "3.141590994895820145416852322544857"_DL);
    BOOST_TEST_EQ(atan2("-835981859.6245254303217756935217128"_DL, "494960112168673.4338021395666511274"_DL), "-0.000001688988342840125745029701458576222"_DL);
    BOOST_TEST_EQ(atan(1E-60_DD), 1.000000000000000E-60_DD);
    BOOST_TEST_EQ(atan(-0.003_DF), -0.002999991_DF);
    BOOST_TEST_EQ(atan(0.5_DD), 0.4636476090008061_DD);
    BOOST_TEST_EQ(atan(-0.99_DL), "-0.7803730800666358988978715172725503"_DL);
    BOOST_TEST_EQ(atan(1_DF), 0.7853982_DF);
    BOOST_TEST_EQ(atan(1E+300_DD), 1.570796326794897_DD);
    BOOST_TEST_EQ(atan2(1_DD, 1E+40_DD), 1.000000000000000E-40_DD);
    BOOST_TEST_EQ(atan2(1_DL, 300_DL), "0.003333320987736624861198407405281840"_DL);
    BOOST_TEST_EQ(atan2(7_DD, 3_DD), 1.165904540509813_DD);
    BOOST_TEST_EQ(atan2(3_DD, 7_DD), 0.4048917862850834_DD);
    BOOST_TEST_EQ(atan2(-1_DF, -3_DF), -2.819842_DF);
    BOOST_TEST_EQ(atan2(2_DD, 2_DD), 0.7853981633974483_DD);
    BOOST_TEST_EQ(atan2(-1234.564356728940_DDF, -1189779.630989923_DDF), -3.140555012768344_DDF);
    BOOST_TEST_EQ(atan2(5.084777E-8_DFF, -1.221514E-12_DFF), 1.570820_DFF);
    BOOST_TEST_EQ(atan("-0.5892067299081487533440595915372664"_DLF), "-0.5324454702390584674534036153852551"_DLF);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(atan(-2.206814_DF), -1.145333_DF);
    BOOST_TEST_EQ(atan(-394952671107993.2_DD), -1.570796326794895_DD);
    BOOST_TEST_EQ(atan("-36664.68652901188993426419044905825"_DL), "-1.570769052594721451683076171619218"_DL);
    BOOST_TEST_EQ(atan2(9.258818E+11_DF, -3.596536E+11_DF), 1.941301_DF);
    BOOST_TEST_EQ(atan2(929.7802769310495_DD, 85804.76753509128_DD), 0.01083557376775061_DD);
    BOOST_TEST_EQ(atan2("8194.175721799908004419907635340890"_DL, "-515.3449464535053570459986576380949"_DL), "1.633605218675387025744968097611059"_DL);
    BOOST_TEST_EQ(atan(1E-60_DD), 9.999999999999999E-61_DD);
    BOOST_TEST_EQ(atan2(1_DD, 1E+40_DD), 9.999999999999999E-41_DD);
    BOOST_TEST_EQ(atan2(1_DF, 3E+20_DF), 3.333333E-21_DF);
    BOOST_TEST_EQ(atan2(5_DL, -2E-30_DL), "1.570796326794896619231321691640151"_DL);
    BOOST_TEST_EQ(atan2(0_DF, -1_DF), 3.141592_DF);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(atan(std::numeric_limits<decimal64_t>::infinity()), 1.570796326794896_DD);
    BOOST_TEST_EQ(atan2(-std::numeric_limits<decimal64_t>::infinity(), 1_DD), -1.570796326794897_DD);
    #endif

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(atan(7.5_DL), "1.438244794498222597961404247935482"_DL);
    BOOST_TEST_EQ(atan2(1_DF, 3E+20_DF), 3.333334E-21_DF);
    BOOST_TEST_EQ(atan2(1E-30_DD, -1_DD), 3.141592653589794_DD);
    BOOST_TEST_EQ(atan2(1_DF, 0_DF), 1.570797_DF);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(atan2(std::numeric_limits<decimal32_t>::infinity(), -std::numeric_limits<decimal32_t>::infinity()), 2.356195_DF);
    BOOST_TEST_EQ(atan2(std::numeric_limits<decimal64_t>::infinity(), std::numeric_limits<decimal64_t>::infinity()), 0.7853981633974484_DD);
    BOOST_TEST_EQ(atan2(1_DD, -std::numeric_limits<decimal64_t>::infinity()), 3.141592653589794_DD);
    #endif

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(atan2(-2_DL, 3E+40_DL), "-6.666666666666666666666666666666666E-41"_DL);

    fesetround(rounding_mode::fe_dec_to_nearest_from_zero);


    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    test_nearest();
    test_directed();

    return boost::report_errors();
}
