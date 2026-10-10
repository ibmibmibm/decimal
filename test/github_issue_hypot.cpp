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

// The expected values are the square roots of the exact sums of squares rounded once, from Python's decimal module
void test_nearest()
{
    BOOST_TEST_EQ(hypot(0.8872412_DF, 88.89712_DF), 88.90155_DF);
    BOOST_TEST_EQ(hypot(-3.5_DF, 0_DF), 3.5_DF);
    BOOST_TEST_EQ(hypot(-0_DD, -2.25_DD), 2.25_DD);
    BOOST_TEST_EQ(hypot(8.991880_DF, -489.1005_DF), 489.1831_DF);
    BOOST_TEST_EQ(hypot(403.3311651415152_DD, 35.19186680910160_DD), 404.8635526494349_DD);
    BOOST_TEST_EQ(hypot(-7292.570170779030_DD, 0.0006151718372943846_DD), 7292.570170779056_DD);
    BOOST_TEST_EQ(hypot("0.6411700099853032910128847337058788"_DL, "-0.7022765464166058658631056914405829"_DL), "0.9509423375533287797104854737288161"_DL);
    BOOST_TEST_EQ(hypot("0.002383928401138997542020258919149313"_DL, "-3896.558403728101209620822754491447"_DL), "3896.558403728830457592099046011541"_DL);
    BOOST_TEST_EQ(hypot(5854.919_DF, 47.23336_DF, -5641.964_DF), 8131.056_DF);
    BOOST_TEST_EQ(hypot(0.6641664868746812_DD, -7.312955409019041_DD, -0.002761570845410468_DD), 7.343053967039997_DD);
    BOOST_TEST_EQ(hypot("-0.08353962925839937517211652709796300"_DL, "2.268382463817469581476097682109652"_DL, "0.1444551254136481821536462822834949"_DL), "2.274512069668903150005847640381864"_DL);
    // The parts of the squares below the last digit add up to 1: halfway cases and a term below 10^-32
    BOOST_TEST_EQ(hypot(54018000_DF, 4190.04_DF, 40037.28_DF), 54018020_DF);
    BOOST_TEST_EQ(hypot(4190.04_DF, -40037.28_DF, 54018000_DF), 54018020_DF);
    BOOST_TEST_EQ(hypot(60000000600000000_DD, 139200014.04_DD, 1334400005.28_DD), 60000000600000020_DD);
    BOOST_TEST_EQ(hypot(4995328455517534_DD, 6.91013748054264E-5_DD, 70677637.591515_DD), 4995328455517535_DD);
    BOOST_TEST_EQ(hypot(8.991880_DFF, -489.1005_DFF), 489.1831_DFF);
    BOOST_TEST_EQ(hypot(403.3311651415152_DDF, 35.19186680910160_DDF), 404.8635526494349_DDF);
    BOOST_TEST_EQ(hypot("0.6411700099853032910128847337058788"_DLF, "-0.7022765464166058658631056914405829"_DLF), "0.9509423375533287797104854737288161"_DLF);

    // Exact values, zeros and negative values, subnormal values and special values
    BOOST_TEST_EQ(hypot(3_DD, 4_DD), 5_DD);
    BOOST_TEST_EQ(hypot(0_DD, -3_DD), 3_DD);
    BOOST_TEST_EQ(hypot(-3_DD, -4_DD, 0_DD), 5_DD);
    BOOST_TEST_EQ(hypot(-3_DL, -4_DL, -12_DL), 13_DL);
    BOOST_TEST_EQ(hypot(3E-398_DD, 4E-398_DD), 5E-398_DD);
    BOOST_TEST_EQ(hypot(1E-398_DD, 1E-398_DD), 1E-398_DD);
    BOOST_TEST_EQ(hypot(std::numeric_limits<decimal64_t>::infinity(), std::numeric_limits<decimal64_t>::quiet_NaN()), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST(isnan(hypot(1_DF, std::numeric_limits<decimal32_t>::quiet_NaN())));
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(hypot(0.3791849285048550_DD, -0.03844413802353041_DD), 0.3811287994282303_DD);
    BOOST_TEST_EQ(hypot(9E+384_DD, 9E+384_DD), 9.999999999999999E+384_DD);
    BOOST_TEST_EQ(hypot(36012000_DF, 2793.36_DF, 26691.52_DF), 36012010_DF);

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(hypot(0.04435630_DF, -43.03082_DF), 43.03085_DF);
    BOOST_TEST_EQ(hypot(-0.02494800724874129_DD, 0.0008362485772738227_DD), 0.02496201869573534_DD);
    BOOST_TEST_EQ(hypot("0.7343205781761127372048851959491316"_DL, "-0.006661145562372408908860666158430723"_DL), "0.7343507897409136793837601300649964"_DL);
    BOOST_TEST_EQ(hypot(1_DD, 1E-20_DD), 1.000000000000001_DD);
    BOOST_TEST_EQ(hypot(1E-398_DD, 1E-398_DD), 2E-398_DD);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(hypot("-181.2072722974621233109684629617206"_DL, "-0.00002677377758681178509196982087484668"_DL), "181.2072722974641012535675086975107"_DL);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    test_nearest();
    test_directed();

    return boost::report_errors();
}
