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

// The expected values are x^y rounded once, from Python's decimal module
void test_nearest()
{
    BOOST_TEST_EQ(pow(0.6433012_DF, -7.624039_DF), 28.88394_DF);
    BOOST_TEST_EQ(pow(0.4602037_DF, -24.41955_DF), 1.700618E+8_DF);
    BOOST_TEST_EQ(pow(29.76225_DF, -6.175466_DF), 7.932813E-10_DF);
    BOOST_TEST_EQ(pow(4.151952_DF, 2.634613_DF), 42.54569_DF);
    BOOST_TEST_EQ(pow(2.855110702918065_DD, 57.89312364130644_DD), 2.384711784845800E+26_DD);
    BOOST_TEST_EQ(pow(6.274273378840336_DD, 3.700058513418584_DD), 893.3700468676490_DD);
    BOOST_TEST_EQ(pow(0.7018756552090100_DD, -60.26699863059102_DD), 1842627634.736790_DD);
    BOOST_TEST_EQ(pow(6.353605224024699_DD, -81.77764720608027_DD), 2.142373526830494E-66_DD);
    BOOST_TEST_EQ(pow(1.811327746040189_DD, 3.46_DD), 7.810336716624140_DD);
    BOOST_TEST_EQ(pow(2.264244524000810_DD, 2.47_DD), 7.527670750896634_DD);
    BOOST_TEST_EQ(pow("70.01743037610289288015126977641270"_DL, "84.41844624783264454010499045632421"_DL), "5.880817842118858824560769070853572E+155"_DL);
    BOOST_TEST_EQ(pow("5.794566780408328697044313144907325"_DL, "-6.126784445420328523321714355288038"_DL), "0.00002114145997294436404045676779599843"_DL);
    BOOST_TEST_EQ(pow("6.551471311573427632363923451955033"_DL, "-6.353198101430488921311896592725542"_DL), "0.000006510848780098593528529454882607523"_DL);
    BOOST_TEST_EQ(pow("5.074669480018405513578494983547775"_DL, "-60.00414451978910529326681927416324"_DL), "4.705553807270951820478930074476548E-43"_DL);
    BOOST_TEST_EQ(pow(4.151952_DFF, 2.634613_DFF), 42.54569_DFF);
    BOOST_TEST_EQ(pow(6.274273378840336_DDF, 3.700058513418584_DDF), 893.3700468676490_DDF);
    BOOST_TEST_EQ(pow("5.794566780408328697044313144907325"_DLF, "-6.126784445420328523321714355288038"_DLF), "0.00002114145997294436404045676779599843"_DLF);

    // x near 1 with a large y, and integer powers
    BOOST_TEST_EQ(pow(0.9999357_DF, 1810.111_DF), 0.8901246_DF);
    BOOST_TEST_EQ(pow(1.000258_DF, -95133.58_DF), 2.197099E-11_DF);
    BOOST_TEST_EQ(pow(0.9998646853092840_DD, -616135.9085800826_DD), 1.623879133980048E+36_DD);
    BOOST_TEST_EQ(pow(0.9992721601885568_DD, 731847.9853306383_DD), 3.810267093839689E-232_DD);
    BOOST_TEST_EQ(pow(1.07_DF, 120), 3357.788_DF);
    BOOST_TEST_EQ(pow(3.7_DD, -25), 6.236729358194125E-15_DD);
    BOOST_TEST_EQ(pow(1.01_DL, 1000), "20959.15563781366006444124578894489"_DL);
    BOOST_TEST_EQ(pow(2.855110702918065_DD, 3), 23.27388350672476_DD);
    BOOST_TEST_EQ(pow(1.234567890123457E-130_DD, 3), 1.88167637E-390_DD);
    BOOST_TEST_EQ(pow(1.1_DF, 3), 1.331_DF);
    BOOST_TEST_EQ(pow(1.000000000000001_DD, 3), 1.000000000000003_DD);
    BOOST_TEST_EQ(pow("1.000000000000000000000000000000001"_DL, 3), "1.000000000000000000000000000000003"_DL);
    BOOST_TEST_EQ(pow(1.01_DL, -1000), "0.00004771184570984531863873560306841899"_DL);
    BOOST_TEST_EQ(pow(0.8_DF, -3), 1.953125_DF);
    BOOST_TEST_EQ(pow(2_DD, -10), 0.0009765625_DD);
    BOOST_TEST_EQ(pow(5.1_DD, 1.97_DD), 24.76927226232869_DD);
    BOOST_TEST_EQ(pow(2E+60_DD, 0.37_DD), 2.048241203535212E+22_DD);
    BOOST_TEST_EQ(pow(2_DD, 1E-20_DD), 1_DD);
    BOOST_TEST_EQ(pow(2_DD, 1.5E-120_DD), 1_DD);
    BOOST_TEST_EQ(pow(5.1_DL, 1.97_DL), "24.76927226232869473551915290362884"_DL);
    BOOST_TEST_EQ(pow(7.3_DL, -0.61_DL), "0.2974219764658053678518153058545516"_DL);
    BOOST_TEST_EQ(pow(1.5_DL, 25.5_DL), "30926.23886477426760604664875174766"_DL);
    BOOST_TEST_EQ(pow(2E+400_DL, 0.37_DL), "1.292352830637492244505565031970707E+148"_DL);

    // Exact values, a negative x, a subnormal x and a huge y
    BOOST_TEST_EQ(pow(2.25_DD, 1.5_DD), 3.375_DD);
    BOOST_TEST_EQ(pow(32_DF, 0.2_DF), 2_DF);
    BOOST_TEST_EQ(pow(1E+20_DL, 0.25_DL), 1E+5_DL);
    BOOST_TEST_EQ(pow(10_DL, 1000_DL), 1E+1000_DL);
    BOOST_TEST_EQ(pow(0.01_DL, 2000), 1E-4000_DL);
    BOOST_TEST_EQ(pow(-2_DD, 3_DD), -8_DD);
    BOOST_TEST_EQ(pow(-1_DL, 1E+20_DL), 1_DL);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isnan(pow(-8_DD, 0.5_DD)));
    #endif
    BOOST_TEST_EQ(pow(1.5E-384_DD, 0.5_DD), 1.224744871391589E-192_DD);
    BOOST_TEST_EQ(pow(1.5E-384_DD, 1), 1.5E-384_DD);
    BOOST_TEST_EQ(pow(0.1_DD, 1000000), 0_DD);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isinf(pow(2_DD, 1.5E+10_DD)));
    #endif
    BOOST_TEST_EQ(pow(2_DD, -1.5E+10_DD), 0_DD);
    BOOST_TEST_EQ(pow(1.000000000000001_DD, 0.5_DD), 1_DD);
    BOOST_TEST_EQ(pow("1.000000000000000000000000000000001"_DL, 0.5_DL), 1_DL);
    BOOST_TEST_EQ(pow("0.99999999999999998"_DL, 0.5_DL), "0.9999999999999999899999999999999999"_DL);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(pow(0.000003077052_DF, -2.037872_DF), 1.707952E+11_DF);
    BOOST_TEST_EQ(pow(6.175466E+15_DF, 0.4032085_DF), 2327728_DF);
    BOOST_TEST_EQ(pow("5470830892553768.575636583018230105"_DL, "-0.4950782376746678656120660229123392"_DL), "1.615966347871851091896748091332516E-8"_DL);
    BOOST_TEST_EQ(pow(2_DL, 0.5_DL), "1.414213562373095048801688724209698"_DL);
    BOOST_TEST_EQ(pow(4_DD, 0.5_DD), 2_DD);
    BOOST_TEST_EQ(pow("3.141592653589793238462643383279503"_DL, -7), "0.0003310936801775667643259528012530679"_DL);
    BOOST_TEST_EQ(pow(8.634450434451396_DD, -7.936486506691E-33_DD), 0.9999999999999999_DD);

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(pow(1.353146907159657E-16_DD, -6.026699863059102_DD), 4.321367699552815E+95_DD);
    BOOST_TEST_EQ(pow(5.473649078247373E+17_DD, 0.5109131778694365_DD), 1155365637.866424_DD);
    BOOST_TEST_EQ(pow("513200.5977674547163479848540467407"_DL, "-0.1836767961627070853012611652679634"_DL), "0.08936197332411062940049168630200383"_DL);
    BOOST_TEST_EQ(pow(7_DF, -3), 0.002915452_DF);
    BOOST_TEST_EQ(pow(0.1_DD, 1000000), 1E-398_DD);
    BOOST_TEST_EQ(pow(-6.274273378840336_DD, 3), -246.9962233269080_DD);
    BOOST_TEST_EQ(pow(1.000001_DF, 7), 1.000008_DF);
    BOOST_TEST_EQ(pow(1.000000000000001_DD, 3), 1.000000000000004_DD);
    BOOST_TEST_EQ(pow(2.855110702918065_DD, -3), 0.04296661533564263_DD);
    BOOST_TEST_EQ(pow("1.000000000000000000000000000000001"_DL, 3), "1.000000000000000000000000000000004"_DL);
    BOOST_TEST_EQ(pow(8.634450434451396_DD, 7.936486506691E-33_DD), 1.000000000000001_DD);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(pow(7.135241E+17_DF, -0.1973060_DF), 0.0003002032_DF);
    BOOST_TEST_EQ(pow(5611077.396074256_DD, -7.819717284820180_DD), 1.676356223809315E-53_DD);
    BOOST_TEST_EQ(pow(9.999999999999999E+127_DD, 3), 9.999999999999997E+383_DD);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

template <typename T>
bool is_neg(const T x, const T want)
{
    return x == want && signbit(x);
}

template <typename T>
bool is_pos(const T x, const T want)
{
    return x == want && !signbit(x);
}

// pow(-0, y) and pow(-inf, y) are negative only for an odd integer y, and pow(+1, y) is 1 also for a NaN y
template <typename T>
void test_special()
{
    const auto inf {std::numeric_limits<T>::infinity()};
    const auto nan {std::numeric_limits<T>::quiet_NaN()};
    const T zero {0};

    BOOST_TEST(is_neg(pow(-zero, T {-1}), -inf));
    BOOST_TEST(is_neg(pow(-zero, T {-3}), -inf));
    BOOST_TEST(is_neg(pow(-zero, -3), -inf));
    BOOST_TEST(is_neg(pow(-zero, T {30U, -1, true}), -inf));
    BOOST_TEST(is_pos(pow(-zero, T {-2}), inf));
    BOOST_TEST(is_pos(pow(-zero, T {25U, -1, true}), inf));
    BOOST_TEST(is_pos(pow(zero, T {-3}), inf));
    BOOST_TEST(is_neg(pow(-zero, T {3}), zero));
    BOOST_TEST(is_pos(pow(-zero, T {25U, -1}), zero));
    BOOST_TEST(is_neg(pow(-inf, T {-3}), zero));
    BOOST_TEST(is_neg(pow(-inf, T {3}), -inf));
    BOOST_TEST(is_pos(pow(-inf, T {25U, -1}), inf));
    BOOST_TEST(is_pos(pow(-inf, T {25U, -1, true}), zero));

    BOOST_TEST_EQ(pow(T {1}, nan), T {1});
    BOOST_TEST_EQ(pow(T {1}, -nan), T {1});
    BOOST_TEST(isnan(pow(T {-1}, nan)));
    BOOST_TEST(isnan(pow(nan, T {1})));
    BOOST_TEST_EQ(pow(T {-1}, inf), T {1});
}

// An integer y outside the range of int keeps its parity
template <typename T>
void test_large()
{
    const auto inf {std::numeric_limits<T>::infinity()};
    const T zero {0};
    const T odd {3000000001U};
    const T even {3000000002U};

    BOOST_TEST(is_neg(pow(-zero, -odd), -inf));
    BOOST_TEST(is_neg(pow(-zero, odd), zero));
    BOOST_TEST(is_pos(pow(-zero, even), zero));
    BOOST_TEST(is_neg(pow(-inf, odd), -inf));
    BOOST_TEST(is_neg(pow(-inf, -odd), zero));
    BOOST_TEST(is_pos(pow(-inf, even), inf));
    BOOST_TEST(is_neg(pow(T {-1}, odd), T {-1}));
    BOOST_TEST(is_pos(pow(T {-1}, even), T {1}));
    BOOST_TEST(is_neg(pow(T {-2}, odd), -inf));
    BOOST_TEST(is_neg(pow(T {5U, -1, true}, -odd), -inf));
    BOOST_TEST(is_pos(pow(T {-2}, T {1U, 20}), inf));
    BOOST_TEST(is_pos(pow(T {5U, -1, true}, T {1U, 20}), zero));
    BOOST_TEST(isnan(pow(T {-2}, T {30000000005U, -1})));
}

int main()
{
    test_nearest();
    test_directed();

    #ifndef BOOST_DECIMAL_FAST_MATH
    test_special<decimal32_t>();
    test_special<decimal64_t>();
    test_special<decimal128_t>();
    test_special<decimal_fast32_t>();
    test_special<decimal_fast64_t>();
    test_special<decimal_fast128_t>();

    test_large<decimal64_t>();
    test_large<decimal128_t>();
    test_large<decimal_fast64_t>();
    test_large<decimal_fast128_t>();
    #endif

    return boost::report_errors();
}
