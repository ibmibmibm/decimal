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

// The expected values are the polynomials rounded once, from exact rational arithmetic
void test_nearest()
{
    BOOST_TEST_EQ(hermite(10, 1.234567890123456_DD), -50861.15164935731_DD);
    BOOST_TEST_EQ(hermite(25, "-3.141592653589793238462643383279503"_DL), "-3.919595891027119179045962031752274E+17"_DL);
    BOOST_TEST_EQ(hermite(60, 9.876543_DF), -6.892005E+70_DF);
    BOOST_TEST_EQ(legendre(30, 0.7071067811865476_DD), -0.06638905244972245_DD);
    BOOST_TEST_EQ(legendre(127, "-0.123456789012345678901234567890123"_DL), "-0.005184979324728096484512022520990776"_DL);
    BOOST_TEST_EQ(laguerre(12, 7.5_DD), 5.767036298652748_DD);
    BOOST_TEST_EQ(laguerre(40, 33.25_DL), "1021021.399327567361779150104972568"_DL);
    BOOST_TEST_EQ(laguerre(3, -2.5_DD), 20.47916666666667_DD);
    BOOST_TEST_EQ(assoc_laguerre(9, 4, 2.71828182845904_DD), -0.3334471777185667_DD);
    BOOST_TEST_EQ(assoc_laguerre(3, 2, 0.001_DF), 9.990002_DF);
    BOOST_TEST_EQ(hermite(10, 1.234567890123456_DDF), -50861.15164935731_DDF);
    BOOST_TEST_EQ(assoc_legendre(1, 1, 0.5_DF), 0.8660254_DF);
    BOOST_TEST_EQ(assoc_legendre(60, 55, 0.3_DD), 2.103727831409394E+92_DD);
    BOOST_TEST_EQ(assoc_legendre(119, 81, "-0.7390575738807128036711446839035489"_DL), "1.164374619515172437698472871467166E+164"_DL);

    // Exact values, exact zeros, and beyond the range
    BOOST_TEST_EQ(legendre(5, 0.5_DF), 0.08984375_DF);
    BOOST_TEST_EQ(hermite(3, 0.5_DD), -5_DD);
    BOOST_TEST_EQ(legendre(3, 0_DD), 0_DD);
    BOOST_TEST_EQ(laguerre(1, 1_DD), 0_DD);
    BOOST_TEST_EQ(assoc_laguerre(1, 3, 4_DL), 0_DL);
    BOOST_TEST_EQ(assoc_legendre(3, 1, 0.6_DD), 0.96_DD);
    BOOST_TEST_EQ(assoc_legendre(2, 3, 0.5_DD), 0_DD);
    BOOST_TEST_EQ(assoc_legendre(4, 2, -1_DD), 0_DD);
    BOOST_TEST_EQ(laguerre(50, 1E+20_DD), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(hermite(3, -std::numeric_limits<decimal32_t>::infinity()), -std::numeric_limits<decimal32_t>::infinity());
    BOOST_TEST_EQ(laguerre(3, std::numeric_limits<decimal32_t>::infinity()), -std::numeric_limits<decimal32_t>::infinity());
    BOOST_TEST_EQ(assoc_laguerre(2, 1, -std::numeric_limits<decimal64_t>::infinity()), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST(isnan(hermite(2, std::numeric_limits<decimal32_t>::quiet_NaN())));
    BOOST_TEST(isnan(laguerre(2, std::numeric_limits<decimal64_t>::quiet_NaN())));
    BOOST_TEST(isnan(assoc_laguerre(2, 1, std::numeric_limits<decimal128_t>::quiet_NaN())));
    BOOST_TEST(isnan(assoc_legendre(2, 1, std::numeric_limits<decimal64_t>::quiet_NaN())));

    // Error intervals over several boundaries: two more words, then a search with the exact integers, which also
    // give the sign and the magnitude next to a root, and the exact zero of L_2^(2)(2)
    BOOST_TEST_EQ(laguerre(100, decimal128_t {1, -100}), 1_DL);
    BOOST_TEST_EQ(laguerre(300, 2.5_DD), -0.3097147065489419_DD);
    BOOST_TEST_EQ(laguerre(1000, 2.5_DD), 0.07368422332292386_DD);
    BOOST_TEST_EQ(assoc_laguerre(2, 2, 2_DD), 0_DD);
    BOOST_TEST(isfinite(laguerre(1000, 2.123456789012345_DD)));
    BOOST_TEST_EQ(hermite(140, 1E+60_DD), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(laguerre(199, 1E+60_DD), -std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(hermite(3, 1E-6000_DL), "-1.2E-5999"_DL);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(legendre(4, 0.5_DD), -0.2890625_DD);
    BOOST_TEST_EQ(legendre(128, 1_DL), 1_DL);
    BOOST_TEST_EQ(assoc_legendre(125, 9, "-0.9928170953853051594606426988320891"_DL), "-1.699935412697470245182432303480512E+18"_DL);
    BOOST_TEST_EQ(assoc_legendre(3, 1, 1E-50_DD), -1.499999999999999_DD);

    // A tiny or huge x, and where the exact integers are too large, the next term decides
    BOOST_TEST_EQ(legendre(2, 1E-300_DD), -0.4999999999999999_DD);
    BOOST_TEST_EQ(laguerre(1, -1E-300_DD), 1.000000000000001_DD);
    BOOST_TEST_EQ(hermite(1, -2.927649194442799E-112_DD), -5.855298388885598E-112_DD);
    BOOST_TEST_EQ(legendre(2, 1E-6000_DL), "-0.4999999999999999999999999999999999"_DL);
    BOOST_TEST_EQ(laguerre(2, -1E-6000_DL), "1.000000000000000000000000000000001"_DL);
    BOOST_TEST_EQ(laguerre(100, decimal128_t {1U, -100, true}), "1.000000000000000000000000000000001"_DL);
    BOOST_TEST_EQ(hermite(3, 1E-6000_DL), "-1.199999999999999999999999999999999E-5999"_DL);
    BOOST_TEST_EQ(legendre(3, 1E-6000_DL), "-1.499999999999999999999999999999999E-6000"_DL);

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(hermite(3, 0.5_DD), -5_DD);
    BOOST_TEST_EQ(legendre(100, 1_DD), 1_DD);
    BOOST_TEST_EQ(legendre(99, -1_DD), -1_DD);
    BOOST_TEST_EQ(laguerre(100, 0_DL), 1_DL);
    BOOST_TEST_EQ(laguerre(1, 2_DD), -1_DD);
    BOOST_TEST_EQ(assoc_laguerre(1, 2, 3_DF), 0_DF);
    BOOST_TEST_EQ(legendre(1, -0.25_DD), -0.25_DD);
    BOOST_TEST_EQ(hermite(3, 1E+200_DL), "7.999999999999999999999999999999999E+600"_DL);
    BOOST_TEST_EQ(assoc_laguerre(5, 3, 1E-300_DD), 55.99999999999999_DD);
    BOOST_TEST_EQ(assoc_legendre(7, 2, -0.8_DD), -7.1959104_DD);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(hermite(30, 1E-30_DF), -2.028432E+20_DF);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    test_nearest();
    test_directed();

    return boost::report_errors();
}
