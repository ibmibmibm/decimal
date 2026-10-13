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

// The expected values are the polynomials rounded once, from exact rational arithmetic
void test_nearest()
{
    // An odd degree and x < 0 where |H_n'(0) x| is next above 1: the sign is that of the lowest term
    BOOST_TEST_EQ(hermite(409, "-3532663936220064713110529972899553e-542"_DL), -1_DL);
    BOOST_TEST_EQ(hermite(409, "3532663936220064713110529972899553e-542"_DL), 1_DL);
    BOOST_TEST_EQ(hermite(10, 1.234567890123456_DD), -50861.15164935731_DD);
    BOOST_TEST_EQ(hermite(25, "-3.141592653589793238462643383279503"_DL), "-3.919595891027119179045962031752274E+17"_DL);
    BOOST_TEST_EQ(hermite(60, 9.876543_DF), -6.892005E+70_DF);
    BOOST_TEST_EQ(legendre(30, 0.7071067811865476_DD), -0.06638905244972245_DD);
    BOOST_TEST_EQ(legendre(127, "-0.123456789012345678901234567890123"_DL), "-0.005184979324728096484512022520990776"_DL);
    BOOST_TEST_EQ(laguerre(12, 7.5_DD), 5.767036298652748_DD);
    BOOST_TEST_EQ(laguerre(40, 33.25_DL), "1021021.399327567361779150104972568"_DL);
    BOOST_TEST_EQ(laguerre(3, -2.5_DD), 20.47916666666667_DD);
    BOOST_TEST_EQ(legendre(85, 0.5_DD), 0.06562026401747063_DD);
    BOOST_TEST_EQ(laguerre(100, 2.5_DD), 0.4257949400323675_DD);
    BOOST_TEST_EQ(laguerre(400, 1.234567890123456E-15_DD), 0.9999999999995062_DD);
    BOOST_TEST_EQ(assoc_laguerre(9, 4, 2.71828182845904_DD), -0.3334471777185667_DD);
    BOOST_TEST_EQ(assoc_laguerre(3, 2, 0.001_DF), 9.990002_DF);
    // The sum of the degree and the order does not fit in 32 bits
    BOOST_TEST_EQ(assoc_laguerre(2U, 4294967295U, 0.5_DF), 9.223372e18_DF);
    BOOST_TEST_EQ(assoc_laguerre(2U, 4294967295U, 0.5_DD), 9.223372036854776e18_DD);
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
    BOOST_TEST_EQ(legendre(3, -1_DL), -1_DL);
    BOOST_TEST_EQ(legendre(7, 0.5_DL), 0.22314453125_DL);
    BOOST_TEST_EQ(legendre(13, 0.999_DD), 0.9110275860123822_DD);
    BOOST_TEST_EQ(assoc_legendre(22, 20, 0.3_DD), 7.327742826737401E+24_DD);
    BOOST_TEST_EQ(hermite(25, 0_DD), 0_DD);
    BOOST_TEST_EQ(legendre(25, -0_DL), 0_DL);
    BOOST_TEST_EQ(assoc_legendre(31, 6, 0_DF), 0_DF);
    BOOST_TEST_EQ(hermite(24, 0_DD), 1295295050649600_DD);
    BOOST_TEST_EQ(legendre(18, 0.123_DD), 0.1215430990459902_DD);
    BOOST_TEST_EQ(laguerre(5, 2.345678901234567_DD), 0.9873245091799757_DD);
    BOOST_TEST_EQ(legendre(17, -0.37_DF), -0.06826659_DF);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(laguerre(50, 1E+20_DD), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(hermite(60, -1E+60_DD), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(hermite(3, -std::numeric_limits<decimal32_t>::infinity()), -std::numeric_limits<decimal32_t>::infinity());
    BOOST_TEST_EQ(laguerre(3, std::numeric_limits<decimal32_t>::infinity()), -std::numeric_limits<decimal32_t>::infinity());
    BOOST_TEST_EQ(assoc_laguerre(2, 1, -std::numeric_limits<decimal64_t>::infinity()), std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST(isnan(hermite(2, std::numeric_limits<decimal32_t>::quiet_NaN())));
    BOOST_TEST(isnan(laguerre(2, std::numeric_limits<decimal64_t>::quiet_NaN())));
    BOOST_TEST(isnan(assoc_laguerre(2, 1, std::numeric_limits<decimal128_t>::quiet_NaN())));
    BOOST_TEST(isnan(assoc_legendre(2, 1, std::numeric_limits<decimal64_t>::quiet_NaN())));
    BOOST_TEST(isnan(legendre(2, 1.5_DD)));
    BOOST_TEST(isnan(assoc_legendre(3, 1, -1.000001_DF)));
    // A degree past the highest one, and the highest one with an overflow in the recurrence
    BOOST_TEST(isnan(hermite(110000, decimal128_t {1, 6000})));
    BOOST_TEST(isnan(laguerre(65536, 2.5_DD)));
    BOOST_TEST(isnan(assoc_laguerre(65536, 3, 2.5_DF)));
    // The same for an infinite x and for L_n(0), which have a value at each lower degree
    BOOST_TEST(isnan(hermite(65536, std::numeric_limits<decimal64_t>::infinity())));
    BOOST_TEST(isnan(laguerre(65536, -std::numeric_limits<decimal64_t>::infinity())));
    BOOST_TEST(isnan(laguerre(65536, 0_DD)));
    BOOST_TEST_EQ(laguerre(65535, 0_DD), 1_DD);
    BOOST_TEST_EQ(hermite(65535, decimal128_t {-1, 6000}), -std::numeric_limits<decimal128_t>::infinity());
    // A value far past the range, where the decimal exponent comes from a binary exponent of millions
    BOOST_TEST_EQ(hermite(5000, decimal128_t {1, 6000}), std::numeric_limits<decimal128_t>::infinity());
    // Laguerre with a large x past the degree where the bounds of the try fit in an int
    BOOST_TEST_EQ(laguerre(500, decimal128_t {1, 6000}), std::numeric_limits<decimal128_t>::infinity());
    BOOST_TEST_EQ(assoc_laguerre(500, 3, decimal128_t {-1, 6000}), std::numeric_limits<decimal128_t>::infinity());
    BOOST_TEST_EQ(laguerre(65535, decimal64_t {1, 300}), -std::numeric_limits<decimal64_t>::infinity());
    #endif
    // H_1 and P_1 keep the sign of a zero
    BOOST_TEST(signbit(hermite(1, -0_DF)));
    BOOST_TEST(signbit(legendre(1, -0_DD)));
    BOOST_TEST(signbit(assoc_legendre(1, 0, -0_DL)));
    BOOST_TEST(!signbit(legendre(1, 0_DD)));

    // Error intervals over several boundaries: two more words, then a search with the exact integers, which also
    // give the sign and the magnitude next to a root, and the exact zero of L_2^(2)(2)
    BOOST_TEST_EQ(laguerre(100, decimal128_t {1, -100}), 1_DL);
    BOOST_TEST_EQ(laguerre(300, 2.5_DD), -0.3097147065489419_DD);
    BOOST_TEST_EQ(laguerre(1000, 2.5_DD), 0.07368422332292386_DD);
    BOOST_TEST_EQ(assoc_laguerre(2, 2, 2_DD), 0_DD);
    BOOST_TEST_EQ(laguerre(1000, 2.123456789012345_DD), -0.2306874989065326_DD);
    // The exact integers do not fit for these, two of them next to a root of the polynomial, with reference values from
    // exact rational arithmetic
    BOOST_TEST_EQ(laguerre(1000, 2.038850628305442_DD), 8.440403839487275E-15_DD);
    BOOST_TEST_EQ(laguerre(1000, 705.6526952668715409384449865450955_DL), 5.597998142965618403696468933222349E+120_DL);
    BOOST_TEST_EQ(hermite(400, 3.170218297556186580538647930045266_DL), -5.534274518357665441073560695687485E+464_DL);
    BOOST_TEST_EQ(laguerre(2000, 2.150102_DF), -0.000009289581_DF);
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
    BOOST_TEST_EQ(laguerre(3, 0.7_DF), -0.4221666_DF);
    BOOST_TEST_EQ(laguerre(65535, decimal128_t {1, 6000}), -(std::numeric_limits<decimal128_t>::max)());
    BOOST_TEST_EQ(assoc_laguerre(501, 3, decimal128_t {-1, 6000}), std::numeric_limits<decimal128_t>::infinity());

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(hermite(409, "-3532663936220064713110529972899553e-542"_DL), "-1000000000000000000000000000000001e-33"_DL);
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
    BOOST_TEST_EQ(laguerre(3, 0.7_DF), -0.4221667_DF);
    // The value is the binomial P(0) less a tiny term, and P(0) takes no integer of the recurrence at this degree
    BOOST_TEST_EQ(laguerre(65535, std::numeric_limits<decimal32_t>::denorm_min()), 0.9999999_DF);
    BOOST_TEST_EQ(laguerre(65535, std::numeric_limits<decimal64_t>::denorm_min()), 0.9999999999999999_DD);
    BOOST_TEST_EQ(laguerre(65535, std::numeric_limits<decimal128_t>::denorm_min()), "0.9999999999999999999999999999999999"_DL);
    BOOST_TEST_EQ(assoc_laguerre(4000, 3, std::numeric_limits<decimal64_t>::denorm_min()), 10682674000.99999_DD);
    // L_1^m(x) = m + 1 - x for an m + 1 that the type holds
    BOOST_TEST_EQ(assoc_laguerre(1, 10000000U, 1E-100_DD), 10000000.99999999_DD);
    BOOST_TEST_EQ(assoc_laguerre(1, 4294967295U, 0_DD), 4294967296_DD);
    BOOST_TEST_EQ(assoc_laguerre(500, 3, decimal128_t {1, 6000}), (std::numeric_limits<decimal128_t>::max)());
    BOOST_TEST_EQ(laguerre(3000, 1E+100_DD), (std::numeric_limits<decimal64_t>::max)());

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(hermite(30, 1E-30_DF), -2.028432E+20_DF);
    BOOST_TEST_EQ(hermite(41, 1E+60_DF), (std::numeric_limits<decimal32_t>::max)());
    // A large x, where the term of R_(k-1) is 2^64 or 2^128 times below the scale of the step
    BOOST_TEST_EQ(hermite(3, 1E+10_DF), 7.999999E+30_DF);
    BOOST_TEST_EQ(hermite(3, 1E+20_DF), 7.999999E+60_DF);
    BOOST_TEST_EQ(laguerre(500, decimal128_t {1, 6000}), (std::numeric_limits<decimal128_t>::max)());
    BOOST_TEST_EQ(assoc_laguerre(65535, 2, decimal128_t {1, 6000}), -(std::numeric_limits<decimal128_t>::max)());

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

// A signaling NaN operand gives a quiet NaN (IEEE 754 6.2)
template <typename T>
void test_snan()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto s : {snan<T>("7"), -snan<T>("7")})
    {
        BOOST_TEST(is_quiet_of(hermite(3U, s), s));
        BOOST_TEST(is_quiet_of(laguerre(3U, s), s));
        BOOST_TEST(is_quiet_of(legendre(3U, s), s));
        BOOST_TEST(is_quiet_of(legendre(200U, s), s));
        BOOST_TEST(is_quiet_of(assoc_legendre(200U, 1U, s), s));
        BOOST_TEST(is_quiet_of(assoc_laguerre(3U, 1U, s), s));
        BOOST_TEST(is_quiet_of(assoc_legendre(3U, 1U, s), s));
    }
    #endif
}

// 10^k from pf_pow10 has an error below 2 units, also where squaring takes |k| far past the table: two words against
// the top two words of three, at the ends of the exponent range of decimal128_t
void test_pow10_bound()
{
    for (const int k : {39, -39, 700, -700, 6144, -6144, 6200, -6200})
    {
        const auto a {detail::poly_detail::pf_pow10<2>(k)};
        const auto b {detail::poly_detail::pf_pow10<3>(k)};
        const auto top {detail::fx_top<2>(b.m)};
        const boost::int128::uint128_t ma {a.m.w[1], a.m.w[0]};
        const boost::int128::uint128_t mb {top.w[1], top.w[0]};
        BOOST_TEST_EQ(a.e, b.e);
        BOOST_TEST((ma > mb ? ma - mb : mb - ma) <= 3U);
    }
}

// The polynomial functions are noexcept, as the other cmath functions
static_assert(noexcept(hermite(2U, decimal64_t {})) && noexcept(legendre(2U, decimal64_t {})) && noexcept(laguerre(2U, decimal64_t {})), "");
static_assert(noexcept(assoc_laguerre(2U, 1U, decimal64_t {})) && noexcept(assoc_legendre(2U, 1U, decimal64_t {})), "");

// Next to zero, the bound of the sum of two terms has a large shift, which must stay below the width of int
void test_bound_shift()
{
    #ifdef __GNUC__
    static_assert(laguerre(3U, 1e-30_DD) == 1_DD, "");
    #endif
    BOOST_TEST_EQ(laguerre(3U, 1e-30_DD), 1_DD);
    BOOST_TEST_EQ(laguerre(5U, 6e-8_DD), 0.9999997000000180_DD);
    BOOST_TEST_EQ(laguerre(5U, 3e-8_DD), 0.9999998500000045_DD);
    BOOST_TEST_EQ(hermite(2U, 1e-30_DD), -2_DD);
}

// Laguerre above degree 22 with a small x, in the form of poly_loop2_laguerre in two and three words
void test_laguerre_diff()
{
    BOOST_TEST_EQ(laguerre(30, 1.5_DF), 0.4560499_DF);
    BOOST_TEST_EQ(laguerre(30, 2.5_DD), -0.4140453741003571_DD);
    BOOST_TEST_EQ(laguerre(40, -3.25_DD), 163519916.7065242_DD);
    BOOST_TEST_EQ(assoc_laguerre(20, 5, 1.5_DD), -260.4354784650168_DD);
    BOOST_TEST_EQ(laguerre(30, 2.5_DL), "-0.4140453741003570988087785940265470"_DL);
    BOOST_TEST_EQ(assoc_laguerre(25, 3, 0.75_DL), "-57.69330531853356818516493312330138"_DL);
    BOOST_TEST_EQ(laguerre(44, -0.375_DL), "416.6277738269558349590927422082194"_DL);
    // Past degree 128 the tries keep diff, so the exact integers seldom run
    BOOST_TEST_EQ(laguerre(300, 1.234567_DF), 0.2382837_DF);
    BOOST_TEST_EQ(laguerre(800, 1.234567_DF), 0.1365404_DF);
    BOOST_TEST_EQ(laguerre(800, 1.23456_DD), 0.1365174855088369_DD);
    BOOST_TEST_EQ(laguerre(1000, 2.5_DD), 0.07368422332292386_DD);
    BOOST_TEST_EQ(assoc_laguerre(500, 3, 2.5_DD), -792.0115958609465_DD);
    BOOST_TEST_EQ(laguerre(800, 1.23456_DL), "0.1365174855088369439097883938418519"_DL);
    BOOST_TEST_EQ(laguerre(2000, 1.23456_DL), "-0.05177876762140553963725069881383381"_DL);
    BOOST_TEST_EQ(assoc_laguerre(1500, 2, 0.75_DL), "263.1761984570871758848655471852810"_DL);
}

// decimal32_t past degree 22 in one word, for the sizes of x where one word decides, and legendre next to n = 16
void test_one_word()
{
    BOOST_TEST_EQ(hermite(40, 0.5_DF), -8.096548e28_DF);
    BOOST_TEST_EQ(hermite(60, 0.25_DF), -2.991731e49_DF);
    BOOST_TEST_EQ(hermite(30, 12.5_DF), 1.854804e41_DF);
    BOOST_TEST_EQ(hermite(24, 0.875_DF), 1.878274e15_DF);
    BOOST_TEST_EQ(legendre(60, 0.05_DF), -0.1019617_DF);
    BOOST_TEST_EQ(legendre(90, 0.03125_DF), 0.07981659_DF);
    BOOST_TEST_EQ(legendre(18, 0.75_DF), 0.2279831_DF);
    BOOST_TEST_EQ(legendre(16, 0.75_DF), 0.03263066_DF);
    BOOST_TEST_EQ(laguerre(20, -1.5_DF), 3866.299_DF);
    BOOST_TEST_EQ(laguerre(28, -0.125_DF), 8.754920_DF);
    BOOST_TEST_EQ(assoc_laguerre(15, 5, -2.5_DF), 1182021_DF);
    BOOST_TEST_EQ(assoc_legendre(40, 2, 0.05_DF), 89.99389_DF);
}

// An approximation just above a power of ten for a value below it, where the boundaries are 10 times closer
void test_resolve_below_power_of_ten()
{
    using boost::int128::uint128_t;
    using namespace boost::decimal::detail;
    const auto p37 {pow10(static_cast<uint128_t>(37))};
    const auto p29 {pow10(static_cast<uint128_t>(29))};
    // L_1(-8.9999987) = 9.9999987 from 10.00000001 with an error bound of 0.0000002
    const poly_detail::poly_out r {p37 + p29, p37, p29 * 20U};
    poly_detail::poly_decimal out {};
    BOOST_TEST(!poly_detail::poly_resolve_digits<false>(poly_detail::poly_kind::laguerre, 1U, 0U, {uint128_t {89999987U}, -7, true}, {r, -36, false}, 7, -101, out));
    BOOST_TEST_EQ(fx_round<decimal32_t>(out.d, out.e10, out.neg), 9.999999_DF);
}

// Next to zero, hermite takes the exact resolver, which GCC 8 does not evaluate at compile time
void test_resolve_constexpr()
{
    #if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 9)
    static_assert(hermite(2U, -0.7071067811865475244008443621048490_DL) == -2.222285911147240187798653227487960e-34_DL, "");
    #endif
    // A factor that the compiler cannot fold keeps the call at run time, also with GCC 8
    volatile int one {1};
    BOOST_TEST_EQ(hermite(2U, -0.7071067811865475244008443621048490_DL * one), -2.222285911147240187798653227487960e-34_DL);
}

int main()
{
    test_pow10_bound();
    test_nearest();
    test_bound_shift();
    test_laguerre_diff();
    test_one_word();
    test_directed();
    test_resolve_below_power_of_ten();
    test_resolve_constexpr();
    test_snan<decimal32_t>();
    test_snan<decimal64_t>();
    test_snan<decimal128_t>();
    test_snan<decimal_fast32_t>();
    test_snan<decimal_fast64_t>();
    test_snan<decimal_fast128_t>();

    return boost::report_errors();
}
