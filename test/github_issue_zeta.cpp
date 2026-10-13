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

// The expected values are zeta(x) rounded once, from MPFR
void test_nearest()
{
    // x >= 0 by Euler-Maclaurin, also next to the pole at 1 where 1 / (x - 1) is an exact division
    BOOST_TEST_EQ(riemann_zeta(2.5_DD), 1.341487257250917_DD);
    BOOST_TEST_EQ(riemann_zeta(0.25_DL), "-0.8132784052618916565214478200735256"_DL);
    BOOST_TEST_EQ(riemann_zeta(1.001_DF), 1000.577_DF);
    BOOST_TEST_EQ(riemann_zeta(0.999_DL), "-999.4228571557887900099207604196947"_DL);
    BOOST_TEST_EQ(riemann_zeta(1.25_DD), 4.595111825842943_DD);
    BOOST_TEST_EQ(riemann_zeta<decimal32_t>(3), 1.202057_DF);
    BOOST_TEST_EQ(riemann_zeta(2.5_DDF), 1.341487257250917_DDF);

    // Next to a rounding boundary, where the first try in one word cannot round
    BOOST_TEST_EQ(riemann_zeta(0.8847344349874664_DD), -8.106859235038839_DD);
    BOOST_TEST_EQ(riemann_zeta(1.836612874317278_DD), 1.829883296033672_DD);
    BOOST_TEST_EQ(riemann_zeta(6.475079494978107_DD), 1.012226563231261_DD);
    BOOST_TEST_EQ(riemann_zeta(5.44106e-8_DF), -0.5000000_DF);

    // x < 0 by the functional equation, with sin(pi x / 2) from the distance to the next even integer
    BOOST_TEST_EQ(riemann_zeta(-0.25_DL), "-0.3204512642285772827904444930548725"_DL);
    BOOST_TEST_EQ(riemann_zeta(-0.75_DD), -0.1336427744365846_DD);
    BOOST_TEST_EQ(riemann_zeta(-5.5_DF), -0.002671458_DF);
    BOOST_TEST_EQ(riemann_zeta("-1.9999999999"_DL), "-3.044845706168144658976147765095535E-12"_DL);
    BOOST_TEST_EQ(riemann_zeta(-3_DD), 0.008333333333333333_DD);
    BOOST_TEST_EQ(riemann_zeta(-199.5_DL), "7.270321237219865278251936806011302E+213"_DL);
    BOOST_TEST_EQ(riemann_zeta(-70.5_DF), 1.213177E+44_DF);
    BOOST_TEST_EQ(riemann_zeta(-1.25_DL), "-0.04890886762685481250801734162719492"_DL);
    BOOST_TEST_EQ(riemann_zeta(-2.75_DD), 0.009147150069079764_DD);
    BOOST_TEST_EQ(riemann_zeta("-1.234567890123456789012345678901234E-6"_DL), "-0.4999988655095229084763728920966713"_DL);
    BOOST_TEST_EQ(riemann_zeta(-5.5_DFF), -0.002671458_DFF);

    // Exact zeros, values beyond the range, and subnormal arguments
    BOOST_TEST_EQ(riemann_zeta(-4_DF), 0_DF);
    BOOST_TEST_EQ(riemann_zeta(-30_DD), 0_DD);
    BOOST_TEST_EQ(riemann_zeta(-300_DL), 0_DL);
    BOOST_TEST_EQ(riemann_zeta(-1E+20_DD), 0_DD);
    BOOST_TEST_EQ(riemann_zeta(-12345670_DF), 0_DF);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(riemann_zeta(-1234567_DF), std::numeric_limits<decimal32_t>::infinity());
    BOOST_TEST_EQ(riemann_zeta(-1001_DD), -std::numeric_limits<decimal64_t>::infinity());
    BOOST_TEST_EQ(riemann_zeta(-12345.5_DL), -std::numeric_limits<decimal128_t>::infinity());
    #endif
    BOOST_TEST_EQ(riemann_zeta(std::numeric_limits<decimal64_t>::denorm_min()), -0.5_DD);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isnan(riemann_zeta(-std::numeric_limits<decimal32_t>::infinity())));
    BOOST_TEST_EQ(riemann_zeta(std::numeric_limits<decimal32_t>::infinity()), 1_DF);
    #endif
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(riemann_zeta(1E-40_DD), -0.5_DD);
    BOOST_TEST_EQ(riemann_zeta(300_DF), 1.000001_DF);
    BOOST_TEST_EQ(riemann_zeta(100_DD), 1.000000000000001_DD);
    BOOST_TEST_EQ(riemann_zeta(0.5_DL), "-1.460354508809586812889499152515298"_DL);

    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(riemann_zeta(-1E-40_DD), -0.5_DD);
    BOOST_TEST_EQ(riemann_zeta(0.75_DF), -3.441286_DF);
    BOOST_TEST_EQ(riemann_zeta(20.5_DD), 1.000000674515618_DD);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

// x is the quiet form of the signaling NaN s, with the same sign and payload
template <typename T>
bool is_quiet_of(const T x, const T s)
{
    return isnan(x) && !issignaling(x) && signbit(x) == signbit(s) && read_payload(x) == read_payload(s);
}

// A signaling NaN operand gives a quiet NaN (IEEE 754 6.2)
template <typename T>
void test_snan()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto s : {snan<T>("7"), -snan<T>("7")})
    {
        BOOST_TEST(is_quiet_of(riemann_zeta(s), s));
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
