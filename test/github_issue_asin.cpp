// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The expected values are asin(x) and acos(x) rounded once, from Python's decimal module
void test_nearest()
{
    BOOST_TEST_EQ(asin(-0.7998671_DF), -0.9270738_DF);
    BOOST_TEST_EQ(asin(1E-60_DD), 1.000000000000000E-60_DD);
    BOOST_TEST_EQ(asin(-0.003_DF), -0.003000005_DF);
    // Where two words leave the result to three, and one word to two
    BOOST_TEST_EQ(asin(0.104_DL), "0.1041883957412958806130176441808407"_DL);
    BOOST_TEST_EQ(acos(0.122_DL), "1.448491640287984430522321981114637"_DL);
    BOOST_TEST_EQ(asin(0.81_DD), 0.9441521151541560_DD);
    BOOST_TEST_EQ(acos(0.713_DD), 0.7770287725707690_DD);
    BOOST_TEST_EQ(acos(0.237_DD), 1.331519622217062_DD);
    BOOST_TEST_EQ(asin(0.116_DD), 0.1162617373233272_DD);
    BOOST_TEST_EQ(asin(0.5_DD), 0.5235987755982989_DD);
    BOOST_TEST_EQ(asin(0.99_DD), 1.429256853470469_DD);
    BOOST_TEST_EQ(asin(0.999_DD), 1.526071239626163_DD);
    BOOST_TEST_EQ(asin(-0.9999_DL), "-1.556654073317383741635081465822095"_DL);
    BOOST_TEST_EQ(asin(0.995_DL), "1.470754613183356589905832917998823"_DL);
    BOOST_TEST_EQ(asin(0.3_DF), 0.3046927_DF);
    BOOST_TEST_EQ(acos(0.003_DD), 1.567796322294878_DD);
    BOOST_TEST_EQ(acos(-0.003_DF), 1.573796_DF);
    BOOST_TEST_EQ(acos(0.99999_DD), 0.004472139681787927_DD);
    BOOST_TEST_EQ(acos(0.99999_DF), 0.004472140_DF);
    BOOST_TEST_EQ(acos("0.9999999999999999999999999999999999"_DL), "1.414213562373095048801688724209698E-17"_DL);
    BOOST_TEST_EQ(acos("0.06537041789163665433003886727831877"_DL), "1.505379261336168285048172435760927"_DL);
    BOOST_TEST_EQ(acos(0.5_DD), 1.047197551196598_DD);
    BOOST_TEST_EQ(acos(-0.5_DF), 2.094395_DF);
    BOOST_TEST_EQ(acos(-0.9999_DD), 3.127450400112280_DD);
    BOOST_TEST_EQ(asin(-0.7998671_DFF), -0.9270738_DFF);
    BOOST_TEST_EQ(acos(0.99999_DDF), 0.004472139681787927_DDF);
    BOOST_TEST_EQ(acos(0.9_DLF), "0.4510268117962624325446446357943518"_DLF);
    BOOST_TEST_EQ(acos(1_DD), 0_DD);
    BOOST_TEST_EQ(acos(0_DF), 1.570796_DF);
    BOOST_TEST_EQ(acos(-0_DD), 1.570796326794897_DD);
    BOOST_TEST_EQ(acos(0_DL), "1.570796326794896619231321691639751"_DL);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(asin(-0.1107962_DF), -0.1110242_DF);
    BOOST_TEST_EQ(acos(0.1929245540962392_DD), 1.376654500811053_DD);
    BOOST_TEST_EQ(asin("-0.9993596794963899153161094199321795"_DL), "-1.535008372107914734058462530021444"_DL);
    BOOST_TEST_EQ(acos("0.6270410123727878600243485968237583"_DL), "0.8930474601568815172728723528362396"_DL);
    BOOST_TEST_EQ(acos(0.9_DL), "0.4510268117962624325446446357943518"_DL);
    BOOST_TEST_EQ(asin(1_DD), 1.570796326794896_DD);
    BOOST_TEST_EQ(acos(-1_DF), 3.141592_DF);

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(acos(0.7244133_DF), 0.7606134_DF);
    BOOST_TEST_EQ(asin(0.75_DL), "0.8480620789814810080529443389984181"_DL);
    BOOST_TEST_EQ(asin(-1_DF), -1.570796_DF);
    BOOST_TEST_EQ(acos(-1_DF), 3.141593_DF);
    BOOST_TEST_EQ(acos(0_DF), 1.570797_DF);

    fesetround(rounding_mode::fe_dec_to_nearest_from_zero);
    BOOST_TEST_EQ(asin(-0.5800143320018901_DD), -0.6187462843400824_DD);

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
        BOOST_TEST(is_quiet_of(asin(s), s));
        BOOST_TEST(is_quiet_of(acos(s), s));
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
