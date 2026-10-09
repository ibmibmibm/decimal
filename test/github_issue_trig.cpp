// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The expected values are sin(x), cos(x) and tan(x) rounded once, from mpmath, MPFR and Python's decimal module
void test_nearest()
{
    // |x| <= pi/4: no reduction
    BOOST_TEST_EQ(sin(6.824108109971738738161997141037691e-01_DL), 6.306657757894080087699016123068339e-1_DL);
    BOOST_TEST_EQ(tan(9.186040435569681927699147918717055e-02_DL), 9.211966239425544881603018873528669e-2_DL);

    // |x| < 10^19: the reduction in binary
    BOOST_TEST_EQ(cos(7.885042103990243198911807988293644e-01_DL), 7.049070669070514707825184701610486e-1_DL);
    BOOST_TEST_EQ(cos(1.071126015241440004653838261529146e+04_DL), -5.083670249067169423296614073587474e-30_DL);
    BOOST_TEST_EQ(tan(1.277144282333587805192310314414510e+00_DL), 3.306939511260400241875095153806795e+0_DL);
    BOOST_TEST_EQ(tan(-9.537305793590876872804276295433355e+10_DL), 3.936100735101734918922338645464199e+0_DL);

    // The reduction in decimal
    BOOST_TEST_EQ(sin(-5.701712728788927995357933841156996e+5879_DL), -4.565990347796265476732945757361854e-1_DL);
    BOOST_TEST_EQ(cos(-9.246464983940743092998069521960153e+2067_DL), 1.564728004882864171853874563553649e-1_DL);
    BOOST_TEST_EQ(sin(9.212725e+65_DF), 9.563225e-9_DF);
    BOOST_TEST_EQ(sin(8.649991364154930e+52_DD), 3.698204338807975e-9_DD);
    BOOST_TEST_EQ(sin(6.974660901835021e+265_DD), -9.824396390797517e-9_DD);

    // decimal64_t next to a rounding boundary, where the first try in one word cannot round, and with 18 and 20 digits
    BOOST_TEST_EQ(sin(0.2367942983969950_DD), 0.2345875906801154_DD);
    BOOST_TEST_EQ(sin(0.9238113093010455_DD), 0.7979048039652705_DD);
    BOOST_TEST_EQ(cos(1.042468791386650_DD), 0.5040896205795680_DD);
    BOOST_TEST_EQ(tan(0.1296559321398583_DD), 0.1303873849787592_DD);
    BOOST_TEST_EQ(tan(4.242340898305579_DD), 1.968401691368704_DD);
    BOOST_TEST_EQ(sin(0.1001_DD), 0.09993291656402304_DD);
    BOOST_TEST_EQ(tan(0.0999_DD), 0.1002336663939089_DD);
    BOOST_TEST_EQ(cos(5.320793723314151e+15_DD), -0.8117436903221579_DD);

    // decimal128_t next to a rounding boundary, where the first try in two words cannot round, and with r too small
    BOOST_TEST_EQ(sin(4593_DL), -0.008459447378284051417990871368510194_DL);
    BOOST_TEST_EQ(cos(74.93_DL), 0.8923713493729135921361696906118304_DL);
    BOOST_TEST_EQ(tan(0.7709_DL), 0.9714160847029317332669884763063492_DL);
    BOOST_TEST_EQ(tan(11.01_DL), -69.31585780828732057320138727212295_DL);
    BOOST_TEST_EQ(cos(1.570796326794896619231321691639751_DL), 4.420985846996875529104874722961539e-34_DL);
}

void test_upward()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(tan(5.792056618670528004052681944417116e-02_DL), 5.798542370527202861717590633225899e-2_DL);
    BOOST_TEST_EQ(sin(4.618948766205371668526340711172817e+02_DL), -8.066879421873760570477739878656002e-2_DL);
    BOOST_TEST_EQ(sin(8.014804287421515201607723180384336e+13_DL), 9.063164340778351549755439430432617e-21_DL);
    BOOST_TEST_EQ(cos(9.311540218536312138490858104372267e+05_DL), 9.634687136646756869115694764091029e-2_DL);
    BOOST_TEST_EQ(tan(7.880488110646563439300105039991022e+03_DL), 5.008716172654112002643963997679707e+0_DL);
    BOOST_TEST_EQ(cos(4.475948463855351000179310105809870e+532_DL), 8.967217740227326360445186954223630e-1_DL);
    BOOST_TEST_EQ(sin(2.000000000061343_DD), 0.9092974268001540_DD);
    BOOST_TEST_EQ(cos(9.662753229975339_DD), -0.9718172673345117_DD);
    BOOST_TEST_EQ(tan(4.572877287488234_DD), 7.121293619524882_DD);
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
        BOOST_TEST(is_quiet_of(sin(s), s));
        BOOST_TEST(is_quiet_of(cos(s), s));
        BOOST_TEST(is_quiet_of(tan(s), s));
    }
    #endif
}

int main()
{
    test_nearest();
    test_upward();
    test_snan<decimal32_t>();
    test_snan<decimal64_t>();
    test_snan<decimal128_t>();
    test_snan<decimal_fast32_t>();
    test_snan<decimal_fast64_t>();
    test_snan<decimal_fast128_t>();

    return boost::report_errors();
}
