// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cmath>
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

#ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

// The long double literal is correctly rounded, so it is the expected value
void test_rounding()
{
    // One multiply or divide by an exact power of ten
    BOOST_TEST_EQ(static_cast<long double>("1.23456789e28"_DD), 1.23456789e28L);
    BOOST_TEST_EQ(static_cast<long double>("1.23456789e-10"_DD), 1.23456789e-10L);

    // decimal32_t does not go through double
    BOOST_TEST_EQ(static_cast<long double>("0.1"_DF), 0.1L);
    BOOST_TEST_EQ(static_cast<long double>("0.1"_DFF), 0.1L);
    BOOST_TEST_EQ(static_cast<long double>("4.366480e70"_DF), 4.366480e70L);

    // The sign of a negative value is kept
    BOOST_TEST_EQ(static_cast<long double>("-3.691692477e-289"_DD), -3.691692477e-289L);
    BOOST_TEST_EQ(static_cast<long double>("-3.691692477e-289"_DDF), -3.691692477e-289L);

    // All 34 digits count
    BOOST_TEST_EQ(static_cast<long double>("2.455594866531838259546e-290"_DL), 2.455594866531838259546e-290L);
    BOOST_TEST_EQ(static_cast<long double>("1.436648e-296"_DD), 1.436648e-296L);
}

// Values next to a halfway point, which need the exact compare
void test_halfway()
{
    BOOST_TEST_EQ(static_cast<long double>("18446744073709551617"_DL), 18446744073709551617.0L);
    BOOST_TEST_EQ(static_cast<long double>("1759218604441600000095367431640625e-20"_DL), 17592186044416.00000095367431640625L);
    BOOST_TEST_EQ(static_cast<long double>("27284841053187847137451171875e-40"_DL), 2.7284841053187847137451171875e-12L);
    BOOST_TEST_EQ(static_cast<long double>("184467440737095516159e-1"_DL), 18446744073709551615.9L);
}

#if BOOST_DECIMAL_LDBL_BITS == 80 || BOOST_DECIMAL_LDBL_BITS == 128

// Exponents beyond the range of double, in steps of the power table
void test_wide()
{
    BOOST_TEST_EQ(static_cast<long double>("8.825512804910594220703691692477000e-2154"_DL), 8.825512804910594220703691692477000e-2154L);
    BOOST_TEST_EQ(static_cast<long double>("3.410453893849661490551322660000000e1290"_DLF), 3.410453893849661490551322660000000e1290L);
    BOOST_TEST_EQ(static_cast<long double>("612583205011400993339800e-854"_DL), 612583205011400993339800e-854L);
    BOOST_TEST_EQ(static_cast<long double>("37052364e1606"_DL), 37052364e1606L);
    BOOST_TEST_EQ(static_cast<long double>("-1.28928e33"_DF), -1.28928e33L);
    BOOST_TEST_EQ(static_cast<long double>("-1.6796764983166403208704e22"_DLF), -1.6796764983166403208704e22L);
    BOOST_TEST_EQ(static_cast<long double>("-1.026932920667407516028613022449664e33"_DLF), -1.026932920667407516028613022449664e33L);
    BOOST_TEST_EQ(static_cast<long double>("-2.26744390569329649286929451390465e-4932"_DLF), -2.26744390569329649286929451390465e-4932L);

    // Next to a halfway point
    BOOST_TEST_EQ(static_cast<long double>("4272846054360183744836159461173983e4170"_DL), 4272846054360183744836159461173983e4170L);
    BOOST_TEST_EQ(static_cast<long double>("9404211713773460935347478161609630e1010"_DL), 9404211713773460935347478161609630e1010L);
    BOOST_TEST_EQ(static_cast<long double>("7707109115044340911414827104497503e-2393"_DL), 7707109115044340911414827104497503e-2393L);
    BOOST_TEST_EQ(static_cast<long double>("7732408295769758696927081992574715e-1077"_DL), 7732408295769758696927081992574715e-1077L);
    BOOST_TEST_EQ(static_cast<long double>("7.504351730091768941268813904223603e4618"_DL), 7.504351730091768941268813904223603e4618L);
    BOOST_TEST_EQ(static_cast<long double>("9.922419968056034721472704337650005e-4615"_DL), 9.922419968056034721472704337650005e-4615L);
    BOOST_TEST_EQ(static_cast<long double>("-6.489884414031349728144312857682782e895"_DL), -6.489884414031349728144312857682782e895L);

    // Each power of ten is 10 times the one before it, to a few units in the last place
    auto previous {static_cast<long double>(decimal128_t {1, -4901})};
    for (int exponent {-4900}; exponent <= 4900; ++exponent)
    {
        const auto power {static_cast<long double>(decimal128_t {1, exponent})};
        if (!BOOST_TEST_LT(std::fabs(power / previous - 10.0L), 1e-17L))
        {
            break;
        }
        previous = power;
    }
}

#else

void test_wide()
{
}

#endif

#if BOOST_DECIMAL_LDBL_BITS == 80

// The ends of the range of the 80-bit format
void test_range()
{
    BOOST_TEST_EQ(static_cast<long double>("1.18973149535723176502e4932"_DL), (std::numeric_limits<long double>::max)());
    // Next to a rounding boundary above 2^16384, so that an exact compare gives the infinity
    BOOST_TEST_EQ(static_cast<long double>("1.189731495357255948308539271593628e4932"_DL), std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("-1.2e4932"_DL), -std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("1e5000"_DL), std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("3.6e-4951"_DL), std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("2e-4951"_DL), std::numeric_limits<long double>::denorm_min());
    // Just below and just above 1.5 times the smallest subnormal value
    BOOST_TEST_EQ(static_cast<long double>("5467799297823711903792608900429129e-4984"_DL), std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("5467799297823711903792608900429130e-4984"_DL), 2 * std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("1e-4940"_DL), 1e-4940L);
    BOOST_TEST_EQ(static_cast<long double>("1e-4900"_DL), 1e-4900L);
    BOOST_TEST_EQ(static_cast<long double>("1e-4952"_DL), 0.0L);
    BOOST_TEST_EQ(static_cast<long double>("1e-6000"_DL), 0.0L);
}

#elif BOOST_DECIMAL_LDBL_BITS == 128

// The ends of the range of the binary128 format
void test_range()
{
    // Just below and just above the midpoint of the largest value and the next power of 2
    BOOST_TEST_EQ(static_cast<long double>("1.189731495357231765085759326628007e4932"_DL), (std::numeric_limits<long double>::max)());
    BOOST_TEST_EQ(static_cast<long double>("1.189731495357231765085759326628008e4932"_DL), std::numeric_limits<long double>::infinity());
    // Next to a rounding boundary above 2^16384, so that an exact compare gives the infinity
    BOOST_TEST_EQ(static_cast<long double>("1.189731495357231765085759327985373e4932"_DL), std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("-1.2e4932"_DL), -std::numeric_limits<long double>::infinity());
    // The largest subnormal value and the next value above the smallest normal value
    BOOST_TEST_EQ(static_cast<long double>("3.362103143112093506262677817321752e-4932"_DL), 3.362103143112093506262677817321752e-4932L);
    BOOST_TEST_EQ(static_cast<long double>("3.362103143112093506262677817321753e-4932"_DL), 3.362103143112093506262677817321753e-4932L);
    BOOST_TEST_EQ(static_cast<long double>("1e-4960"_DL), 1e-4960L);
    BOOST_TEST_EQ(static_cast<long double>("3.3e-4966"_DL), std::numeric_limits<long double>::denorm_min());
    // Just below and just above 1.5 times the smallest subnormal value
    BOOST_TEST_EQ(static_cast<long double>("9712762679157037666386658437341469e-4999"_DL), std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("9712762679157037666386658437341470e-4999"_DL), 2 * std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("3.2e-4966"_DL), 0.0L);
    // Values that need the exact compare, a power of 10 at the end of the exact range, and a value below the maximum
    BOOST_TEST_EQ(static_cast<long double>("3.523266e-469"_DL), 3.523266e-469L);
    BOOST_TEST_EQ(static_cast<long double>("-1.38781597797754629010400373425272e-1531"_DL), -1.38781597797754629010400373425272e-1531L);
    BOOST_TEST_EQ(static_cast<long double>("3.8468e53"_DL), 3.8468e53L);
    BOOST_TEST_EQ(static_cast<long double>("5.0867e4925"_DL), 5.0867e4925L);
}

#else

void test_range()
{
}

#endif

int main()
{
    test_rounding();
    test_halfway();
    test_wide();
    test_range();

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
