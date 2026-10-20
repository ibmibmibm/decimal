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

    // Next to a halfway point
    BOOST_TEST_EQ(static_cast<long double>("4272846054360183744836159461173983e4170"_DL), 4272846054360183744836159461173983e4170L);
    BOOST_TEST_EQ(static_cast<long double>("9404211713773460935347478161609630e1010"_DL), 9404211713773460935347478161609630e1010L);
    BOOST_TEST_EQ(static_cast<long double>("7707109115044340911414827104497503e-2393"_DL), 7707109115044340911414827104497503e-2393L);
    BOOST_TEST_EQ(static_cast<long double>("7732408295769758696927081992574715e-1077"_DL), 7732408295769758696927081992574715e-1077L);
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
    BOOST_TEST_EQ(static_cast<long double>("-1.2e4932"_DL), -std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("1e5000"_DL), std::numeric_limits<long double>::infinity());
    BOOST_TEST_EQ(static_cast<long double>("3.6e-4951"_DL), std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("2e-4951"_DL), std::numeric_limits<long double>::denorm_min());
    BOOST_TEST_EQ(static_cast<long double>("1e-4940"_DL), 1e-4940L);
    BOOST_TEST_EQ(static_cast<long double>("1e-4900"_DL), 1e-4900L);
    BOOST_TEST_EQ(static_cast<long double>("1e-4952"_DL), 0.0L);
    BOOST_TEST_EQ(static_cast<long double>("1e-6000"_DL), 0.0L);
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
