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

// An infinity and a zero keep the sign of the decimal value
template <typename T>
void test_sign()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(static_cast<double>(-std::numeric_limits<T>::infinity()), -std::numeric_limits<double>::infinity());
    BOOST_TEST_EQ(static_cast<float>(-std::numeric_limits<T>::infinity()), -std::numeric_limits<float>::infinity());
    #endif

    BOOST_TEST(std::signbit(static_cast<double>(T(0U, 0, true))));
    BOOST_TEST(std::signbit(static_cast<float>(T(0U, 0, true))));
}

// A value above the largest binary value is infinity, not 0
void test_overflow()
{
    BOOST_TEST_EQ(static_cast<double>("1e309"_DD), std::numeric_limits<double>::infinity());
    BOOST_TEST_EQ(static_cast<double>("-2e400"_DL), -std::numeric_limits<double>::infinity());
    BOOST_TEST_EQ(static_cast<float>("3.4028236e38"_DD), std::numeric_limits<float>::infinity());
    BOOST_TEST_EQ(static_cast<float>("3.4028235e38"_DD), (std::numeric_limits<float>::max)());
    BOOST_TEST_EQ(static_cast<float>("1e39"_DF), std::numeric_limits<float>::infinity());
    BOOST_TEST_EQ(static_cast<double>("1e309"_DDF), std::numeric_limits<double>::infinity());
    BOOST_TEST_EQ(static_cast<double>("2e308"_DD), std::numeric_limits<double>::infinity());
}

// A value in the subnormal range of the binary type is subnormal, not 0
void test_subnormal()
{
    BOOST_TEST_EQ(static_cast<double>("5e-324"_DD), std::numeric_limits<double>::denorm_min());
    BOOST_TEST_EQ(static_cast<double>("2.470328229206233e-324"_DD), std::numeric_limits<double>::denorm_min());
    BOOST_TEST_EQ(static_cast<double>("2.470328229206232e-324"_DD), 0.0);
    BOOST_TEST_EQ(static_cast<double>("1e-320"_DD), 1e-320);
    BOOST_TEST_EQ(static_cast<double>("1e-310"_DL), 1e-310);
    BOOST_TEST_EQ(static_cast<float>("1e-45"_DF), std::numeric_limits<float>::denorm_min());
    BOOST_TEST_EQ(static_cast<float>("7e-46"_DF), 0.0F);
    BOOST_TEST_EQ(static_cast<double>("5e-324"_DLF), std::numeric_limits<double>::denorm_min());
    BOOST_TEST_EQ(static_cast<double>("1e-342"_DL), 0.0);
    BOOST_TEST_EQ(static_cast<double>("1e-400"_DL), 0.0);
    BOOST_TEST_EQ(static_cast<double>("2.2250738585072013e-308"_DL), (std::numeric_limits<double>::min)());
}

// A tie rounds to the even value
void test_ties()
{
    BOOST_TEST_EQ(static_cast<double>("4e23"_DF), 4e23);
    BOOST_TEST_EQ(static_cast<double>("9007199254740993"_DD), 9007199254740992.0);
    BOOST_TEST_EQ(static_cast<double>("9007199254740995"_DD), 9007199254740996.0);
    BOOST_TEST_EQ(static_cast<float>("16777217"_DD), 16777216.0F);
    BOOST_TEST_EQ(static_cast<float>("16777219"_DD), 16777220.0F);

    // (2^53 + 1) / 2^25 has 34 digits and is halfway between 2^28 and the next double
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953125e-25"_DL), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953126e-25"_DL), 268435456.00000006);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953124e-25"_DL), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953125e-25"_DLF), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("9007199254740991.9"_DL), 9007199254740992.0);
    BOOST_TEST_EQ(static_cast<float>("0.5"_DF), 0.5F);
    BOOST_TEST_EQ(static_cast<double>("0.1"_DD), 0.1);
}

// The 19 leading digits of these values do not decide the rounding, so an exact compare does
void test_wide_significand()
{
    BOOST_TEST_EQ(static_cast<double>("1234567890123456789012345678901234e-20"_DL), 12345678901234.568);
    BOOST_TEST_EQ(static_cast<double>("1112856818313216103296202118680802e265"_DL), 1.112856818313216e+298);
    BOOST_TEST_EQ(static_cast<double>("1193075813432793838432928155755821e-308"_DL), 1.1930758134327938e-275);
    BOOST_TEST_EQ(static_cast<double>("3880824273833896262359420750490712e-35"_DL), 0.03880824273833896);
    BOOST_TEST_EQ(static_cast<double>("6099573892221131429277060931893378e-351"_DL), 6.09957e-318);
    BOOST_TEST_EQ(static_cast<double>("6099573892221131429277060931893379e-351"_DL), 6.099576e-318);
    BOOST_TEST_EQ(static_cast<float>("1058826583522949505694339754652399e-63"_DL), 1.05882654e-30F);
    BOOST_TEST_EQ(static_cast<float>("1058826583522949505694339754652400e-63"_DL), 1.05882663e-30F);
    BOOST_TEST_EQ(static_cast<float>("3402820728259992105141534601966245e5"_DL), 3.40282063e+38F);
    BOOST_TEST_EQ(static_cast<float>("3402820728259992105141534601966246e5"_DLF), 3.40282083e+38F);
}

int main()
{
    test_sign<decimal32_t>();
    test_sign<decimal64_t>();
    test_sign<decimal128_t>();
    test_sign<decimal_fast32_t>();
    test_sign<decimal_fast64_t>();
    test_sign<decimal_fast128_t>();

    test_overflow();
    test_subnormal();
    test_ties();
    test_wide_significand();

    return boost::report_errors();
}
