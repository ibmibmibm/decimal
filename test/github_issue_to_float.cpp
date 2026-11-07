// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cmath>
#include <cstring>
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
    BOOST_TEST_EQ(static_cast<double>(T(0U, 0, true)), 0.0);
    BOOST_TEST_EQ(static_cast<double>(T(0U, 0)), 0.0);
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

    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(std::isnan(static_cast<double>(std::numeric_limits<decimal64_t>::signaling_NaN())));
    #endif
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

    // A coefficient above 64 bits next to half the smallest subnormal value, and far below it
    BOOST_TEST_EQ(static_cast<double>("2470328229206232720882843964341106e-357"_DL), 0.0);
    BOOST_TEST_EQ(static_cast<double>("2470328229206232720882843964341107e-357"_DL),
                  std::numeric_limits<double>::denorm_min());
    BOOST_TEST_EQ(static_cast<double>("1234567890123456789012345678901234e-400"_DL), 0.0);
    BOOST_TEST(std::signbit(static_cast<double>("-1234567890123456789012345678901234e-400"_DL)));
    BOOST_TEST_EQ(static_cast<float>("1234567890123456789012345678901234e-110"_DL), 0.0F);
}

// A tie rounds to the even value
void test_ties()
{
    BOOST_TEST_EQ(static_cast<double>("4e23"_DF), 4e23);
    BOOST_TEST_EQ(static_cast<double>("9007199254740993"_DD), 9007199254740992.0);
    BOOST_TEST_EQ(static_cast<double>("9007199254740995"_DD), 9007199254740996.0);
    BOOST_TEST_EQ(static_cast<float>("16777217"_DD), 16777216.0F);
    BOOST_TEST_EQ(static_cast<float>("16777219"_DD), 16777220.0F);
    // Next to the midpoint 1 + 2^-24, where the product with the rounded 10^-15 cannot tell the side
    BOOST_TEST_EQ(static_cast<float>("1000000059604645e-15"_DD), 1.00000011920928955078125F);
    BOOST_TEST_EQ(static_cast<float>("1000000059604644e-15"_DD), 1.0F);

    // (2^53 + 1) / 2^25 has 34 digits and is halfway between 2^28 and the next double
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953125e-25"_DL), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953126e-25"_DL), 268435456.00000006);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953124e-25"_DL), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953125e-25"_DLF), 268435456.0);
    BOOST_TEST_EQ(static_cast<double>("9007199254740991.9"_DL), 9007199254740992.0);
    BOOST_TEST_EQ(static_cast<float>("0.5"_DF), 0.5F);
    BOOST_TEST_EQ(static_cast<double>("0.1"_DD), 0.1);

    // The decimal rounding mode does not change the conversion: 0.3 rounds down and -0.3 up to the nearest double
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(static_cast<double>("0.3"_DD), 0.3);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953124e-25"_DL), 268435456.0);
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(static_cast<double>("-0.3"_DD), -0.3);
    BOOST_TEST_EQ(static_cast<double>("2684354560000000298023223876953126e-25"_DL), 268435456.00000006);
    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

// The 19 leading digits of these values do not decide the rounding, so an exact compare does
void test_wide_significand()
{
    BOOST_TEST_EQ(static_cast<double>("1234567890123456789012345678901234e-20"_DL), 12345678901234.568);
    // The top 64 bits are all ones, and the value is just below the midpoint 5^23 * 2^135
    BOOST_TEST_EQ(static_cast<double>("5192296858534827628530496329220095e23"_DL), 5.192296858534827e+56);
    BOOST_TEST_EQ(static_cast<double>("1112856818313216103296202118680802e265"_DL), 1.112856818313216e+298);
    BOOST_TEST_EQ(static_cast<double>("1193075813432793838432928155755821e-308"_DL), 1.1930758134327938e-275);
    BOOST_TEST_EQ(static_cast<double>("3880824273833896262359420750490712e-35"_DL), 0.03880824273833896);
    BOOST_TEST_EQ(static_cast<double>("6099573892221131429277060931893378e-351"_DL), 6.09957e-318);
    BOOST_TEST_EQ(static_cast<double>("6099573892221131429277060931893379e-351"_DL), 6.099576e-318);
    BOOST_TEST_EQ(static_cast<float>("1058826583522949505694339754652399e-63"_DL), 1.05882654e-30F);
    BOOST_TEST_EQ(static_cast<float>("1058826583522949505694339754652400e-63"_DL), 1.05882663e-30F);
    BOOST_TEST_EQ(static_cast<float>("1175494420887210724209590083408724e-71"_DL), 1.17549435e-38F);
    BOOST_TEST_EQ(static_cast<float>("1175494420887210724209590083408725e-71"_DLF), 1.17549449e-38F);
    BOOST_TEST_EQ(static_cast<float>("7006492321624085354618647916449580e-79"_DL), 0.0F);
    BOOST_TEST_EQ(static_cast<float>("7006492321624085354618647916449581e-79"_DL), 1e-45F);
    BOOST_TEST_EQ(static_cast<float>("3402820728259992105141534601966245e5"_DL), 3.40282063e+38F);
    BOOST_TEST_EQ(static_cast<float>("3402820728259992105141534601966246e5"_DLF), 3.40282083e+38F);
}

// The coefficient of a short value of a fast type ends in zeros, and the last value has only its low bits clear
void test_fast_low_zeros()
{
    BOOST_TEST_EQ(static_cast<float>(19.99_DLF), 19.99F);
    BOOST_TEST_EQ(static_cast<double>(19.99_DLF), 19.99);
    BOOST_TEST_EQ(static_cast<double>(1e30_DLF), 1e30);
    BOOST_TEST_EQ(static_cast<double>(1e-320_DLF), 1e-320);
    BOOST_TEST_EQ(static_cast<double>(9007199254740993_DLF), 9007199254740992.0);
    BOOST_TEST_EQ(static_cast<double>("1267650600228229401496703205376"_DLF), 1267650600228229401496703205376.0);

    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    BOOST_TEST_EQ(static_cast<long double>(19.99_DLF), 19.99L);
    BOOST_TEST_EQ(static_cast<long double>("1267650600228229401496703205376"_DLF), 1267650600228229401496703205376.0L);
    #endif
}

// to_float rounds a C++23 type in its own format, and the 16-bit types round to infinity at the midpoint above the
// largest value. A tie there goes to infinity, as the largest value has an odd significand.
void test_extended_types()
{
    #ifdef BOOST_DECIMAL_HAS_FLOAT32
    BOOST_TEST_EQ((to_float<decimal128_t, std::float32_t>("1.000000059604644775390626"_DL)), 0x1.000002p+0f32);
    BOOST_TEST_EQ(static_cast<std::float32_t>("1.000000059604644775390626"_DL), 0x1.000002p+0f32);
    #endif

    #ifdef BOOST_DECIMAL_HAS_FLOAT64
    BOOST_TEST_EQ((to_float<decimal128_t, std::float64_t>("1e-320"_DL)), static_cast<std::float64_t>(1e-320));
    #endif

    // binary128 also rounds where long double has another format. The bits compare, as GCC can not print the type.
    #ifdef BOOST_DECIMAL_HAS_STDFLOAT128
    const auto same {[](const std::float128_t a, const std::float128_t b)
                     { return std::memcmp(&a, &b, sizeof(a)) == 0; }};
    BOOST_TEST(same(to_float<decimal128_t, std::float128_t>("0.1"_DL), 0.1f128));
    BOOST_TEST(same(to_float<decimal128_t, std::float128_t>("-1.189731495357231765085759326628007e4932"_DL),
                    -std::numeric_limits<std::float128_t>::max()));
    BOOST_TEST(same(to_float<decimal128_t, std::float128_t>("1e4933"_DL),
                    std::numeric_limits<std::float128_t>::infinity()));
    BOOST_TEST(same(to_float<decimal128_t, std::float128_t>("6.5e-4966"_DL),
                    std::numeric_limits<std::float128_t>::denorm_min()));
    #endif

    #if defined(BOOST_DECIMAL_HAS_FLOAT16) && !defined(BOOST_DECIMAL_FAST_MATH)
    constexpr auto float16_inf {std::numeric_limits<std::float16_t>::infinity()};
    BOOST_TEST_EQ((to_float<decimal64_t, std::float16_t>(65519.99_DD)), static_cast<std::float16_t>(65504));
    BOOST_TEST_EQ((to_float<decimal64_t, std::float16_t>(65520_DD)), float16_inf);
    BOOST_TEST_EQ((to_float<decimal64_t, std::float16_t>(-65520_DD)), -float16_inf);
    BOOST_TEST_EQ(static_cast<std::float16_t>(65520.01_DD), float16_inf);
    #endif

    #if defined(BOOST_DECIMAL_HAS_BRAINFLOAT16) && !defined(BOOST_DECIMAL_FAST_MATH)
    constexpr auto bfloat16_max {(std::numeric_limits<std::bfloat16_t>::max)()};
    BOOST_TEST_EQ((to_float<decimal128_t, std::bfloat16_t>("3396177529230460055269227039016280e5"_DL)), bfloat16_max);
    BOOST_TEST_EQ(static_cast<std::bfloat16_t>("3396177529230460055269227039016281e5"_DL),
                  std::numeric_limits<std::bfloat16_t>::infinity());
    #endif
}

// __float128 rounds in binary128 also where numeric_limits has nothing for it, as with libc++
#if defined(BOOST_DECIMAL_HAS_FLOAT128) && !defined(BOOST_DECIMAL_LDBL_IS_FLOAT128)
auto float128_high(const __float128 x) -> std::uint64_t
{
    std::uint64_t words[2] {};
    std::memcpy(words, &x, sizeof(words));
    return words[BOOST_DECIMAL_ENDIAN_LITTLE_BYTE ? 1 : 0];
}

void test_float128()
{
    const decimal128_t large {1, 400};
    const auto x {to_float<decimal128_t, __float128>(large)};
    BOOST_TEST_EQ(decimal128_t {x}, large);
    BOOST_TEST_EQ(float128_high(x) >> 48U, UINT64_C(0x452F));
    BOOST_TEST_EQ(float128_high(to_float<decimal128_t, __float128>(-std::numeric_limits<decimal128_t>::infinity())),
                  UINT64_C(0xFFFF000000000000));
    BOOST_TEST_EQ(float128_high(to_float<decimal128_t, __float128>(std::numeric_limits<decimal128_t>::quiet_NaN())),
                  UINT64_C(0x7FFF800000000000));
    BOOST_TEST_EQ(float128_high(to_float<decimal128_t, __float128>(decimal128_t {1, 4933})), UINT64_C(0x7FFF000000000000));
}
#endif

// A raw encoding with a coefficient of 10^34 is not canonical, and IEEE 754 reads it as zero
void test_noncanonical()
{
    const auto x {from_bits(decimal128_t::significand_type {UINT64_C(0x3041ED09BEAD87C0), UINT64_C(0x378D8E6400000000)})};
    const auto positive {static_cast<double>(x)};
    const auto negative {static_cast<double>(-x)};
    BOOST_TEST(std::fpclassify(positive) == FP_ZERO && !std::signbit(positive));
    BOOST_TEST(std::fpclassify(negative) == FP_ZERO && std::signbit(negative));
    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    BOOST_TEST(std::fpclassify(static_cast<long double>(x)) == FP_ZERO);
    #endif
    BOOST_TEST(std::fpclassify(static_cast<float>(-x)) == FP_ZERO && std::signbit(static_cast<float>(-x)));
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
    test_fast_low_zeros();
    test_extended_types();
    test_noncanonical();
    #if defined(BOOST_DECIMAL_HAS_FLOAT128) && !defined(BOOST_DECIMAL_LDBL_IS_FLOAT128)
    test_float128();
    #endif

    return boost::report_errors();
}
