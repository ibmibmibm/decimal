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

#ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

// A direction that rounds to x in the decimal type is still above or below x
void test_direction_rounds_to_x()
{
    BOOST_TEST_EQ(nexttoward("1"_DF, 1.0L + 1e-12L), "1.000001"_DF);
    BOOST_TEST_EQ(nexttoward("1"_DF, 1.0L - 1e-12L), "0.9999999"_DF);
    BOOST_TEST_EQ(nexttoward("1"_DD, 1.0L + 4e-16L), "1.000000000000001"_DD);
    BOOST_TEST_EQ(nexttoward("1"_DFF, 1.0L - 1e-12L), "0.9999999"_DFF);
    BOOST_TEST_EQ(nexttoward("1e10"_DF, 1e10L - 1e-3L), "9.999999e9"_DF);

    // 10^-25 and 10^20 times the significand of x
    BOOST_TEST_EQ(nexttoward("1e-10"_DD, 1e-10L + 4e-26L), "1.000000000000001e-10"_DD);
    BOOST_TEST_EQ(nexttoward("1e-10"_DD, 1e-10L - 4e-26L), "9.999999999999999e-11"_DD);
    BOOST_TEST_EQ(nexttoward("1e35"_DD, 1e35L + 4e19L), "1.000000000000001e35"_DD);

    // 10^34 times the significand of x, which is too large for the compare in 128 bits
    BOOST_TEST_EQ(nexttoward("1e40"_DF, 1e40L + 1e28L), "1.000001e40"_DF);
    BOOST_TEST_EQ(nexttoward("1e40"_DF, 1e40L - 1e28L), "9.999999e39"_DF);

    // The binary value of 0.1L is above 0.1
    BOOST_TEST_EQ(nexttoward("0.1"_DF, 0.1L), "0.1000001"_DF);
    BOOST_TEST_EQ(nexttoward("-0.1"_DF, -0.1L), "-0.1000001"_DF);

    // The binary values are next to x, where 10^-28 and 10^28 times its significand end the compare in 128 bits
    #if BOOST_DECIMAL_LDBL_BITS == 80
    BOOST_TEST_EQ(nexttoward("6e-22"_DF, 6e-22L), "6.000001e-22"_DF);
    BOOST_TEST_EQ(nexttoward("3e34"_DF, 3e34L), "2.999999e34"_DF);
    BOOST_TEST_EQ(nexttoward("-3.3e44"_DLF, -3.3e44L), "-3.299999999999999999999999999999999e44"_DLF);
    #endif

    // A constant expression reads the bits with bit_cast, which Clang does not do for a long double
    #if defined(BOOST_DECIMAL_HAS_CONSTEXPR_BITCAST) && !defined(BOOST_DECIMAL_NO_CONSTEVAL_DETECTION)
    static_assert(nexttoward("1"_DF, 1.0L + 1e-12L) == "1.000001"_DF, "the direction is above 1");
    #endif
}

// The shortest decimal string of a long double can be on the other side of x
void test_shortest_string_other_side()
{
    #if BOOST_DECIMAL_LDBL_BITS == 80
    // y = 1 - 2^-64 is 0.99999999999999999994578989137572477...
    const long double y {1.0L - 1.0L / 18446744073709551616.0L};
    BOOST_TEST_EQ(nexttoward("0.9999999999999999999457898913757248"_DL, y), "0.9999999999999999999457898913757247"_DL);
    BOOST_TEST_EQ(nexttoward("0.9999999999999999999457898913757247"_DL, y), "0.9999999999999999999457898913757248"_DL);
    #endif
}

// An equal direction gives x, with the sign of direction for a zero
void test_equal()
{
    BOOST_TEST_EQ(nexttoward("0.5"_DF, 0.5L), "0.5"_DF);

    const auto neg_zero {nexttoward("0"_DF, -0.0L)};
    BOOST_TEST_EQ(neg_zero, "0"_DF);
    BOOST_TEST(signbit(neg_zero));

    BOOST_TEST_EQ(nexttoward("0"_DD, 1e-300L), std::numeric_limits<decimal64_t>::denorm_min());
    BOOST_TEST_EQ(nexttoward(std::numeric_limits<decimal32_t>::denorm_min(), 0.0L), "0"_DF);
}

// A direction of the other sign, and a direction whose exact compare needs the large big integers
void test_far_direction()
{
    BOOST_TEST_EQ(nexttoward("-1"_DF, 1.0L), "-0.9999999"_DF);
    BOOST_TEST_EQ(nexttoward("1"_DF, -1.0L), "0.9999999"_DF);
    BOOST_TEST_EQ(nexttoward("-2e23"_DF, -2.3e22L), "-1.999999e23"_DF);
    BOOST_TEST_EQ(nexttoward("-5e3"_DF, -2.2e6L), "-5.000001e3"_DF);
    BOOST_TEST_EQ(nexttoward("1e-101"_DF, std::numeric_limits<long double>::denorm_min()), "0"_DF);

    #if BOOST_DECIMAL_LDBL_BITS == 80
    BOOST_TEST_EQ(nexttoward("1e4000"_DL, 1e4000L), "9.999999999999999999999999999999999e3999"_DL);
    // Subnormal directions, whose exponent is that of the smallest normal number
    BOOST_TEST_EQ(nexttoward("3e-4937"_DL, 3e-4937L), "3.000000000000000000000000000000001e-4937"_DL);
    BOOST_TEST_EQ(nexttoward("3.18725502e-4949"_DL, 87 * std::numeric_limits<long double>::denorm_min()),
                  "3.187255019999999999999999999999999e-4949"_DL);
    #elif BOOST_DECIMAL_LDBL_BITS == 128
    BOOST_TEST_EQ(nexttoward("1e4000"_DL, 1e4000L), "1.000000000000000000000000000000001e4000"_DL);
    #endif
}

// The neighbor of a normal x toward to, with exact arithmetic only: GCC with -mlong-double-64 or -mlong-double-128
// keeps the library functions of the host format. A step toward zero from a power of two is half as large.
long double long_double_step(const long double x, const long double to)
{
    const long double magnitude {x < 0 ? -x : x};
    long double power {1};
    while (power > magnitude)
    {
        power /= 2;
    }
    while (power * 2 <= magnitude)
    {
        power *= 2;
    }
    const long double unit {power * std::numeric_limits<long double>::epsilon()};
    const bool inward {(to > x) != (x > 0)};
    const long double step {inward && !(magnitude > power) ? unit / 2 : unit};
    return to > x ? x + step : x - step;
}

// The side of x on which the nearest long double lies, from exact arithmetic. Its neighbor on the other side of x is a
// step the other way.
template <typename T>
void check_side(const T x, const long double nearest, const int side)
{
    const T up {nextafter(x, std::numeric_limits<T>::infinity())};
    const T down {nextafter(x, -std::numeric_limits<T>::infinity())};
    BOOST_TEST_EQ(nexttoward(x, nearest), side > 0 ? up : side < 0 ? down : x);
    BOOST_TEST_EQ(nexttoward(-x, -nearest), side > 0 ? -up : side < 0 ? -down : -x);
    if (side >= 0)
    {
        BOOST_TEST_EQ(nexttoward(x, long_double_step(nearest, -std::numeric_limits<long double>::infinity())), down);
    }
    if (side <= 0)
    {
        BOOST_TEST_EQ(nexttoward(x, long_double_step(nearest, std::numeric_limits<long double>::infinity())), up);
    }
}

// Values next to the exponent bounds of the exact compares: 27 for 128 bits, and 550 for the 128-bit decimal types.
// The side of each value is for the x87, binary128 and binary64 formats of long double, in this order.
void test_compare_bounds()
{
    #if BOOST_DECIMAL_LDBL_BITS == 80
    constexpr std::size_t format {0};
    #elif BOOST_DECIMAL_LDBL_BITS == 128
    constexpr std::size_t format {1};
    #else
    constexpr std::size_t format {2};
    #endif

    struct narrow_case { decimal64_t x; long double nearest; int side[3]; };
    const narrow_case narrow_cases[] {
        {"1234567890123456e26"_DD, 1234567890123456e26L, {-1, 0, 1}},
        {"1234567890123456e27"_DD, 1234567890123456e27L, {-1, 0, -1}},
        {"1234567890123456e28"_DD, 1234567890123456e28L, {1, 0, 1}},
        {"1234567890123456e-26"_DD, 1234567890123456e-26L, {-1, -1, 1}},
        {"1234567890123456e-27"_DD, 1234567890123456e-27L, {-1, 1, 1}},
        {"1234567890123456e-28"_DD, 1234567890123456e-28L, {-1, -1, 1}},
        {"9999999999999999e26"_DD, 9999999999999999e26L, {-1, 1, -1}},
        {"9999999999999999e27"_DD, 9999999999999999e27L, {-1, -1, -1}},
        {"9999999999999999e28"_DD, 9999999999999999e28L, {1, 1, -1}},
        {"9999999999999999e-26"_DD, 9999999999999999e-26L, {-1, 1, 1}},
        {"9999999999999999e-27"_DD, 9999999999999999e-27L, {-1, -1, 1}},
        {"9999999999999999e-28"_DD, 9999999999999999e-28L, {1, -1, 1}},
        {"1000000000000000e26"_DD, 1000000000000000e26L, {-1, 0, 1}},
        {"1000000000000000e27"_DD, 1000000000000000e27L, {1, 0, 1}},
        {"1000000000000000e28"_DD, 1000000000000000e28L, {-1, 0, 1}},
        {"1000000000000000e-26"_DD, 1000000000000000e-26L, {1, 1, -1}},
        {"1000000000000000e-27"_DD, 1000000000000000e-27L, {1, 1, -1}},
        {"1000000000000000e-28"_DD, 1000000000000000e-28L, {-1, 1, 1}},
    };
    for (const auto& c : narrow_cases)
    {
        check_side(c.x, c.nearest, c.side[format]);
    }

    // Binary64 cannot hold these directions
    #if BOOST_DECIMAL_LDBL_BITS != 64
    struct wide_case { decimal128_t x; long double nearest; int side[2]; };
    const wide_case wide_cases[] {
        {"1234567890123456789012345678901234e549"_DL, 1234567890123456789012345678901234e549L, {-1, 1}},
        {"1234567890123456789012345678901234e550"_DL, 1234567890123456789012345678901234e550L, {1, 1}},
        {"1234567890123456789012345678901234e551"_DL, 1234567890123456789012345678901234e551L, {-1, -1}},
        {"1234567890123456789012345678901234e-549"_DL, 1234567890123456789012345678901234e-549L, {-1, -1}},
        {"1234567890123456789012345678901234e-550"_DL, 1234567890123456789012345678901234e-550L, {-1, -1}},
        {"1234567890123456789012345678901234e-551"_DL, 1234567890123456789012345678901234e-551L, {1, -1}},
        {"9999999999999999999999999999999999e549"_DL, 9999999999999999999999999999999999e549L, {1, 1}},
        {"9999999999999999999999999999999999e550"_DL, 9999999999999999999999999999999999e550L, {1, 1}},
        {"9999999999999999999999999999999999e551"_DL, 9999999999999999999999999999999999e551L, {1, -1}},
        {"9999999999999999999999999999999999e-549"_DL, 9999999999999999999999999999999999e-549L, {1, 1}},
        {"9999999999999999999999999999999999e-550"_DL, 9999999999999999999999999999999999e-550L, {1, -1}},
        {"9999999999999999999999999999999999e-551"_DL, 9999999999999999999999999999999999e-551L, {1, 1}},
        {"1000000000000000000000000000000000e549"_DL, 1000000000000000000000000000000000e549L, {-1, 1}},
        {"1000000000000000000000000000000000e550"_DL, 1000000000000000000000000000000000e550L, {1, -1}},
        {"1000000000000000000000000000000000e551"_DL, 1000000000000000000000000000000000e551L, {1, -1}},
        {"1000000000000000000000000000000000e-549"_DL, 1000000000000000000000000000000000e-549L, {1, -1}},
        {"1000000000000000000000000000000000e-550"_DL, 1000000000000000000000000000000000e-550L, {1, 1}},
        {"1000000000000000000000000000000000e-551"_DL, 1000000000000000000000000000000000e-551L, {-1, 1}},
    };
    for (const auto& c : wide_cases)
    {
        check_side(c.x, c.nearest, c.side[format]);
    }
    #endif
}

void test_non_finite()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto inf {std::numeric_limits<long double>::infinity()};

    BOOST_TEST(isinf(nexttoward((std::numeric_limits<decimal32_t>::max)(), inf)));
    BOOST_TEST_EQ(nexttoward(std::numeric_limits<decimal32_t>::infinity(), 1.0L), (std::numeric_limits<decimal32_t>::max)());
    BOOST_TEST(isinf(nexttoward(std::numeric_limits<decimal32_t>::infinity(), inf)));
    BOOST_TEST_EQ(nexttoward("1"_DF, -inf), "0.9999999"_DF);
    BOOST_TEST(isnan(nexttoward("1"_DF, std::numeric_limits<long double>::quiet_NaN())));

    #if BOOST_DECIMAL_LDBL_BITS == 80
    // An unnormal number, with the leading bit clear, is a NaN as for the hardware
    const unsigned char unnormal_bytes[10] {0, 0, 0, 0, 0, 0, 0, 0x40, 0xFF, 0x3F};
    long double unnormal {};
    std::memcpy(&unnormal, unnormal_bytes, sizeof(unnormal_bytes));
    BOOST_TEST(isnan(nexttoward("1"_DF, unnormal)));
    #endif

    // A decimal128_t value can be beyond the long double range
    BOOST_TEST_EQ(nexttoward("1e5000"_DL, inf), "1.000000000000000000000000000000001e5000"_DL);
    BOOST_TEST_EQ(nexttoward("1e5000"_DL, 1.0L), "9.999999999999999999999999999999999e4999"_DL);
    #endif
}

int main()
{
    test_direction_rounds_to_x();
    test_shortest_string_other_side();
    test_equal();
    test_far_direction();
    test_compare_bounds();
    test_non_finite();

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
