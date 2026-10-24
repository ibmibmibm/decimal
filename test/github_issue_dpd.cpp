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

// The DPD encodings come from the Intel Decimal Floating-Point Math Library
void test_encodings()
{
    BOOST_TEST_EQ(to_dpd(std::numeric_limits<decimal32_t>::max()), UINT32_C(0x77F3FCFF));
    BOOST_TEST_EQ(to_dpd(1234.567_DF), UINT32_C(0x2624D2E7));
    BOOST_TEST_EQ(to_dpd(std::numeric_limits<decimal64_t>::max()), UINT64_C(0x77FCFF3FCFF3FCFF));
    BOOST_TEST_EQ(to_dpd(-std::numeric_limits<decimal64_t>::denorm_min()), UINT64_C(0x8000000000000001));
    BOOST_TEST(to_dpd(std::numeric_limits<decimal128_t>::max()) == (boost::int128::uint128_t{UINT64_C(0x77FFCFF3FCFF3FCF), UINT64_C(0xF3FCFF3FCFF3FCFF)}));
    BOOST_TEST(to_dpd(8E6111_DL) == (boost::int128::uint128_t{UINT64_C(0x43FFC00000000000), UINT64_C(8)}));

    BOOST_TEST_EQ(from_dpd<decimal32_t>(UINT32_C(0x2624D2E7)), 1234.567_DF);
    BOOST_TEST_EQ(from_dpd<decimal64_t>(UINT64_C(0x77FCFF3FCFF3FCFF)), std::numeric_limits<decimal64_t>::max());
    BOOST_TEST_EQ(from_dpd<decimal128_t>(boost::int128::uint128_t{UINT64_C(0x43FFC00000000000), UINT64_C(8)}), 8E6111_DL);
}

// The payload of a NaN is in declets, as the digits of a finite value are (IEEE 754 3.5.2)
template <typename T, typename Bits>
void test_nan(const char* payload, const Bits quiet, const Bits negative_signaling)
{
    const auto q {nan<T>(payload)};
    const auto s {-snan<T>(payload)};
    BOOST_TEST(to_dpd(q) == quiet);
    BOOST_TEST(to_dpd(s) == negative_signaling);

    const auto q2 {from_dpd<T>(quiet)};
    const auto s2 {from_dpd<T>(negative_signaling)};
    BOOST_TEST(isnan(q2) && !issignaling(q2) && !signbit(q2));
    BOOST_TEST(issignaling(s2) && signbit(s2));
    BOOST_TEST_EQ(read_payload(q2), read_payload(q));
    BOOST_TEST_EQ(read_payload(s2), read_payload(s));
}

// An infinity and a noncanonical NaN payload give the canonical encoding, and from_dpd reads every declet
void test_noncanonical()
{
    BOOST_TEST_EQ(to_dpd(from_bid<decimal32_t>(UINT32_C(0xF9FFFFFF))), UINT32_C(0xF8000000));
    BOOST_TEST_EQ(to_dpd(from_bid<decimal32_t>(UINT32_C(0xFE0F4240))), UINT32_C(0xFE000000));
    BOOST_TEST_EQ(to_dpd(from_bid<decimal64_t>(UINT64_C(0x7C038D7EA4C68000))), UINT64_C(0x7C00000000000000));
    BOOST_TEST(to_dpd(from_bid<decimal128_t>(boost::int128::uint128_t{UINT64_C(0x7C00314DC6448D93), UINT64_C(0x38C15B0A00000000)})) ==
               (boost::int128::uint128_t{UINT64_C(0x7C00000000000000), UINT64_C(0)}));

    const auto alternate {from_dpd<decimal32_t>(UINT32_C(0x7C0003FF))};
    BOOST_TEST_EQ(read_payload(alternate), 999U);
    BOOST_TEST_EQ(to_dpd(alternate), UINT32_C(0x7C0000FF));
}

// The 24 alternate declets decode to the 8 numbers whose digits are all 8 or 9
void test_alternate_declets()
{
    const auto zero {to_dpd(0_DF)};
    int alternates {};
    for (std::uint32_t declet {}; declet < 1024U; ++declet)
    {
        const auto val {from_dpd<decimal32_t>(zero | declet)};
        const auto canonical {to_dpd(val)};
        BOOST_TEST_EQ(to_bid(from_dpd<decimal32_t>(canonical)), to_bid(val));
        if (canonical != (zero | declet))
        {
            ++alternates;
            const auto n {static_cast<unsigned>(val)};
            BOOST_TEST(n / 100U >= 8U && n / 10U % 10U >= 8U && n % 10U >= 8U);
        }
    }
    BOOST_TEST_EQ(alternates, 24);
}

// A round trip keeps the sign of zero and the exponent, which numeric equality does not show
template <typename T>
void test_representation(const T one, const T one_point_zero)
{
    const T negative_zero {-T{0}};
    BOOST_TEST(to_bid(from_dpd<T>(to_dpd(negative_zero))) == to_bid(negative_zero));
    BOOST_TEST(to_bid(from_dpd<T>(to_dpd(one))) == to_bid(one));
    BOOST_TEST(to_bid(from_dpd<T>(to_dpd(one_point_zero))) == to_bid(one_point_zero));
    BOOST_TEST(to_dpd(one) != to_dpd(one_point_zero));
}

int main()
{
    test_encodings();
    test_alternate_declets();
    test_representation(1.0_DF, 1.00_DF);
    test_representation(1.0_DD, 1.00_DD);
    test_representation(1.0_DL, 1.00_DL);

    #ifndef BOOST_DECIMAL_FAST_MATH
    test_noncanonical();
    test_nan<decimal32_t>("123", UINT32_C(0x7C0000A3), UINT32_C(0xFE0000A3));
    test_nan<decimal32_t>("999999", UINT32_C(0x7C03FCFF), UINT32_C(0xFE03FCFF));
    test_nan<decimal64_t>("123", UINT64_C(0x7C000000000000A3), UINT64_C(0xFE000000000000A3));
    test_nan<decimal64_t>("999999999999999", UINT64_C(0x7C00FF3FCFF3FCFF), UINT64_C(0xFE00FF3FCFF3FCFF));
    test_nan<decimal128_t>("123", boost::int128::uint128_t{UINT64_C(0x7C00000000000000), UINT64_C(0xA3)},
                           boost::int128::uint128_t{UINT64_C(0xFE00000000000000), UINT64_C(0xA3)});
    test_nan<decimal128_t>("456123000000000000789",
                           boost::int128::uint128_t{UINT64_C(0x7C00000000000025), UINT64_C(0x628C0000000003CF)},
                           boost::int128::uint128_t{UINT64_C(0xFE00000000000025), UINT64_C(0x628C0000000003CF)});
    test_nan<decimal128_t>("999999999999999999999999999999999",
                           boost::int128::uint128_t{UINT64_C(0x7C000FF3FCFF3FCF), UINT64_C(0xF3FCFF3FCFF3FCFF)},
                           boost::int128::uint128_t{UINT64_C(0xFE000FF3FCFF3FCF), UINT64_C(0xF3FCFF3FCFF3FCFF)});
    #endif

    return boost::report_errors();
}
