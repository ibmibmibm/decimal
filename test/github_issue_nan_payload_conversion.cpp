// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// A conversion keeps the sign and the payload of a NaN (IEEE 754 6.2.3)
template <typename To, typename From>
void test_conversion()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    constexpr bool same_precision {std::numeric_limits<To>::digits10 == std::numeric_limits<From>::digits10};

    for (const auto x : {nan<From>("123"), -nan<From>("123"), snan<From>("45"), -snan<From>("45")})
    {
        const auto y {static_cast<To>(x)};

        BOOST_TEST(isnan(y));
        BOOST_TEST_EQ(signbit(y), signbit(x));
        BOOST_TEST_EQ(read_payload(y), static_cast<typename To::significand_type>(read_payload(x)));
        BOOST_TEST_EQ(issignaling(y), same_precision && issignaling(x));
    }
    #endif
}

template <typename From>
void test_all_targets()
{
    test_conversion<decimal32_t, From>();
    test_conversion<decimal64_t, From>();
    test_conversion<decimal128_t, From>();
    test_conversion<decimal_fast32_t, From>();
    test_conversion<decimal_fast64_t, From>();
    test_conversion<decimal_fast128_t, From>();
}

// A payload that the target cannot hold gives the NaN without a payload
void test_wide_payload()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto x {static_cast<decimal32_t>(nand128("100000000000000000000"))};

    BOOST_TEST(isnan(x));
    BOOST_TEST_EQ(read_payload(x), 0U);

    const auto y {static_cast<decimal_fast64_t>(nand128f("100000000000000000000"))};

    BOOST_TEST(isnan(y));
    BOOST_TEST_EQ(read_payload(y), 0U);
    #endif
}

// The NaN operand of an operation between two types gives the payload of the result
void test_mixed_arithmetic()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(read_payload(nand32("5") + decimal64_t{1}), 5U);
    BOOST_TEST_EQ(read_payload(decimal64_t{1} - nand32("5")), 5U);
    BOOST_TEST_EQ(read_payload(nand64("5") * decimal128_t{2}), 5U);
    BOOST_TEST_EQ(read_payload(decimal_fast128_t{2} / nand32f("5")), 5U);
    #endif
}

// The largest payload is 10^(p - 1) - 1, and the trailing significand field holds it (IEEE 754 3.5.2)
template <typename T, typename Bits>
void test_payload_limit(const Bits nan_bits, const Bits ignored_bit)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    using sig = typename T::significand_type;
    const auto largest {static_cast<sig>(boost::decimal::detail::pow10(static_cast<sig>(std::numeric_limits<T>::digits10 - 1)) - 1U)};

    BOOST_TEST_EQ(read_payload(boost::decimal::detail::write_payload<T, false>(largest)), largest);
    BOOST_TEST_EQ(read_payload(boost::decimal::detail::write_payload<T, false>(static_cast<sig>(largest + 1U))), 0U);
    BOOST_TEST_EQ(read_payload(from_bits(static_cast<Bits>(nan_bits | ignored_bit | 123U))), 123U);
    BOOST_TEST_EQ(read_payload(from_bits(static_cast<Bits>(nan_bits | (largest + 1U)))), 0U);
    #else
    static_cast<void>(nan_bits);
    static_cast<void>(ignored_bit);
    #endif
}

// A payload above the limit of a narrower target gives the NaN without a payload
void test_narrow_payload()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(read_payload(static_cast<decimal32_t>(nand64("999999"))), 999999U);
    BOOST_TEST_EQ(read_payload(static_cast<decimal32_t>(nand64("1000000"))), 0U);
    #endif
}

// A signaling NaN of an other type is selected as for operands of one type
void test_mixed_signaling()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto a {nand64("1") + -snand32("2")};
    BOOST_TEST(isnan(a) && !issignaling(a) && signbit(a) && read_payload(a) == 2U);

    const auto b {-snand32("2") + nand64("1")};
    BOOST_TEST(isnan(b) && !issignaling(b) && signbit(b) && read_payload(b) == 2U);
    #endif
}

int main()
{
    test_all_targets<decimal32_t>();
    test_all_targets<decimal64_t>();
    test_all_targets<decimal128_t>();
    test_all_targets<decimal_fast32_t>();
    test_all_targets<decimal_fast64_t>();
    test_all_targets<decimal_fast128_t>();

    test_wide_payload();
    test_mixed_arithmetic();
    test_payload_limit<decimal32_t>(UINT32_C(0x7C000000), UINT32_C(0x00100000));
    test_payload_limit<decimal64_t>(UINT64_C(0x7C00000000000000), UINT64_C(0x0004000000000000));
    test_payload_limit<decimal128_t>(boost::int128::uint128_t {UINT64_C(0x7C00000000000000), 0U}, boost::int128::uint128_t {UINT64_C(0x0000400000000000), 0U});
    test_narrow_payload();
    test_mixed_signaling();

    return boost::report_errors();
}
