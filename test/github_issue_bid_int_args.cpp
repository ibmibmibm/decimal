// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// from_bid and from_dpd take any integer of 32 or 64 bits, such as unsigned long long

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <type_traits>
#include <utility>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The overloads for other integers take only a target that the conversion supports
template <typename T, typename Integer, typename = void>
struct has_from_bid : std::false_type {};

template <typename T, typename Integer>
struct has_from_bid<T, Integer, decltype(static_cast<void>(from_bid<T>(std::declval<Integer>())))> : std::true_type {};

template <typename T, typename Integer, typename = void>
struct has_from_dpd : std::false_type {};

template <typename T, typename Integer>
struct has_from_dpd<T, Integer, decltype(static_cast<void>(from_dpd<T>(std::declval<Integer>())))> : std::true_type {};

static_assert(!has_from_bid<void, int>::value && !has_from_bid<void, long long>::value, "no decimal target");
static_assert(!has_from_dpd<void, int>::value && !has_from_dpd<void, long long>::value, "no decimal target");
static_assert(!has_from_dpd<decimal64_t, int>::value && !has_from_dpd<decimal128_t, long long>::value, "other width");
static_assert(has_from_bid<decimal128_t, long long>::value && has_from_dpd<decimal_fast32_t, int>::value, "supported");

int main()
{
    const unsigned long long bid64 {0x31C0000000000001ULL};
    const unsigned long long dpd64 {0x2238000000000001ULL};

    static_assert(std::is_same<decltype(from_bid(bid64)), decimal64_t>::value, "default type");
    static_assert(std::is_same<decltype(from_dpd(dpd64)), decimal64_t>::value, "default type");
    static_assert(std::is_same<decltype(from_bid(0x32800001)), decimal32_t>::value, "default type");
    static_assert(std::is_same<decltype(from_dpd(0x22500001)), decimal32_t>::value, "default type");

    BOOST_TEST_EQ(from_bid(bid64), 1_DD);
    BOOST_TEST_EQ(from_bid<decimal64_t>(bid64), 1_DD);
    BOOST_TEST_EQ(from_bid<decimal_fast64_t>(bid64), 1_DDF);
    BOOST_TEST_EQ(from_bid<decimal128_t>(bid64), 1_DL);
    BOOST_TEST_EQ(from_bid<decimal64_t>(0x31C0000000000002LL), 2_DD);
    BOOST_TEST_EQ(from_bid(0x32800001), 1_DF);
    BOOST_TEST_EQ(from_bid<decimal_fast32_t>(0x32800003), 3_DFF);

    BOOST_TEST_EQ(from_dpd(dpd64), 1_DD);
    BOOST_TEST_EQ(from_dpd<decimal64_t>(dpd64), 1_DD);
    BOOST_TEST_EQ(from_dpd<decimal_fast64_t>(dpd64), 1_DDF);
    BOOST_TEST_EQ(from_dpd(0x22500001), 1_DF);
    BOOST_TEST_EQ(from_dpd<decimal_fast32_t>(0x22500001), 1_DFF);

    constexpr auto one {from_bid(0x31C0000000000001LL)};
    static_assert(one == 1_DD, "constexpr");

    return boost::report_errors();
}
