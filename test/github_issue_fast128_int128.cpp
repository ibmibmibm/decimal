// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

using namespace boost::decimal::literals;

int main()
{
    constexpr std::uint64_t pow10_19 {UINT64_C(10000000000000000000)};
    const boost::int128::uint128_t expected {boost::int128::uint128_t{pow10_19} * pow10_19};
    const boost::int128::int128_t signed_expected {expected};

    const auto value {1e38_DLF};

    BOOST_TEST(static_cast<boost::int128::uint128_t>(value) == expected);
    BOOST_TEST(static_cast<boost::int128::int128_t>(value) == signed_expected);

    // std::numeric_limits of an old libstdc++ knows __int128 only with the GNU extensions
    #if defined(BOOST_DECIMAL_HAS_INT128) && !defined(__STRICT_ANSI__)

    BOOST_TEST(static_cast<boost::decimal::detail::builtin_uint128_t>(value) == static_cast<boost::decimal::detail::builtin_uint128_t>(expected));
    BOOST_TEST(static_cast<boost::decimal::detail::builtin_int128_t>(value) == static_cast<boost::decimal::detail::builtin_int128_t>(signed_expected));

    #endif

    return boost::report_errors();
}
