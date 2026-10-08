// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// https://github.com/boostorg/decimal/issues/NNNN
//
// The conversion to an integer compared the value with the limits of the integer type converted
// to the decimal type, and these limits can round. It also used a 64-bit integer for 128-bit types.

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cerrno>
#include <cstdint>
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename Int, typename T>
void test_out_of_range(const T value, const Int largest = (std::numeric_limits<Int>::max)())
{
    errno = 0;
    BOOST_TEST(static_cast<Int>(value) == largest);
    BOOST_TEST_EQ(errno, ERANGE);
}

template <typename Int, typename T>
void test_in_range(const T value, const Int expected)
{
    errno = 0;
    BOOST_TEST(static_cast<Int>(value) == expected);
    BOOST_TEST_EQ(errno, 0);
}

// INT_MAX and LLONG_MAX round up in decimal32_t and decimal64_t, so the next value overflowed.
void test_rounded_limits()
{
    test_out_of_range<int>(2147484e3_DF);
    test_out_of_range<int>(-2147484e3_DF);
    test_out_of_range<int>(2147484e3_DFF);
    test_out_of_range<long long>(9223372036854776e3_DD);
    test_out_of_range<long long>(-9223372036854776e3_DD);
    test_out_of_range<long long>(9223372036854776e3_DDF);
    test_in_range<long long>(9223372036854775e3_DD, 9223372036854775000LL);
    test_in_range<int>(2147483e3_DF, 2147483000);
}

// A value between a limit and the next integer away from zero truncates to the limit.
void test_truncation_to_limits()
{
    test_in_range<int>(2147483647.5_DD, 2147483647);
    test_in_range<int>(-2147483648.9_DD, (std::numeric_limits<int>::min)());
    test_in_range<long long>(9223372036854775807.9_DL, (std::numeric_limits<long long>::max)());
    test_in_range<long long>(-9223372036854775808.5_DL, (std::numeric_limits<long long>::min)());
    test_in_range<long long>(-9223372036854775808.5_DLF, (std::numeric_limits<long long>::min)());
    test_in_range<unsigned>(4294967295.5_DD, 4294967295U);
    test_in_range<unsigned>(-0.5_DF, 0U);
    test_out_of_range<int>(2147483648_DD);
    test_out_of_range<int>(-2147483649_DD);
    test_out_of_range<unsigned>(-1_DF);
}

// The int128 class is not std::is_signed, so a negative value lost its sign.
void test_int128_class()
{
    test_in_range<boost::int128::int128_t>(-5_DL, boost::int128::int128_t {-5});
    test_in_range<boost::int128::int128_t>(3_DL, boost::int128::int128_t {3});
    test_in_range<boost::int128::int128_t>(-3e19_DL, -boost::int128::int128_t {boost::int128::uint128_t {UINT64_C(10000000000000000000)} * 3U});
    test_in_range<boost::int128::int128_t>(-1e30_DL, boost::int128::int128_t {-1000000000000000LL} * 1000000000000000LL);

    #ifdef BOOST_DECIMAL_HAS_INT128
    test_in_range<boost::int128::int128_t>(-5_DLF, boost::int128::int128_t {-5});
    test_in_range<boost::int128::int128_t>(3_DLF, boost::int128::int128_t {3});
    test_in_range<boost::int128::int128_t>(-3e19_DLF, -boost::int128::int128_t {boost::int128::uint128_t {UINT64_C(10000000000000000000)} * 3U});
    test_in_range<boost::int128::uint128_t>(1e30_DLF, boost::int128::uint128_t {1000000000000000U} * 1000000000000000U);
    #endif
}

// A 64-bit integer cannot hold 1e30, so the result was the largest value.
void test_builtin_int128()
{
    #ifdef BOOST_DECIMAL_HAS_INT128
    using i128 = detail::builtin_int128_t;
    using u128 = detail::builtin_uint128_t;
    const auto e30 {static_cast<i128>(1000000000000000LL) * 1000000000000000LL};

    test_in_range<i128>(1e30_DF, e30);
    test_in_range<i128>(-1e30_DF, -e30);
    test_in_range<i128>(1e30_DD, e30);
    test_in_range<i128>(-1e30_DDF, -e30);
    test_in_range<u128>(1e30_DD, static_cast<u128>(e30));
    test_in_range<i128>(-1e30_DL, -e30);
    test_in_range<u128>(1e30_DL, static_cast<u128>(e30));
    test_out_of_range<i128>(1701412e32_DF, BOOST_DECIMAL_INT128_MAX);
    test_out_of_range<u128>(3402824e32_DF, BOOST_DECIMAL_UINT128_MAX);
    test_out_of_range<u128>(3402823669209385e23_DD, BOOST_DECIMAL_UINT128_MAX);
    test_out_of_range<u128>(-1e30_DD, BOOST_DECIMAL_UINT128_MAX);
    #endif
}

// The smallest value of a signed type is in range, and bool has the largest value 1
template <typename Int, typename T>
void test_helper(const T value, const Int expected, const int expected_errno)
{
    errno = 0;
    BOOST_TEST((to_integral<T, Int>(value) == expected));
    BOOST_TEST_EQ(errno, expected_errno);
}

void test_narrow_types()
{
    test_helper<signed char>(-128.5_DF, static_cast<signed char>(-128), 0);
    test_helper<signed char>(-129_DF, static_cast<signed char>(127), ERANGE);
    test_helper<short>(-32768_DD, static_cast<short>(-32768), 0);
    test_helper<bool>(decimal64_t {1}, true, 0);
    test_helper<bool>(decimal64_t {2}, true, ERANGE);
}

// A zero is in range with every exponent
template <typename T>
void test_zero_with_exponent()
{
    test_in_range<int>(T{0, 30}, 0);
    test_in_range<unsigned>(T{0, 45}, 0U);
    test_in_range<long long>(-T{0, 90}, 0LL);
}

int main()
{
    test_zero_with_exponent<decimal32_t>();
    test_zero_with_exponent<decimal64_t>();
    test_zero_with_exponent<decimal128_t>();
    test_zero_with_exponent<decimal_fast32_t>();
    test_zero_with_exponent<decimal_fast64_t>();
    test_zero_with_exponent<decimal_fast128_t>();
    test_in_range<boost::int128::uint128_t>(decimal128_t{0, 60}, boost::int128::uint128_t{0U});

    test_rounded_limits();
    test_truncation_to_limits();
    test_int128_class();
    test_narrow_types();
    test_builtin_int128();

    return boost::report_errors();
}
