// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// An alignment character after the closing brace of a replacement field is text,
// and a brace is not a fill character

#define FMT_HEADER_ONLY
#include <boost/decimal.hpp>
#include <boost/decimal/fmt_format.hpp>
#include <boost/core/lightweight_test.hpp>
#include <string>

using namespace boost::decimal;
using namespace boost::decimal::literals;

#ifdef BOOST_DECIMAL_HAS_FORMAT_SUPPORT

template <typename T>
void test_std()
{
    const T a {1.5_DD};
    const T b {2.5_DD};

    BOOST_TEST_EQ(std::format("{}<{}", a, b), "1.5<2.5");
    BOOST_TEST_EQ(std::format("{}>{}", a, b), "1.5>2.5");
    BOOST_TEST_EQ(std::format("{}^{}", a, b), "1.5^2.5");
    BOOST_TEST_EQ(std::format("{:}<{:}", a, b), "1.5<2.5");
    BOOST_TEST_EQ(std::format("{0}>{1}<", a, b), "1.5>2.5<");
    BOOST_TEST_EQ(std::format("{:.2f}<{:5}>", a, b), "1.50<  2.5>");
    BOOST_TEST(std::format(L"{}<{}", a, b) == L"1.5<2.5");

    // The fill of a field stays as it was
    BOOST_TEST_EQ(std::format("{:*<5}<{:<5}|", a, b), "1.5**<2.5  |");

    BOOST_TEST_THROWS(static_cast<void>(std::vformat("{:{<5}", std::make_format_args(a))), std::format_error);
    BOOST_TEST_THROWS(static_cast<void>(std::vformat("{:}>5}", std::make_format_args(a))), std::format_error);
}

#endif

#ifdef BOOST_DECIMAL_HAS_FMTLIB_SUPPORT

template <typename T>
void test_fmt()
{
    const T a {1.5_DD};
    const T b {2.5_DD};

    BOOST_TEST_EQ(fmt::format("{}<{}", a, b), "1.5<2.5");
    BOOST_TEST_EQ(fmt::format("{}>{}", a, b), "1.5>2.5");
    BOOST_TEST_EQ(fmt::format("{}^{}", a, b), "1.5^2.5");
    BOOST_TEST_EQ(fmt::format("{:}<{:}", a, b), "1.5<2.5");
    BOOST_TEST_EQ(fmt::format("{0}>{1}<", a, b), "1.5>2.5<");
    BOOST_TEST_EQ(fmt::format("{:*<5}<{:<5}|", a, b), "1.5**<2.5  |");

    BOOST_TEST_THROWS(static_cast<void>(fmt::format(fmt::runtime("{:{<5}"), a)), fmt::format_error);
    BOOST_TEST_THROWS(static_cast<void>(fmt::format(fmt::runtime("{:}>5}"), a)), fmt::format_error);
}

#endif

int main()
{
    #ifdef BOOST_DECIMAL_HAS_FORMAT_SUPPORT
    test_std<decimal32_t>();
    test_std<decimal64_t>();
    test_std<decimal128_t>();
    test_std<decimal_fast32_t>();
    test_std<decimal_fast64_t>();
    test_std<decimal_fast128_t>();
    #endif

    #ifdef BOOST_DECIMAL_HAS_FMTLIB_SUPPORT
    test_fmt<decimal32_t>();
    test_fmt<decimal64_t>();
    test_fmt<decimal128_t>();
    test_fmt<decimal_fast32_t>();
    test_fmt<decimal_fast64_t>();
    test_fmt<decimal_fast128_t>();
    #endif

    return boost::report_errors();
}
