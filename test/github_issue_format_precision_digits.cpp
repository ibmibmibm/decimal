// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// A precision of a format specification needs digits after the point

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
    const T x {1.5_DD};

    #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
    for (const auto* text : {"{:.}", "{:.f}", "{:.e}", "{:5.}", "{:.L}", "{:.g}"})
    {
        BOOST_TEST_THROWS(static_cast<void>(std::vformat(text, std::make_format_args(x))), std::format_error);
    }
    #endif

    BOOST_TEST_EQ(std::format("{:.0f}", x), "2");
    BOOST_TEST_EQ(std::format("{:.1f}", x), "1.5");
    BOOST_TEST_EQ(std::format("{:5.00f}|", x), "    2|");
}

#endif

#ifdef BOOST_DECIMAL_HAS_FMTLIB_SUPPORT

template <typename T>
void test_fmt()
{
    const T x {1.5_DD};

    #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
    for (const auto* text : {"{:.}", "{:.f}", "{:.e}", "{:5.}", "{:.L}", "{:.g}"})
    {
        BOOST_TEST_THROWS(static_cast<void>(fmt::format(fmt::runtime(text), x)), fmt::format_error);
    }
    #endif

    BOOST_TEST_EQ(fmt::format("{:.0f}", x), "2");
    BOOST_TEST_EQ(fmt::format("{:.1f}", x), "1.5");
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
