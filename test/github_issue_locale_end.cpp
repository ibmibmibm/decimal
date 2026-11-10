// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// strtod and operator>> remove the digit separators of the locale between the digits of the integer
// part, give the end in the input, and operator>> fails on a grouping that num_put does not give

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <iterator>
#include <locale>
#include <sstream>
#include <string>

using namespace boost::decimal;

#ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS

template <typename T>
void check(const char* str, const T value, const int used, const bool grouping_ok = true)
{
    char* end {};
    BOOST_TEST_EQ(strtod<T>(str, &end), value);
    BOOST_TEST_EQ(end - str, used);

    const std::string s {str};
    const std::wstring ws (s.begin(), s.end());
    wchar_t* wend {};
    BOOST_TEST_EQ(wcstod<T>(ws.c_str(), &wend), value);
    BOOST_TEST_EQ(wend - ws.c_str(), used);

    std::istringstream is {s};
    T read {};
    is >> read;
    BOOST_TEST_EQ(read, value);
    BOOST_TEST_EQ(is.fail(), !grouping_ok);
    const std::string rest {std::istreambuf_iterator<char>(is.rdbuf()), std::istreambuf_iterator<char>()};
    BOOST_TEST_EQ(rest, s.substr(static_cast<std::size_t>(used)));
}

template <typename T>
void check_none(const char* str)
{
    char* end {};
    const auto value {strtod<T>(str, &end)};
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isnan(value));
    #else
    static_cast<void>(value);
    #endif
    BOOST_TEST_EQ(end, str);

    const std::string s {str};
    std::istringstream is {s};
    T read {};
    is >> read;
    BOOST_TEST(is.fail());
    const std::string rest {std::istreambuf_iterator<char>(is.rdbuf()), std::istreambuf_iterator<char>()};
    BOOST_TEST_EQ(rest, s);
}

template <typename T>
void test_en()
{
    check("1,234,567xyz", T{1234567}, 9);
    check("12,345.5,6", T{123455, -1}, 8);
    check("12,34.5,6", T{12345, -1}, 7, false);
    check("-1,000e3 rest", T{-1, 6}, 8);
    check("1,234", T{1234}, 5);
    check("1,234,", T{1234}, 5);
    check("  1,234x", T{1234}, 7);
    check("1e3,5", T{1000}, 3);
    check("1,,2", T{1}, 1);
    check("1,2,3", T{123}, 5, false);
    check("1234,567", T{1234567}, 8, false);
    check("1,234,56", T{123456}, 8, false);
    check_none<T>(",5");
}

template <typename T>
void test_de()
{
    check("1.234,5x", T{12345, -1}, 7);
    check("1.234.567,25.5", T{123456725, -2}, 12);
    check("1.23,5", T{1235, -1}, 6, false);
    check("1.,5", T{1}, 1);
}

template <typename T>
void test()
{
    try
    {
        std::locale::global(std::locale("en_US.UTF-8"));
        test_en<T>();
    }
    catch (const std::runtime_error&)
    {
        std::cerr << "en_US.UTF-8 is not installed" << std::endl;
    }

    try
    {
        std::locale::global(std::locale("de_DE.UTF-8"));
        test_de<T>();
    }
    catch (const std::runtime_error&)
    {
        std::cerr << "de_DE.UTF-8 is not installed" << std::endl;
    }

    std::locale::global(std::locale::classic());
}

int main()
{
    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    return boost::report_errors();
}

#else

int main()
{
    return 0;
}

#endif
