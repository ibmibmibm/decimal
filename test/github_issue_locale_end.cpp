// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// strtod and operator>> remove the digit separators of the locale, and give the end in the input

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
void check(const char* str, const T value, const int used)
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
    const std::string rest {std::istreambuf_iterator<char>(is.rdbuf()), std::istreambuf_iterator<char>()};
    BOOST_TEST_EQ(rest, s.substr(static_cast<std::size_t>(used)));
}

template <typename T>
void test_en()
{
    check("1,234,567xyz", T{1234567}, 9);
    check("12,34.5,6", T{12345, -1}, 7);
    check("-1,000e3 rest", T{-1, 6}, 8);
    check("1,234", T{1234}, 5);
    check("1,234,", T{1234}, 5);
}

template <typename T>
void test_de()
{
    check("1.234,5x", T{12345, -1}, 7);
    check("1.234.567,25.5", T{123456725, -2}, 12);
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
