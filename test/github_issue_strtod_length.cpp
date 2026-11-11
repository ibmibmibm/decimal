// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// strtod and wcstod read only the characters that can belong to the number, so that a loop with
// endptr over a long text does not copy the rest of the text for each number

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdlib>
#include <iostream>
#include <locale>
#include <new>
#include <stdexcept>
#include <string>

using namespace boost::decimal;
using namespace boost::decimal::literals;

static std::size_t allocations {};

// Both array new forms use malloc, so that delete[] can free either, and a failure ends the test
void* operator new[](std::size_t size)
{
    void* p {std::malloc(size)};
    if (p == nullptr)
    {
        std::abort();
    }
    return p;
}

void* operator new[](std::size_t size, const std::nothrow_t&) noexcept
{
    ++allocations;
    return std::malloc(size);
}

void operator delete[](void* p) noexcept
{
    std::free(p);
}

void operator delete[](void* p, std::size_t) noexcept
{
    std::free(p);
}

template <typename T>
void test_text()
{
    std::string text;
    for (int i {}; i < 1000; ++i)
    {
        text += "1.5 -2.5e1 ";
    }
    const std::wstring wtext (text.begin(), text.end());

    allocations = 0U;
    T sum {};
    const char* p {text.c_str()};
    for (char* end {};; p = end)
    {
        const auto value {strtod<T>(p, &end)};
        if (end == p)
        {
            break;
        }
        sum += value;
    }
    BOOST_TEST_EQ(sum, T{-23500});
    BOOST_TEST_EQ(p - text.c_str(), static_cast<std::ptrdiff_t>(text.size() - 1U));

    sum = T{};
    const wchar_t* wp {wtext.c_str()};
    for (wchar_t* end {};; wp = end)
    {
        const auto value {wcstod<T>(wp, &end)};
        if (end == wp)
        {
            break;
        }
        sum += value;
    }
    BOOST_TEST_EQ(sum, T{-23500});
    BOOST_TEST_EQ(wp - wtext.c_str(), static_cast<std::ptrdiff_t>(wtext.size() - 1U));
    BOOST_TEST_EQ(allocations, 0U);
}

// A number longer than the stack buffer still reads
template <typename T>
void test_long(const char* separator)
{
    std::string text {"1"};
    for (int i {}; i < 600; ++i)
    {
        text += separator;
    }
    text += "e-1800 2";
    const std::wstring wtext (text.begin(), text.end());

    char* end {};
    BOOST_TEST_EQ(strtod<T>(text.c_str(), &end), T{1});
    BOOST_TEST_EQ(end - text.c_str(), static_cast<std::ptrdiff_t>(text.size() - 2U));

    wchar_t* wend {};
    BOOST_TEST_EQ(wcstod<T>(wtext.c_str(), &wend), T{1});
    BOOST_TEST_EQ(wend - wtext.c_str(), static_cast<std::ptrdiff_t>(wtext.size() - 2U));
}

template <typename T>
void test()
{
    test_text<T>();
    test_long<T>("000");

    #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
    try
    {
        std::locale::global(std::locale("de_DE.UTF-8"));
        test_long<T>(".000");
        BOOST_TEST_EQ(strtod<T>("1,5", nullptr), 1.5_DF);
    }
    catch (const std::runtime_error&)
    {
        std::cerr << "de_DE.UTF-8 is not available" << std::endl;
    }
    std::locale::global(std::locale::classic());
    #endif

    BOOST_TEST_EQ(strtod<T>("2.5", nullptr), 2.5_DF);
}

// A locale with ',' as the point and '.' between groups of three digits
class point_comma : public std::numpunct<char>
{
private:
    auto do_decimal_point() const -> char override { return ','; }
    auto do_thousands_sep() const -> char override { return '.'; }
    auto do_grouping() const -> std::string override { return "\3"; }
};

// A '.' that is not between two digits is not the point of that locale, and ends the number
void test_separator_as_point()
{
    std::locale::global(std::locale(std::locale::classic(), new point_comma));

    char* end {};
    const char* str {"1.e2"};
    BOOST_TEST_EQ(strtod64(str, &end), decimal64_t{1});
    BOOST_TEST_EQ(end - str, 1);

    str = "1.234,5.6";
    BOOST_TEST_EQ(strtod64(str, &end), (decimal64_t{12345, -1}));
    BOOST_TEST_EQ(end - str, 7);

    str = "-.5";
    static_cast<void>(strtod64(str, &end));
    BOOST_TEST_EQ(end - str, 0);

    // A wide character that is not a char is not a digit
    wchar_t* wide_end {};
    const wchar_t* wide {L"1.234.\u0135"};
    BOOST_TEST_EQ(wcstod64(wide, &wide_end), decimal64_t{1234});
    BOOST_TEST_EQ(wide_end - wide, 5);

    std::locale::global(std::locale::classic());
}

int main()
{
    test_separator_as_point();

    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    return boost::report_errors();
}
