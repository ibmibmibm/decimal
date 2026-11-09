// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// from_chars, strtod and nan read nan(n-char-sequence) as C does

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstring>

using namespace boost::decimal;

template <typename T>
void check(const char* str, const int used, const unsigned payload, const bool signaling)
{
    T value;
    const auto r {from_chars(str, str + std::strlen(str), value)};
    BOOST_TEST(r.ec == std::errc());
    BOOST_TEST_EQ(r.ptr - str, used);
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isnan(value));
    BOOST_TEST_EQ(issignaling(value), signaling);
    BOOST_TEST(read_payload(value) == payload);
    #else
    static_cast<void>(payload);
    static_cast<void>(signaling);
    #endif

    char* end {};
    static_cast<void>(strtod<T>(str, &end));
    BOOST_TEST_EQ(end - str, used);
}

template <typename T>
void test()
{
    check<T>("nan(5)", 6, 5U, false);
    check<T>("nan(snan7)", 10, 7U, true);
    check<T>("nan5x", 4, 5U, false);
    check<T>("snan5x", 5, 5U, true);

    // Letters, digits and underscores, up to the closing parenthesis
    check<T>("nan()", 5, 0U, false);
    check<T>("nan(_)", 6, 0U, false);
    check<T>("nan(a_b)", 8, 0U, false);
    check<T>("nan(12abc)", 10, 0U, false);
    check<T>("nan(0x10)", 9, 0U, false);
    check<T>("nan(snanx)", 10, 0U, false);
    check<T>("nan(index)", 10, 0U, false);
    check<T>("nan(1234567890123456789012345678901234567890)", 45, 0U, false);

    // Without the closing parenthesis only "nan" is read
    check<T>("nan(5", 3, 0U, false);
    check<T>("nan(5 )", 3, 0U, false);
    check<T>("nan(snan", 3, 0U, false);
    check<T>("nan(ind", 3, 0U, false);

    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(read_payload(nan<T>("12")) == 12U);
    BOOST_TEST(read_payload(nan<T>("12abc")) == 0U);
    BOOST_TEST(read_payload(snan<T>("5 ")) == 0U);
    BOOST_TEST(issignaling(snan<T>("5 ")));
    #endif
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
