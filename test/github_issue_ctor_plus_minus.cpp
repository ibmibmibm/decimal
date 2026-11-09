// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The string constructors reject a minus sign after a plus sign, as strtod and operator>> do

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <stdexcept>
#include <string>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename T>
void test(const T seven)
{
    BOOST_TEST_EQ(T{"+7"}, seven);
    BOOST_TEST_EQ(T{std::string("+7")}, seven);
    BOOST_TEST_EQ(T{"-7"}, -seven);

    #if !defined(BOOST_DECIMAL_DISABLE_EXCEPTIONS) || !defined(BOOST_DECIMAL_FAST_MATH)
    for (const char* str : {"+-7", "+-inf", "+"})
    {
        #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
        BOOST_TEST_THROWS(T{str}, std::runtime_error);
        BOOST_TEST_THROWS(T{std::string(str)}, std::runtime_error);
        #else
        BOOST_TEST(isnan(T{str}));
        BOOST_TEST(isnan(T{std::string(str)}));
        #endif
    }
    #endif
}

int main()
{
    test(7_DF);
    test(7_DD);
    test(7_DL);
    test(7_DFF);
    test(7_DDF);
    test(7_DLF);

    return boost::report_errors();
}
