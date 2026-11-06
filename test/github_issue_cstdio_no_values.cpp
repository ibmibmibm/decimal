// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cerrno>
#include <cstring>

using namespace boost::decimal;
using namespace boost::decimal::literals;

int main()
{
    char buffer[64] {};

    // The format alone is a valid call, as for std::snprintf
    BOOST_TEST_EQ(boost::decimal::snprintf(buffer, sizeof(buffer), "100%% text"), 9);
    BOOST_TEST_CSTR_EQ(buffer, "100% text");

    errno = 0;
    BOOST_TEST_EQ(boost::decimal::snprintf(buffer, sizeof(buffer), "%Hf"), -1);
    BOOST_TEST_EQ(errno, EINVAL);

    // The common type of the values does not change
    BOOST_TEST_EQ(boost::decimal::snprintf(buffer, sizeof(buffer), "%.30DDf", 1.2345678901234567890123456789_DL), 32);
    BOOST_TEST_CSTR_EQ(buffer, "1.234567890123456789012345678900");

    BOOST_TEST_EQ(boost::decimal::snprintf(buffer, sizeof(buffer), "%.1Df %.1DDf", 1.5_DF, 2.5_DLF), 7);
    BOOST_TEST_CSTR_EQ(buffer, "1.5 2.5");

    return boost::report_errors();
}
