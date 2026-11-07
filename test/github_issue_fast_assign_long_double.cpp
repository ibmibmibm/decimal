// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Without long double support, a fast type can not take a long double, also not in operator=.
// The macro gives that configuration here, as on a platform with an unknown long double.

#define BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <type_traits>

using namespace boost::decimal;

static_assert(!std::is_assignable<decimal_fast32_t&, long double>::value, "decimal_fast32_t");
static_assert(!std::is_assignable<decimal_fast64_t&, long double>::value, "decimal_fast64_t");
static_assert(!std::is_assignable<decimal_fast128_t&, long double>::value, "decimal_fast128_t");
static_assert(std::is_assignable<decimal_fast64_t&, double>::value, "double");

int main()
{
    decimal_fast64_t value {};
    value = 1.5;
    BOOST_TEST_EQ(value, (decimal_fast64_t{15, -1}));

    return boost::report_errors();
}
