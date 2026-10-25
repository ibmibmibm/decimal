// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/decimal/cfloat.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

int main()
{
    BOOST_TEST_EQ(BOOST_DECIMAL_DEC32_MIN, (std::numeric_limits<decimal32_t>::min)());
    BOOST_TEST_EQ(BOOST_DECIMAL_DEC64_MIN, (std::numeric_limits<decimal64_t>::min)());
    BOOST_TEST_EQ(BOOST_DECIMAL_DEC128_MIN, (std::numeric_limits<decimal128_t>::min)());

    return boost::report_errors();
}
