// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// The fast types flush a sub-normal result to zero
template <typename T>
void test_fast()
{
    // These members were deprecated in C++23
    #if ((!defined(_MSC_VER) && (__cplusplus <= 202002L)) || (defined(_MSC_VER) && (_MSVC_LANG <= 202002L)))
    BOOST_TEST(std::numeric_limits<T>::has_denorm == std::denorm_absent);
    BOOST_TEST(!std::numeric_limits<T>::has_denorm_loss);
    #endif
    BOOST_TEST_EQ(std::numeric_limits<T>::denorm_min(), (std::numeric_limits<T>::min)());
    BOOST_TEST_EQ((std::numeric_limits<T>::min)() / T{10}, T{0});
}

int main()
{
    test_fast<decimal_fast32_t>();
    test_fast<decimal_fast64_t>();
    test_fast<decimal_fast128_t>();

    return boost::report_errors();
}
