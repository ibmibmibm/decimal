// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// frexp of an infinity gives the same infinity (C F.10.3.4). The exponent 0 is a choice of this library.
template <typename T>
void test()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto x : {std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity()})
    {
        int e {1};
        const auto m {frexp(x, &e)};
        BOOST_TEST(isinf(m));
        BOOST_TEST_EQ(signbit(m), signbit(x));
        BOOST_TEST_EQ(e, 0);
    }
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
