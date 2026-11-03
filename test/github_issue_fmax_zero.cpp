// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;

// fmax and fmin order -0 below +0, as IEEE 754-2019 9.6 does. C F.10.9.2 calls this desirable, not required
template <typename T1, typename T2>
void test()
{
    const T1 pz1 {0};
    const T2 pz2 {0};
    const T1 nz1 {-pz1};
    const T2 nz2 {-pz2};

    BOOST_TEST(!signbit(fmax(pz1, nz2)));
    BOOST_TEST(!signbit(fmax(nz1, pz2)));
    BOOST_TEST(signbit(fmax(nz1, nz2)));
    BOOST_TEST(!signbit(fmax(pz1, pz2)));

    BOOST_TEST(signbit(fmin(pz1, nz2)));
    BOOST_TEST(signbit(fmin(nz1, pz2)));
    BOOST_TEST(signbit(fmin(nz1, nz2)));
    BOOST_TEST(!signbit(fmin(pz1, pz2)));

    BOOST_TEST_EQ(fmax(T1{-1}, T2{2}), T2{2});
    BOOST_TEST_EQ(fmin(T1{-1}, T2{2}), T1{-1});
}

int main()
{
    test<decimal32_t, decimal32_t>();
    test<decimal64_t, decimal64_t>();
    test<decimal128_t, decimal128_t>();
    test<decimal_fast32_t, decimal_fast32_t>();
    test<decimal_fast64_t, decimal_fast64_t>();
    test<decimal_fast128_t, decimal_fast128_t>();
    test<decimal32_t, decimal128_t>();
    test<decimal_fast64_t, decimal32_t>();

    return boost::report_errors();
}
