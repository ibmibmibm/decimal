// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>
#include <type_traits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The subnormal values are multiples of 10^etiny, so the next value is one such unit away
void test_subnormal()
{
    BOOST_TEST_EQ(nextafter("-6.93e-97"_DF, "1"_DF), "-6.9299e-97"_DF);
    BOOST_TEST_EQ(nextafter("5.6e-98"_DF, "1"_DF), "5.601e-98"_DF);
    BOOST_TEST_EQ(nextafter("1e-390"_DD, "1"_DD), "1.00000001e-390"_DD);
    BOOST_TEST_EQ(nextafter("1e-390"_DD, "0"_DD), "9.9999999e-391"_DD);
    BOOST_TEST_EQ(nextafter("1e-6170"_DL, "1"_DL), "1.000001e-6170"_DL);
    BOOST_TEST_EQ(nextafter("-1e-6170"_DL, "-1"_DL), "-1.000001e-6170"_DL);
}

// The largest subnormal value, with p - 1 digits, and the smallest normal value are neighbors
void test_normal_boundary()
{
    BOOST_TEST_EQ(nextafter("9.99999999999999e-384"_DD, "0"_DD), "9.99999999999998e-384"_DD);
    BOOST_TEST_EQ(nextafter("9.99999999999999e-384"_DD, "1"_DD), "1e-383"_DD);
    BOOST_TEST_EQ(nextafter("1e-383"_DD, "0"_DD), "9.99999999999999e-384"_DD);
    BOOST_TEST_EQ(nextafter("9.99999e-96"_DF, "1"_DF), "1e-95"_DF);
    BOOST_TEST_EQ(nextafter("1e-95"_DF, "-1"_DF), "9.99999e-96"_DF);
}

// The smallest subnormal value steps to a zero with the sign of x
template <typename T>
void test_to_zero()
{
    const auto tiny {std::numeric_limits<T>::denorm_min()};

    const auto pos {nextafter(tiny, T(0))};
    BOOST_TEST_EQ(pos, T(0));
    BOOST_TEST(!signbit(pos));

    const auto neg {nextafter(-tiny, T(0))};
    BOOST_TEST_EQ(neg, T(0));
    BOOST_TEST(signbit(neg));
}

void test_normal()
{
    BOOST_TEST_EQ(nextafter("1.5"_DF, "2"_DF), "1.500001"_DF);
    BOOST_TEST_EQ(nextafter("1"_DD, "0"_DD), "0.9999999999999999"_DD);
    BOOST_TEST_EQ(nextafter("1.5"_DFF, "2"_DFF), "1.500001"_DFF);
    BOOST_TEST_EQ(nextafter("1"_DLF, "0"_DLF), "0.9999999999999999999999999999999999"_DLF);

    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isinf(nextafter((std::numeric_limits<decimal64_t>::max)(), std::numeric_limits<decimal64_t>::infinity())));
    #endif
}

// For two types, the step is one unit of the common type
void test_mixed_types()
{
    BOOST_TEST_EQ(nextafter("1"_DF, "2"_DD), "1.000000000000001"_DD);
    BOOST_TEST_EQ(nextafter("1"_DD, "0"_DL), "0.9999999999999999999999999999999999"_DL);
    BOOST_TEST_EQ(nextafter("1e-101"_DF, "0"_DD), "9.999999999999999e-102"_DD);
}

// Behavior that all six types share
template <typename T>
void test_shared()
{
    const auto tiny {std::numeric_limits<T>::denorm_min()};
    const auto min_normal {(std::numeric_limits<T>::min)()};

    // A zero steps to the zero of the direction, and two equal cohorts give the direction
    BOOST_TEST(signbit(nextafter(T(0), -T(0))));
    BOOST_TEST(!signbit(nextafter(-T(0), T(0))));
    BOOST_TEST_EQ(nextafter(T(10, -1), T(1)), T(1));
    BOOST_TEST_EQ(nextafter(T(0), -T(1)), -tiny);

    // The lower boundary steps out and back to the same value
    BOOST_TEST_EQ(nextafter(nextafter(tiny, T(1)), T(0)), tiny);
    BOOST_TEST_EQ(nextafter(nextafter(min_normal, T(0)), T(1)), min_normal);
    BOOST_TEST_LT(nextafter(min_normal, T(0)), min_normal);
    BOOST_TEST_GT(nextafter(-min_normal, T(0)), -min_normal);

    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto inf {std::numeric_limits<T>::infinity()};
    BOOST_TEST_EQ(nextafter(inf, inf), inf);
    BOOST_TEST_EQ(nextafter(-inf, -inf), -inf);
    BOOST_TEST_EQ(nextafter(inf, T(0)), (std::numeric_limits<T>::max)());
    BOOST_TEST_EQ(nextafter(-inf, T(0)), -(std::numeric_limits<T>::max)());
    #endif
}

// The result type and the step do not depend on the order of the arguments
void test_mixed_order()
{
    static_assert(std::is_same<decltype(nextafter("1"_DF, "2"_DD)), decimal64_t>::value, "Wrong type");
    static_assert(std::is_same<decltype(nextafter("2"_DD, "1"_DF)), decimal64_t>::value, "Wrong type");

    BOOST_TEST_EQ(nextafter("2"_DD, "1"_DF), "1.999999999999999"_DD);
    BOOST_TEST_EQ(nextafter("1"_DL, "2"_DD), "1.000000000000000000000000000000001"_DL);
    BOOST_TEST_EQ(nextafter("1"_DF, "2"_DFF), static_cast<decltype(nextafter("1"_DF, "2"_DFF))>("1.000001"_DF));

    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST_EQ(nextafter(std::numeric_limits<decimal32_t>::infinity(), std::numeric_limits<decimal64_t>::infinity()),
                  std::numeric_limits<decimal64_t>::infinity());
    #endif
}

int main()
{
    test_subnormal();
    test_normal_boundary();
    test_to_zero<decimal32_t>();
    test_to_zero<decimal64_t>();
    test_to_zero<decimal128_t>();
    test_normal();
    test_mixed_types();
    test_mixed_order();

    test_shared<decimal32_t>();
    test_shared<decimal64_t>();
    test_shared<decimal128_t>();
    test_shared<decimal_fast32_t>();
    test_shared<decimal_fast64_t>();
    test_shared<decimal_fast128_t>();

    return boost::report_errors();
}
