// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <initializer_list>
#include <limits>

using namespace boost::decimal;

// The step out from the largest finite value is exact, so it is infinity in each rounding mode
template <typename T>
void test()
{
    const T max_value {(std::numeric_limits<T>::max)()};
    const T inf {std::numeric_limits<T>::infinity()};

    BOOST_TEST_EQ(nextafter(max_value, inf), inf);
    BOOST_TEST_EQ(nextafter(-max_value, -inf), -inf);
    BOOST_TEST_LT(nextafter(max_value, T{0}), max_value);
    BOOST_TEST_GT(nextafter(-max_value, T{0}), -max_value);

    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    BOOST_TEST_EQ(nexttoward(max_value, std::numeric_limits<long double>::infinity()), inf);
    BOOST_TEST_EQ(nexttoward(-max_value, -std::numeric_limits<long double>::infinity()), -inf);
    #endif
}

// A finite direction past the largest finite value also steps to infinity, also with fast math
void test_finite_direction()
{
    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    const decimal32_t max_value {(std::numeric_limits<decimal32_t>::max)()};
    const auto inf_bits {to_bits(std::numeric_limits<decimal32_t>::infinity())};

    BOOST_TEST_EQ(to_bits(nexttoward(max_value, 1e100L)), inf_bits);
    BOOST_TEST_EQ(to_bits(-nexttoward(-max_value, -1e100L)), inf_bits);
    #endif
}

// Runs the check in each rounding mode
template <typename Check>
void in_all_modes(const Check check)
{
    check();

    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    for (const auto mode : {rounding_mode::fe_dec_downward, rounding_mode::fe_dec_upward,
                            rounding_mode::fe_dec_toward_zero, rounding_mode::fe_dec_to_nearest_from_zero})
    {
        fesetround(mode);
        check();
    }
    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    in_all_modes(test<decimal32_t>);
    in_all_modes(test<decimal64_t>);
    in_all_modes(test<decimal128_t>);
    in_all_modes(test<decimal_fast32_t>);
    in_all_modes(test<decimal_fast64_t>);
    in_all_modes(test<decimal_fast128_t>);
    #endif

    in_all_modes(test_finite_direction);

    return boost::report_errors();
}
