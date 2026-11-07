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

template <typename T>
void test_all_modes()
{
    test<T>();

    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    for (const auto mode : {rounding_mode::fe_dec_downward, rounding_mode::fe_dec_upward,
                            rounding_mode::fe_dec_toward_zero, rounding_mode::fe_dec_to_nearest_from_zero})
    {
        fesetround(mode);
        test<T>();
    }
    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

int main()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    test_all_modes<decimal32_t>();
    test_all_modes<decimal64_t>();
    test_all_modes<decimal128_t>();
    test_all_modes<decimal_fast32_t>();
    test_all_modes<decimal_fast64_t>();
    test_all_modes<decimal_fast128_t>();
    #endif

    return boost::report_errors();
}
