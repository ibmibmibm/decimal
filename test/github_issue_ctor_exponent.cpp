// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <climits>
#include <cstdint>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename T>
void test_wide_exponents()
{
    const auto inf {std::numeric_limits<T>::infinity()};
    const T zero {0U, 0};

    BOOST_TEST_EQ((T{1U, INT64_C(1) << 31}), inf);
    BOOST_TEST_EQ((T{1U, INT64_C(1) << 32}), inf);
    BOOST_TEST_EQ((T{1U, INT64_MAX}), inf);
    BOOST_TEST_EQ((T{1U, UINT32_MAX}), inf);
    BOOST_TEST_EQ((T{1U, UINT64_MAX}), inf);
    BOOST_TEST_EQ((T{1U, INT_MAX}), inf);

    BOOST_TEST_EQ((T{1U, -(INT64_C(1) << 31) - 1}), zero);
    BOOST_TEST_EQ((T{1U, -(INT64_C(1) << 32)}), zero);
    BOOST_TEST_EQ((T{1U, INT64_MIN}), zero);
    BOOST_TEST_EQ((T{1U, INT_MIN}), zero);

    BOOST_TEST_EQ((T{-1, INT64_C(1) << 32}), -inf);
    BOOST_TEST_EQ((T{1U, boost::int128::int128_t{INT64_C(1) << 40}}), inf);
    BOOST_TEST_EQ((T{0U, INT64_C(1) << 32}), zero);
    BOOST_TEST_EQ((T{5U, UINT64_C(2)}), (T{500U, 0}));
    BOOST_TEST_EQ((T{5U, INT64_C(-2)}), (T{500U, -4}));
}

template <typename T>
void test_narrow_exponents()
{
    BOOST_TEST_EQ((T{12U, static_cast<unsigned char>(3)}), (T{12000U, 0}));
    BOOST_TEST_EQ((T{12U, static_cast<unsigned short>(3)}), (T{12000U, 0}));
    BOOST_TEST_EQ((T{12U, static_cast<signed char>(-3)}), (T{12U, -3}));
    BOOST_TEST_EQ((T{12U, static_cast<short>(-3)}), (T{12U, -3}));
    BOOST_TEST_EQ((T{12U, true}), (T{120U, 0}));
    BOOST_TEST_EQ((T{12U, static_cast<signed char>(-90)}), (T{12U, -90}));
    BOOST_TEST_EQ((T{12U, static_cast<unsigned short>(65535)}), std::numeric_limits<T>::infinity());
    BOOST_TEST_EQ((T{12U, static_cast<short>(-32768)}), (T{0U, 0}));
}

int main()
{
    test_wide_exponents<decimal32_t>();
    test_wide_exponents<decimal64_t>();
    test_wide_exponents<decimal128_t>();
    test_wide_exponents<decimal_fast32_t>();
    test_wide_exponents<decimal_fast64_t>();
    test_wide_exponents<decimal_fast128_t>();

    test_narrow_exponents<decimal32_t>();
    test_narrow_exponents<decimal64_t>();
    test_narrow_exponents<decimal128_t>();
    test_narrow_exponents<decimal_fast32_t>();
    test_narrow_exponents<decimal_fast64_t>();
    test_narrow_exponents<decimal_fast128_t>();

    BOOST_TEST_EQ((decimal64_t{125U, static_cast<signed char>(-2)}), 1.25_DD);

    return boost::report_errors();
}
