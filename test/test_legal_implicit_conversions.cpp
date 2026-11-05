// Copyright 2025 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

template <typename From, typename To>
void test_implicit()
{
    const From from_val {2, 1};
    const To to_val = from_val;

    BOOST_TEST_EQ(from_val, to_val);
}

template <typename From, typename To>
void test_explicit()
{
    const From from_val {2, 1};
    const To to_val = static_cast<To>(from_val);

    BOOST_TEST_EQ(from_val, to_val);
}

// A zero of either sign at any stored exponent is false. The smallest value, an infinity and both NaNs are true.
template <typename T>
void test_bool()
{
    BOOST_TEST(!static_cast<bool>(T{0, 0}));
    BOOST_TEST(!static_cast<bool>(T{0U, -5, true}));
    BOOST_TEST(!static_cast<bool>(T{0, 20}));
    BOOST_TEST(static_cast<bool>(std::numeric_limits<T>::denorm_min()));
    BOOST_TEST(static_cast<bool>(-std::numeric_limits<T>::denorm_min()));
    BOOST_TEST(static_cast<bool>(-std::numeric_limits<T>::infinity()));
    BOOST_TEST(static_cast<bool>(std::numeric_limits<T>::quiet_NaN()));
    BOOST_TEST(static_cast<bool>(std::numeric_limits<T>::signaling_NaN()));
}

int main()
{
    test_bool<boost::decimal::decimal32_t>();
    test_bool<boost::decimal::decimal_fast32_t>();
    test_bool<boost::decimal::decimal64_t>();
    test_bool<boost::decimal::decimal_fast64_t>();
    test_bool<boost::decimal::decimal128_t>();
    test_bool<boost::decimal::decimal_fast128_t>();

    test_explicit<boost::decimal::decimal32_t, boost::decimal::decimal32_t>();
    test_explicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal32_t>();
    test_explicit<boost::decimal::decimal64_t, boost::decimal::decimal32_t>();
    test_explicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal32_t>();
    test_explicit<boost::decimal::decimal128_t, boost::decimal::decimal32_t>();
    test_explicit<boost::decimal::decimal_fast128_t, boost::decimal::decimal32_t>();

    test_explicit<boost::decimal::decimal32_t, boost::decimal::decimal_fast32_t>();
    test_explicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast32_t>();
    test_explicit<boost::decimal::decimal64_t, boost::decimal::decimal_fast32_t>();
    test_explicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal_fast32_t>();
    test_explicit<boost::decimal::decimal128_t, boost::decimal::decimal_fast32_t>();
    test_explicit<boost::decimal::decimal_fast128_t, boost::decimal::decimal_fast32_t>();

    test_explicit<boost::decimal::decimal32_t, boost::decimal::decimal_fast64_t>();
    test_explicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast64_t>();
    test_explicit<boost::decimal::decimal64_t, boost::decimal::decimal_fast64_t>();
    test_explicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal_fast64_t>();
    test_explicit<boost::decimal::decimal128_t, boost::decimal::decimal_fast64_t>();
    test_explicit<boost::decimal::decimal_fast128_t, boost::decimal::decimal_fast64_t>();

    test_explicit<boost::decimal::decimal32_t, boost::decimal::decimal128_t>();
    test_explicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal128_t>();
    test_explicit<boost::decimal::decimal64_t, boost::decimal::decimal128_t>();
    test_explicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal128_t>();
    test_explicit<boost::decimal::decimal128_t, boost::decimal::decimal128_t>();
    test_explicit<boost::decimal::decimal_fast128_t, boost::decimal::decimal128_t>();

    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal32_t>();
    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal_fast32_t>();
    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal64_t>();
    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal_fast64_t>();
    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal128_t>();
    test_implicit<boost::decimal::decimal32_t, boost::decimal::decimal_fast128_t>();

    test_implicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast32_t>();
    test_implicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal64_t>();
    test_implicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast64_t>();
    test_implicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal128_t>();
    test_implicit<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast128_t>();

    test_implicit<boost::decimal::decimal64_t, boost::decimal::decimal64_t>();
    test_implicit<boost::decimal::decimal64_t, boost::decimal::decimal_fast64_t>();
    test_implicit<boost::decimal::decimal64_t, boost::decimal::decimal128_t>();
    test_implicit<boost::decimal::decimal64_t, boost::decimal::decimal_fast128_t>();

    test_implicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal_fast64_t>();
    test_implicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal128_t>();
    test_implicit<boost::decimal::decimal_fast64_t, boost::decimal::decimal_fast128_t>();

    test_implicit<boost::decimal::decimal128_t, boost::decimal::decimal128_t>();
    test_implicit<boost::decimal::decimal128_t, boost::decimal::decimal_fast128_t>();

    test_implicit<boost::decimal::decimal_fast128_t, boost::decimal::decimal_fast128_t>();

    return boost::report_errors();
}
