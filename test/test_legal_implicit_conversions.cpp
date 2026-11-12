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

// A type with more digits holds each value of a type with fewer digits, and an IEEE type
// holds each value of the fast type of its size
template <typename From, typename To>
void test_wider()
{
    const From values[] {From {0, 0}, From {0, 5}, From {0, -5}, From {1, 0}, From {100, -2}, From {1234567, -3},
                         From {9999999, 20}, From {1, -30}, std::numeric_limits<From>::max(),
                         std::numeric_limits<From>::min(), std::numeric_limits<From>::denorm_min(),
                         std::numeric_limits<From>::epsilon()};

    for (const From value : values)
    {
        for (const From from_val : {value, -value})
        {
            const To to_val {from_val};

            BOOST_TEST_EQ(from_val, to_val);
            BOOST_TEST_EQ(signbit(from_val), signbit(to_val));
            BOOST_TEST_EQ(static_cast<From>(to_val), from_val);
        }
    }

    for (const From from_val : {std::numeric_limits<From>::infinity(), -std::numeric_limits<From>::infinity()})
    {
        const To to_val {from_val};

        BOOST_TEST(isinf(to_val));
        BOOST_TEST_EQ(signbit(from_val), signbit(to_val));
    }

    BOOST_TEST(isnan(static_cast<To>(std::numeric_limits<From>::quiet_NaN())));
    BOOST_TEST(isnan(static_cast<To>(std::numeric_limits<From>::signaling_NaN())));
}

int main()
{
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

    test_wider<boost::decimal::decimal32_t, boost::decimal::decimal64_t>();
    test_wider<boost::decimal::decimal32_t, boost::decimal::decimal128_t>();
    test_wider<boost::decimal::decimal64_t, boost::decimal::decimal128_t>();
    test_wider<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast64_t>();
    test_wider<boost::decimal::decimal_fast32_t, boost::decimal::decimal_fast128_t>();
    test_wider<boost::decimal::decimal_fast64_t, boost::decimal::decimal_fast128_t>();
    test_wider<boost::decimal::decimal_fast32_t, boost::decimal::decimal64_t>();
    test_wider<boost::decimal::decimal_fast32_t, boost::decimal::decimal128_t>();
    test_wider<boost::decimal::decimal_fast64_t, boost::decimal::decimal128_t>();
    test_wider<boost::decimal::decimal32_t, boost::decimal::decimal_fast64_t>();
    test_wider<boost::decimal::decimal32_t, boost::decimal::decimal_fast128_t>();
    test_wider<boost::decimal::decimal64_t, boost::decimal::decimal_fast128_t>();
    test_wider<boost::decimal::decimal_fast32_t, boost::decimal::decimal32_t>();
    test_wider<boost::decimal::decimal_fast64_t, boost::decimal::decimal64_t>();
    test_wider<boost::decimal::decimal_fast128_t, boost::decimal::decimal128_t>();

    return boost::report_errors();
}
