// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A compound assignment with any other decimal type computes in the promoted type
// and converts the result to the left type

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename Lhs, typename Rhs>
void test_pair()
{
    const Lhs a {75, -1};
    const Rhs b {25, -1};

    Lhs r {a};
    BOOST_TEST_EQ(&(r += b), &r);
    BOOST_TEST_EQ(r, Lhs{10});

    r = a;
    BOOST_TEST_EQ(&(r -= b), &r);
    BOOST_TEST_EQ(r, Lhs{5});

    r = a;
    BOOST_TEST_EQ(&(r *= b), &r);
    BOOST_TEST_EQ(r, (Lhs{1875, -2}));

    r = a;
    BOOST_TEST_EQ(&(r /= b), &r);
    BOOST_TEST_EQ(r, Lhs{3});
}

template <typename Lhs>
void test_lhs()
{
    test_pair<Lhs, decimal32_t>();
    test_pair<Lhs, decimal64_t>();
    test_pair<Lhs, decimal128_t>();
    test_pair<Lhs, decimal_fast32_t>();
    test_pair<Lhs, decimal_fast64_t>();
    test_pair<Lhs, decimal_fast128_t>();
}

// The sum rounds to the digits of the left type, after the wider sum
void test_narrowing()
{
    auto d32 {1_DF};
    d32 += 1.23456789_DD;
    BOOST_TEST_EQ(d32, 2.234568_DF);

    auto f32 {1_DFF};
    f32 += 1.23456789_DL;
    BOOST_TEST_EQ(f32, 2.234568_DFF);

    auto d64 {1_DD};
    d64 -= 1.2345678901234567890_DLF;
    BOOST_TEST_EQ(d64, -0.2345678901234568_DD);

    auto f64 {1_DDF};
    f64 *= 1.2345678901234567890_DL;
    BOOST_TEST_EQ(f64, 1.234567890123457_DDF);

    auto d128 {1_DL};
    d128 /= 3_DLF;
    BOOST_TEST_EQ(d128, 1_DL / 3_DL);

    auto f128 {1_DLF};
    f128 /= 3_DD;
    BOOST_TEST_EQ(f128, 1_DLF / 3_DLF);
}

int main()
{
    test_lhs<decimal32_t>();
    test_lhs<decimal64_t>();
    test_lhs<decimal128_t>();
    test_lhs<decimal_fast32_t>();
    test_lhs<decimal_fast64_t>();
    test_lhs<decimal_fast128_t>();
    test_narrowing();

    return boost::report_errors();
}
