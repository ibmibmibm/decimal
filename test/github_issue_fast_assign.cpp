// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The fast types assign from an integer and from a binary float, as the IEEE types do

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename T>
constexpr auto assign_int(const int v) noexcept -> T
{
    T x {};
    x = v;
    return x;
}

template <typename T>
void test()
{
    static_assert(assign_int<T>(-3) == T{-3}, "constexpr assignment");

    T x {};
    BOOST_TEST_EQ(&(x = 5), &x);
    BOOST_TEST_EQ(x, T{5});

    x = -5LL;
    BOOST_TEST_EQ(x, T{-5});

    x = 7U;
    BOOST_TEST_EQ(x, T{7});

    x = true;
    BOOST_TEST_EQ(x, T{1});

    x = false;
    BOOST_TEST_EQ(x, T{0});

    BOOST_TEST_EQ(&(x = 1.5F), &x);
    BOOST_TEST_EQ(x, (T{15, -1}));

    x = -0.25;
    BOOST_TEST_EQ(x, (T{-25, -2}));

    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    x = 0.5L;
    BOOST_TEST_EQ(x, (T{5, -1}));
    #endif

    #ifndef BOOST_DECIMAL_FAST_MATH
    x = std::numeric_limits<double>::infinity();
    BOOST_TEST(isinf(x));

    x = std::numeric_limits<double>::quiet_NaN();
    BOOST_TEST(isnan(x));
    #endif
}

int main()
{
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    decimal_fast32_t f32 {};
    f32 = 123456789;
    BOOST_TEST_EQ(f32, 1.234568e8_DFF);

    decimal_fast64_t f64 {};
    f64 = 12345678901234567890ULL;
    BOOST_TEST_EQ(f64, 1.234567890123457e19_DDF);

    return boost::report_errors();
}
