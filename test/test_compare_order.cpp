// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The order of values at the limits of the comparison paths: signs, cohorts, the smallest and
// largest values, and exponent gaps around the size of the products

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

template <typename T>
void test_pair(const T& low, const T& high)
{
    BOOST_TEST(low < high);
    BOOST_TEST(low <= high);
    BOOST_TEST(high > low);
    BOOST_TEST(high >= low);
    BOOST_TEST(low != high);
    BOOST_TEST(!(high < low));
    BOOST_TEST(!(high <= low));
    BOOST_TEST(!(low > high));
    BOOST_TEST(!(low >= high));
    BOOST_TEST(!(low == high));

    // Negation reverses the order
    BOOST_TEST(-high < -low);
    BOOST_TEST(-high <= -low);
    BOOST_TEST(!(-low < -high));

    // Different signs
    if (high > T{0})
    {
        BOOST_TEST(-high < high);
        BOOST_TEST(!(high < -high));
    }
}

template <typename T>
void test_same(const T& a, const T& b)
{
    BOOST_TEST(a == b);
    BOOST_TEST(a <= b);
    BOOST_TEST(a >= b);
    BOOST_TEST(!(a != b));
    BOOST_TEST(!(a < b));
    BOOST_TEST(!(a > b));
    BOOST_TEST(-a == -b);
    BOOST_TEST(!(-a < -b));
    BOOST_TEST(!(-b < -a));
}

template <typename T>
void test_limits()
{
    using limits = std::numeric_limits<T>;

    // The smallest normal value has the smallest biased exponent in the fast types
    test_pair(T{0}, (limits::min)());
    test_pair((limits::min)(), T{1});
    test_pair(T{1}, (limits::max)());
    test_pair(limits::lowest(), (limits::min)());
    BOOST_TEST(-(limits::min)() < (limits::min)());
    BOOST_TEST(!((limits::min)() < -(limits::min)()));
    BOOST_TEST(-(limits::min)() <= (limits::min)());
    BOOST_TEST(!((limits::min)() <= -(limits::min)()));
    BOOST_TEST(limits::lowest() < (limits::max)());

    // Zeros are equal whatever their signs and exponents are
    test_same(T{0}, -T{0});
    test_same(T(0, 5), T(0, -5));
    test_same(T(0, 5), -T(0, -5));
    BOOST_TEST(-T{0} == 0);
    BOOST_TEST(T{5} == 5);
    BOOST_TEST(T{-5} == -5);
    BOOST_TEST(T{5} != -5);
    BOOST_TEST(T{-5} != 5);
    const T zero_cohort {0, 5};
    BOOST_TEST(-T{1} < zero_cohort);
    BOOST_TEST(zero_cohort < T{1});
    BOOST_TEST(-T{1} < -zero_cohort);

    // Infinities
    test_pair((limits::max)(), limits::infinity());
    BOOST_TEST(-limits::infinity() < limits::lowest());
    BOOST_TEST(-limits::infinity() < limits::infinity());
    BOOST_TEST(!(limits::infinity() < limits::infinity()));
    BOOST_TEST(limits::infinity() <= limits::infinity());
    BOOST_TEST(!(limits::infinity() < -limits::infinity()));

    // NaN is not ordered
    BOOST_TEST(!(limits::quiet_NaN() < T{1}));
    BOOST_TEST(!(T{1} < limits::quiet_NaN()));
    BOOST_TEST(!(limits::quiet_NaN() <= T{1}));
    BOOST_TEST(!(T{1} <= limits::quiet_NaN()));
    BOOST_TEST(!(limits::quiet_NaN() == limits::quiet_NaN()));
}

// Exponent gaps from 0 to 40, on the two sides of equality
template <typename T>
void test_gaps()
{
    using sig_type = typename T::significand_type;
    const int digits {std::numeric_limits<T>::digits10};

    for (int gap {0}; gap <= 40; ++gap)
    {
        const T large {sig_type{1U}, gap};
        const T small_value {sig_type{9U}, gap - 1};
        test_pair(small_value, large);
        test_pair(T{sig_type{1U}, 0}, T{sig_type{2U}, gap});

        if (gap < digits)
        {
            sig_type scaled {1U};
            for (int i {0}; i < gap; ++i)
            {
                scaled = static_cast<sig_type>(scaled * 10U);
            }

            test_same(T{sig_type{1U}, gap}, T{scaled, 0});
            BOOST_TEST((T{sig_type{1U}, gap} != -T{scaled, 0}));
            BOOST_TEST((-T{sig_type{1U}, gap} != T{scaled, 0}));
            test_pair(T{scaled, 0}, T{static_cast<sig_type>(scaled + 1U), 0});
            test_pair(T{sig_type{1U}, gap}, T{static_cast<sig_type>(scaled + 1U), 0});
            if (gap > 0)
            {
                test_pair(T{static_cast<sig_type>(scaled - 1U), 0}, T{sig_type{1U}, gap});
            }
        }
    }
}

int main()
{
    test_limits<decimal32_t>();
    test_limits<decimal64_t>();
    test_limits<decimal128_t>();
    test_limits<decimal_fast32_t>();
    test_limits<decimal_fast64_t>();
    test_limits<decimal_fast128_t>();

    test_gaps<decimal32_t>();
    test_gaps<decimal64_t>();
    test_gaps<decimal128_t>();
    test_gaps<decimal_fast32_t>();
    test_gaps<decimal_fast64_t>();
    test_gaps<decimal_fast128_t>();

    return boost::report_errors();
}
