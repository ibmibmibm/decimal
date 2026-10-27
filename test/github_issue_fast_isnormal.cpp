// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The fast types have no subnormal values: each finite non-zero value is normal,
// and gives the same result as the IEEE type of the same format
template <typename Fast, typename Ieee>
void test(const Fast x, const Ieee y)
{
    BOOST_TEST(isnormal(x));
    BOOST_TEST(isnormal(-x));
    BOOST_TEST_EQ(fpclassify(x), FP_NORMAL);

    BOOST_TEST_EQ(static_cast<Ieee>(nextafter(x, Fast{1})), nextafter(y, Ieee{1}));
    if (x != std::numeric_limits<Fast>::min())
    {
        BOOST_TEST_EQ(static_cast<Ieee>(nextafter(x, Fast{0})), nextafter(y, Ieee{0}));
    }
    BOOST_TEST_EQ(expm1(x), x);
    BOOST_TEST_EQ(static_cast<Ieee>(log10(x)), log10(y));
    BOOST_TEST_EQ(static_cast<Ieee>(tgamma(x)), tgamma(y));
    BOOST_TEST_EQ(pow(x, Fast{2}), Fast{0});
}

int main()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    BOOST_TEST(isnormal(std::numeric_limits<decimal_fast32_t>::min()));
    BOOST_TEST(isnormal(std::numeric_limits<decimal_fast64_t>::min()));
    BOOST_TEST(isnormal(std::numeric_limits<decimal_fast128_t>::min()));

    test(1e-95_dff, 1e-95_df);
    test(1.234567e-90_dff, 1.234567e-90_df);
    test(9.999999e-90_dff, 9.999999e-90_df);
    test(1e-383_ddf, 1e-383_dd);
    test(1.234567e-374_ddf, 1.234567e-374_dd);
    test(9.999999999999999e-369_ddf, 9.999999999999999e-369_dd);
    test(1e-6143_dlf, 1e-6143_dl);
    test(1.234567e-6134_dlf, 1.234567e-6134_dl);
    test(9.999999999999999999999999999999999e-6111_dlf, 9.999999999999999999999999999999999e-6111_dl);
    #endif

    return boost::report_errors();
}
