// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>

using namespace boost::decimal;

// Defined in github_issue_fesetround_2.cpp, another translation unit
auto other_unit_rounding_mode() noexcept -> rounding_mode;
auto other_unit_third() noexcept -> decimal32_t;

int main()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    using namespace boost::decimal::literals;

    boost::decimal::fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST(other_unit_rounding_mode() == rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(other_unit_third(), "0.3333334"_DF);

    boost::decimal::fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST(other_unit_rounding_mode() == rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(other_unit_third(), "0.3333333"_DF);

    boost::decimal::fesetround(rounding_mode::fe_dec_default);
    #endif

    return boost::report_errors();
}
