// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>

using namespace boost::decimal;

auto other_unit_rounding_mode() noexcept -> rounding_mode
{
    return boost::decimal::fegetround();
}

auto other_unit_third() noexcept -> decimal32_t
{
    return decimal32_t{1} / decimal32_t{3};
}
