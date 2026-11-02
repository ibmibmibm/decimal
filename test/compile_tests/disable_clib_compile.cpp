// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_DECIMAL_DISABLE_CLIB
#include <boost/decimal.hpp>

// A constant expression can not hold an active assertion that fails
constexpr auto no_assertions() noexcept -> int
{
    BOOST_DECIMAL_ASSERT(false);
    BOOST_DECIMAL_ASSERT_MSG(false, "Not active");
    return 0;
}

static_assert(no_assertions() == 0, "Assertions must be off without the C library");

int main()
{
    return 0;
}
