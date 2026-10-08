// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The same checks with BOOST_DECIMAL_DETAIL_INT128_NO_BUILTIN_INT128. The int128 classes then have no conversion
// to a built-in 128-bit integer, but the decimal types still convert to one.

#define BOOST_DECIMAL_DETAIL_INT128_NO_BUILTIN_INT128

#include "github_issue_to_integral.cpp"
