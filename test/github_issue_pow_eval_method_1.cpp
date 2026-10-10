// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The pow checks with BOOST_DECIMAL_DEC_EVAL_METHOD 1. pow uses no wider type, so it rounds once and the results do
// not change. Each cmath function returns the type of its argument.

#define BOOST_DECIMAL_DEC_EVAL_METHOD 1

#include "github_issue_pow.cpp"
#include <type_traits>

static_assert(std::is_same<decltype(pow(decimal32_t {}, decimal32_t {})), decimal32_t>::value, "pow");
static_assert(std::is_same<decltype(exp(decimal32_t {})), decimal32_t>::value, "exp");
static_assert(std::is_same<decltype(hypot(decimal32_t {}, decimal32_t {})), decimal32_t>::value, "hypot");
static_assert(std::is_same<decltype(lgamma(decimal_fast32_t {})), decimal_fast32_t>::value, "lgamma");
static_assert(std::is_same<decltype(beta(decimal64_t {}, decimal64_t {})), decimal64_t>::value, "beta");
static_assert(std::is_same<decltype(legendre(2U, decimal64_t {})), decimal64_t>::value, "legendre");
static_assert(std::is_same<decltype(atan2(decimal_fast64_t {}, decimal_fast64_t {})), decimal_fast64_t>::value,
              "atan2");
