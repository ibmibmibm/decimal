// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The check for the AVX2 compile test: the flag is valid for this compiler and target, and u256.hpp uses AVX2
#include <boost/decimal/detail/config.hpp>

#if !defined(__AVX2__) || defined(BOOST_DECIMAL_NO_CONSTEVAL_DETECTION)
#  error "The u256 AVX2 path is not active"
#endif

int main()
{
    return 0;
}
