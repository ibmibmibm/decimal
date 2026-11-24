// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// u256.hpp is the only include, so it must include the AVX2 declarations itself
#include <boost/decimal/detail/u256.hpp>

#if !defined(__AVX2__) || defined(BOOST_DECIMAL_NO_CONSTEVAL_DETECTION)
#  error "The u256 AVX2 path is not active"
#endif

int main()
{
    return 0;
}
