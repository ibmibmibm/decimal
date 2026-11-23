// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The fields of a record that rounding_grid writes and bound_mutants reads.
// See README.md in this directory.

#ifndef BOOST_DECIMAL_TOOLS_RECORD_FIELD_HPP
#define BOOST_DECIMAL_TOOLS_RECORD_FIELD_HPP

#include <limits>
#include <string>

// A special value as the records write it: nan, inf or -inf. false for each other field, also for a second sign.
template <typename T>
auto parse_special(const std::string& s, T& x) -> bool
{
    if (s != "nan" && s != "inf" && s != "-inf")
    {
        return false;
    }
    const T inf {std::numeric_limits<T>::infinity()};
    x = s == "nan" ? std::numeric_limits<T>::quiet_NaN() : s == "inf" ? inf : -inf;
    return true;
}

#endif // BOOST_DECIMAL_TOOLS_RECORD_FIELD_HPP
