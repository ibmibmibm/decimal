// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_BIT_LAYOUTS_HPP
#define BOOST_DECIMAL_DETAIL_BIT_LAYOUTS_HPP

#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cfloat>
#endif

#if LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
#  define BOOST_DECIMAL_LDBL_BITS 80
#elif LDBL_MANT_DIG == 113 && LDBL_MAX_EXP == 16384
#  define BOOST_DECIMAL_LDBL_BITS 128
#elif LDBL_MANT_DIG == 53 && LDBL_MAX_EXP == 1024
#  define BOOST_DECIMAL_LDBL_BITS 64
#else // Unsupported long double representation
#  define BOOST_DECIMAL_LDBL_BITS 0
#  define BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
#endif

#endif // BOOST_DECIMAL_DETAIL_BIT_LAYOUTS_HPP
