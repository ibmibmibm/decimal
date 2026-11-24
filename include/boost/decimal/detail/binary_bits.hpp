// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_BINARY_BITS_HPP
#define BOOST_DECIMAL_DETAIL_BINARY_BITS_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/int128.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <cstring>
#include <limits>
#endif

// The bits of the binary formats wider than double, for the conversions and for nexttoward

namespace boost {
namespace decimal {
namespace detail {

// The x87 and binary128 formats have a 15-bit biased exponent with this bias, and its largest value marks a NaN or an
// infinity
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int wide_exponent_bias {16383};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int wide_max_biased_exponent {0x7FFF};

// Reads inside the 10 bytes of an x87 number, as one past them waits on the store. GCC can widen a read of bytes 8 and
// 9, and GCC 16 narrows one of bytes 2 to 9 to them. The callers take a reference, as Clang copies a long double.
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto wide_bits(const Float& val) noexcept -> int128::uint128_t
{
    static_assert(std::numeric_limits<Float>::digits != 64 || BOOST_DECIMAL_ENDIAN_LITTLE_BYTE,
                  "The x87 bits are read in the little-endian layout of x86");
    BOOST_DECIMAL_IF_CONSTEXPR (std::numeric_limits<Float>::digits == 64)
    {
        // A constant expression reads the bits with bit_cast, as memcpy is not one
        #ifdef BOOST_DECIMAL_HAS_CONSTEXPR_BITCAST
        if (BOOST_DECIMAL_IS_CONSTANT_EVALUATED(val))
        {
            const auto bits {bit_cast<x87_bits<Float>>(val)};
            return {bits.high, bits.low};
        }
        #endif

        std::uint64_t low {};
        std::memcpy(&low, &val, sizeof(low));
        #if defined(__GNUC__) && !defined(__clang__)
        std::uint64_t high {};
        std::memcpy(&high, reinterpret_cast<const unsigned char*>(&val) + 2, sizeof(high));
        return {high >> 48U, low};
        #else
        std::uint16_t high {};
        std::memcpy(&high, reinterpret_cast<const unsigned char*>(&val) + sizeof(low), sizeof(high));
        return {high, low};
        #endif
    }
    return bit_cast<int128::uint128_t>(val);
}

// The biased exponent and the sign of an x87 or binary128 number, from its bits
template <typename Float>
constexpr auto wide_biased_exponent(const int128::uint128_t& bits) noexcept -> int
{
    return static_cast<int>((std::numeric_limits<Float>::digits == 64 ? bits.high : bits.high >> 48U) & static_cast<std::uint64_t>(wide_max_biased_exponent));
}

template <typename Float>
constexpr auto wide_sign(const int128::uint128_t& bits) noexcept -> bool
{
    return ((bits.high >> (std::numeric_limits<Float>::digits == 64 ? 15U : 63U)) & 1U) != 0U;
}

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_BINARY_BITS_HPP
