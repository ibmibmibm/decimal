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

// An x87 number with the largest exponent, or an unnormal one with the integer bit clear, is an infinity or a NaN
constexpr auto x87_not_finite(const int biased, const std::uint64_t significand) noexcept -> bool
{
    return biased == wide_max_biased_exponent || (biased != 0 && significand < (UINT64_C(1) << 63U));
}

// An x87 or binary128 number as significand * 2^exponent, or an infinity or a NaN. x87 keeps the leading bit of the
// significand, and binary128 implies it for a normal number.
struct wide_fields
{
    int128::uint128_t significand;
    int exponent;
    bool negative;
    bool infinite;
    bool nan;
};

template <typename Float>
constexpr auto decode_wide(const int128::uint128_t& bits) noexcept -> wide_fields
{
    constexpr bool x87 {std::numeric_limits<Float>::digits == 64};
    constexpr int fraction_bits {x87 ? 63 : 112};
    const int128::uint128_t leading {x87 ? int128::uint128_t {UINT64_C(1) << 63U} : int128::uint128_t {1U} << 112U};
    const auto biased {wide_biased_exponent<Float>(bits)};
    const auto fraction {x87 ? int128::uint128_t {bits.low & (leading.low - 1U)} : bits & (leading - 1U)};
    const auto significand {x87 ? int128::uint128_t {bits.low} : biased == 0 ? fraction : fraction | leading};
    const bool not_finite {x87 ? x87_not_finite(biased, bits.low) : biased == wide_max_biased_exponent};
    const bool infinite {biased == wide_max_biased_exponent && fraction == 0U && (!x87 || bits.low != 0U)};
    const auto exponent {(biased == 0 ? 1 : biased) - wide_exponent_bias - fraction_bits};
    return {significand, exponent, wide_sign<Float>(bits), infinite, not_finite && !infinite};
}

// binary64 has 52 fraction bits and an 11-bit biased exponent with bias 1023, whose largest value marks a NaN or an
// infinity. decode_binary64 gives its fields as decode_wide does for the wider formats.
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int binary64_fraction_bits {52};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int binary64_exponent_bias {1023};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t binary64_max_biased_exponent {0x7FF};

constexpr auto decode_binary64(const std::uint64_t bits) noexcept -> wide_fields
{
    const auto biased {(bits >> binary64_fraction_bits) & binary64_max_biased_exponent};
    const auto fraction {bits & ((UINT64_C(1) << binary64_fraction_bits) - 1U)};
    const auto significand {biased == 0U ? fraction : fraction | (UINT64_C(1) << binary64_fraction_bits)};
    const auto unbiased {(biased == 0U ? 1 : static_cast<int>(biased)) - binary64_exponent_bias};
    const auto exponent {unbiased - binary64_fraction_bits};
    const bool special {biased == binary64_max_biased_exponent};
    return {int128::uint128_t {significand}, exponent, (bits >> 63U) != 0U, special && fraction == 0U,
            special && fraction != 0U};
}

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_BINARY_BITS_HPP
