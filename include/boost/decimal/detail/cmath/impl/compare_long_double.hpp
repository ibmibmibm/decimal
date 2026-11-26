// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_COMPARE_LONG_DOUBLE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_COMPARE_LONG_DOUBLE_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/binary_bits.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/fast_float/big_uint.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/type_traits.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#endif

namespace boost {
namespace decimal {
namespace detail {

#ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

// A finite long double as significand * 2^exponent, or else an infinity or a NaN, also for binary64
using long_double_split = wide_fields;

// A template, as Clang rejects the memcpy of the x87 bits in a constexpr function that is not one
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto split_long_double(const Float& val) noexcept -> long_double_split
{
    #if BOOST_DECIMAL_LDBL_BITS == 64
    return decode_binary64(bit_cast<std::uint64_t>(val));
    #else
    return decode_wide<Float>(wide_bits(val));
    #endif
}

// A function call, because MSVC warns about a constant condition in C++14
template <typename DecimalType>
constexpr auto has_wide_exponent() noexcept -> bool { return decimal_val_v<DecimalType> >= 128; }

// The order test leaves this many bits for the error of the two bit length estimates below, which is under 5.2 bits
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int bit_length_margin {8};

// 5^27 < 2^63, so for |exp| <= 27 both values fit in 128 bits when sig and the significand of y have one word each
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int max_one_word_exp {27};

// Each type below 128 bits has |exp| < bias_v + precision_v, so it never needs the largest compare
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int max_narrow_exp {550};

// The largest |exp| of a decimal value within bit_length_margin of the long double range, as sig is below 2^113
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int max_long_double_exp {
    (std::numeric_limits<long double>::digits - std::numeric_limits<long double>::min_exponent + 113 + bit_length_margin) *
    512 / 1701 + 1};

// Words for an exact compare with |exp| <= max_exp. The integer without a shift is sig or the significand of y, below
// 2^113, times at most 5^|exp|, and 1189 / 512 > log2(5). The shifted one is within 2^(bit_length_margin + 5) of it.
constexpr auto exact_compare_words(const int max_exp) noexcept -> std::size_t
{
    return static_cast<std::size_t>((113 + 1 + max_exp * 1189 / 512 + bit_length_margin + 5) / 64 + 1);
}

// The compares take the sizes of the big integer that the conversions to binary use: 24 words, or 200 words
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::size_t long_double_compare_words {
    exact_compare_words(max_long_double_exp) <= 24U ? 24U : 200U};

static_assert(exact_compare_words(max_narrow_exp) <= 24U &&
              exact_compare_words(max_long_double_exp) <= long_double_compare_words,
              "The exact compares need more words");

// The sign of y - sig * 10^exp for a value of DecimalType: sig below 2^113 and |exp| below 6300. The magnitudes compare
// first by their bit lengths. Values within a few binades of each other compare as integers, in 128 bits if exp is small.
template <typename DecimalType, typename SigType>
BOOST_DECIMAL_FORCE_INLINE constexpr auto compare_long_double(const long_double_split& y, const SigType sig, const int exp, const bool is_neg) noexcept -> int
{
    static_assert(has_wide_exponent<DecimalType>() || bias_v<DecimalType> + precision_v<DecimalType> <= max_narrow_exp,
                  "A type below 128 bits needs the largest compare");

    if (y.significand == 0U || sig == 0U)
    {
        if (y.significand == 0U && sig == 0U)
        {
            return 0;
        }
        return y.significand != 0U ? (y.negative ? -1 : 1) : (is_neg ? 1 : -1);
    }

    if (y.negative != is_neg)
    {
        return y.negative ? -1 : 1;
    }

    const int128::uint128_t wide_sig {sig};
    const int y_bits {y.exponent + 128 - countl_zero(y.significand)};

    // 1701 / 512 is log2(10) to within 3.4e-4, so x_bits is off by less than 4.2 bits for |exp| below 6300.
    // y_bits is off by less than 1 bit.
    const int x_bits {128 - countl_zero(wide_sig) + exp * 1701 / 512};

    int order {};
    if (y_bits > x_bits + bit_length_margin)
    {
        order = 1;
    }
    else if (y_bits < x_bits - bit_length_margin)
    {
        order = -1;
    }
    else if (wide_sig.high == 0U && y.significand.high == 0U && exp >= -max_one_word_exp && exp <= max_one_word_exp)
    {
        // 5^|exp| has at most 63 bits, so both values are an integer below 2^127 times a power of 2
        const auto shift {static_cast<unsigned>(exp < 0 ? -exp : exp)};
        const auto pow5 {static_cast<std::uint64_t>(pow10(int128::uint128_t {shift}) >> shift)};
        auto x {exp < 0 ? wide_sig : wide_sig * pow5};
        auto y_sig {exp < 0 ? y.significand * pow5 : y.significand};
        const int x_zeros {countl_zero(x)};
        const int y_zeros {countl_zero(y_sig)};
        x <<= x_zeros;
        y_sig <<= y_zeros;
        if (y.exponent - y_zeros != exp - x_zeros)
        {
            order = y.exponent - y_zeros > exp - x_zeros ? 1 : -1;
        }
        else
        {
            order = y_sig > x ? 1 : (y_sig < x ? -1 : 0);
        }
    }
    else if (!has_wide_exponent<DecimalType>() || (exp >= -max_narrow_exp && exp <= max_narrow_exp))
    {
        order = -fast_float::compare_exact(exp, wide_sig, y.significand, y.exponent);
    }
    else
    {
        // |exp| <= max_long_double_exp here, as else the bit lengths differ by more than bit_length_margin
        order = -fast_float::compare_exact<long_double_compare_words>(exp, wide_sig, y.significand, y.exponent);
    }

    return is_neg ? -order : order;
}

#endif // BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_COMPARE_LONG_DOUBLE_HPP
