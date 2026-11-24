// Copyright 2026 Pavel P
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP
#define BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/add_impl.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/binary_bits.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/dragonbox/dragonbox.hpp>
#include <boost/decimal/detail/fast_float/compute_float64.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/remove_trailing_zeros.hpp>
#include <boost/decimal/detail/utilities.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace binary_to_decimal {

// Selects the constructor of a decimal type that has the code of the conversion
struct by_reference {};

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4324) // Structure was padded due to alignment specifier
#endif

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpadded"
#endif

// The value is significand * 10^exponent, plus a nonzero fraction of the last digit when inexact
template <typename UInt>
struct scaled_value
{
    UInt significand;
    int exponent;
    // A word, as GCC returned a bool through the stack: it cleared the padding bytewise, and the word load stalled
    std::uint32_t inexact;
};

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

// x * 10^k for Dragonbox's k, 18 or 19 digits for a normal double, exact as Dragonbox proves its cache precision.
// It is forced inline, as GCC keeps it out of line in from_double_narrow.
BOOST_DECIMAL_FORCE_INLINE constexpr auto dragonbox_scaled(const std::uint64_t two_fc, const int binary_exponent) noexcept
    -> scaled_value<std::uint64_t>
{
    using binary64 = dragonbox::binary64_format;

    const auto minus_k {dragonbox::floor_log10_pow2(binary_exponent) - binary64::kappa};
    const auto cache {binary64::get_cache(minus_k)};
    const auto beta {binary_exponent + dragonbox::floor_log2_pow10(-minus_k)};

    const auto high {dragonbox::umul128(two_fc, cache.high)};
    const auto low {dragonbox::umul128(two_fc, cache.low)};
    const auto middle {high.low + low.high};
    const auto top {high.high + static_cast<std::uint64_t>(middle < low.high)};

    const auto integer_part {(top << beta) | (middle >> (64 - beta))};
    const auto is_integer {((middle << beta) | (low.low >> (64 - beta))) == 0U};
    return {integer_part, minus_k, static_cast<std::uint32_t>(!is_integer)};
}

// The same for a float from the float table, with 8 to 10 digits for a normal float
BOOST_DECIMAL_FORCE_INLINE constexpr auto dragonbox_scaled(const std::uint32_t two_fc, const int binary_exponent) noexcept
    -> scaled_value<std::uint64_t>
{
    using binary32 = dragonbox::binary32_format;

    const auto minus_k {dragonbox::floor_log10_pow2(binary_exponent) - binary32::kappa};
    const auto beta {binary_exponent + dragonbox::floor_log2_pow10(-minus_k)};
    const auto product {dragonbox::umul128(two_fc, binary32::get_cache(minus_k))};
    return {static_cast<std::uint64_t>(product >> (64 - beta)), minus_k,
            static_cast<std::uint32_t>(static_cast<std::uint32_t>(product.low >> (32 - beta)) != 0U)};
}

// The binary type whose table converts Float: float for 24 significand bits, else double
template <typename Float>
using narrow_binary = std::conditional_t<std::numeric_limits<Float>::digits == 24, float, double>;

template <typename Binary>
using binary_uint = std::conditional_t<sizeof(Binary) == 4, std::uint32_t, std::uint64_t>;

// dragonbox_scaled of the bits of a float or a double
template <typename Binary>
BOOST_DECIMAL_FORCE_INLINE constexpr auto scaled_bits(const binary_uint<Binary> bits) noexcept
    -> scaled_value<std::uint64_t>
{
    constexpr int fraction_bits {std::numeric_limits<Binary>::digits - 1};
    constexpr int offset {std::numeric_limits<Binary>::max_exponent - 2 + std::numeric_limits<Binary>::digits};
    const auto exponent_bits {static_cast<int>(bits >> fraction_bits) & (2 * std::numeric_limits<Binary>::max_exponent - 1)};
    const auto fraction {bits & ((binary_uint<Binary> {1U} << fraction_bits) - 1U)};
    const auto significand {exponent_bits == 0 ? fraction : fraction | (binary_uint<Binary> {1U} << fraction_bits)};
    return dragonbox_scaled(significand * 2U, exponent_bits == 0 ? 1 - offset : exponent_bits - offset);
}

// x * 10^k exactly with 37 or 38 digits (log10 x is in [estimate, estimate + 1.31)), for decimal128_t and subnormals
BOOST_DECIMAL_CXX20_CONSTEXPR auto exact_scaled(const std::uint64_t significand, const int binary_exponent) noexcept
    -> scaled_value<int128::uint128_t>
{
    const auto estimate {dragonbox::floor_log10_pow2(binary_exponent + 63 - countl_zero(significand))};

    const auto k {36 - estimate};

    std::uint64_t words[24];
    fast_float::big_uint value {words, 24U, int128::uint128_t {significand}};
    bool inexact {};
    if (k >= 0)
    {
        value.mul_pow5(k);
        const auto shift {binary_exponent + k};
        if (shift >= 0)
        {
            value.shl(shift);
        }
        else
        {
            inexact = value.shr(-shift);
        }
    }
    else
    {
        value.shl(binary_exponent);
        auto pow10_left {-k};
        for (; pow10_left >= 19; pow10_left -= 19)
        {
            inexact = value.divide(UINT64_C(10000000000000000000)) != 0U || inexact;
        }
        inexact = value.divide(pow10(static_cast<std::uint64_t>(pow10_left))) != 0U || inexact;
    }
    return {value.to_uint128(), -k, inexact};
}

// Constant divisors, which compilers turn into multiplications
constexpr auto divide_by_pow10(std::uint64_t& value, int count) noexcept -> void
{
    for (; count >= 4; count -= 4)
    {
        value /= 10000U;
    }
    for (; count > 0; --count)
    {
        value /= 10U;
    }
}

constexpr auto divide_by_pow10(int128::uint128_t& value, const int count) noexcept -> void
{
    if (count > 0)
    {
        value = impl::divmod_pow10_dispatch(value, count, pow10(static_cast<int128::uint128_t>(count))).quotient;
    }
}

template <typename DecimalType, typename UInt>
constexpr auto to_decimal_type(UInt significand, int exponent, const bool inexact, const bool sign) noexcept
    -> DecimalType
{
    if (inexact)
    {
        // There is at least one digit past the precision
        const auto dropped {num_digits(significand) - precision_v<DecimalType>};
        if (exponent + dropped + precision_v<DecimalType> - 1 >= std::numeric_limits<DecimalType>::min_exponent10)
        {
            divide_by_pow10(significand, dropped - 1);
            exponent += dropped - 1 + fenv_round<DecimalType>(significand, sign, true);
            return pack_in_range<DecimalType>(significand, exponent, sign);
        }

        // A subnormal result is rounded by the constructor, from all the digits rounded to odd: they can no longer
        // read as exactly zero or half, so it rounds them as it would the exact value
        if (significand % 5U == 0U)
        {
            ++significand;
        }
        return DecimalType {significand, exponent, sign};
    }

    // Exact: IEEE 754 5.4.2 prefers exponent 0, so the coefficient keeps the zeros above it as far as precision allows
    using coefficient_type = std::conditional_t<(precision_v<DecimalType> > 19), int128::uint128_t, UInt>;
    const auto trimmed {remove_trailing_zeros(significand)};
    coefficient_type coefficient {trimmed.trimmed_number};
    exponent += static_cast<int>(trimmed.number_of_removed_zeros);
    const auto room {precision_v<DecimalType> - num_digits(coefficient)};
    if (room < 0)
    {
        return DecimalType {coefficient, exponent, sign};
    }
    if (exponent > 0)
    {
        const auto moved {exponent < room ? exponent : room};
        coefficient *= pow10(static_cast<coefficient_type>(moved));
        exponent -= moved;
    }
    return pack_in_range<DecimalType>(coefficient, exponent, sign);
}

template <typename DecimalType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_binary64(const double val) noexcept -> DecimalType
{
    const auto bits {bit_cast<std::uint64_t>(val)};
    const auto sign {(bits >> 63U) != 0U};
    const auto exponent_bits {static_cast<int>((bits >> 52U) & 0x7FFU)};
    auto significand {bits & ((UINT64_C(1) << 52U) - 1U)};

    if (exponent_bits == 0x7FF)
    {
        return DecimalType {val};
    }
    if (exponent_bits == 0 && significand == 0U)
    {
        return DecimalType {0U, 0, sign};
    }

    auto binary_exponent {-1074};
    if (exponent_bits != 0)
    {
        significand |= UINT64_C(1) << 52U;
        binary_exponent = exponent_bits - 1075;
    }

    constexpr auto enough_digits {precision_v<DecimalType> < 19 ?
                                  pow10(static_cast<std::uint64_t>(precision_v<DecimalType>)) : UINT64_MAX};

    // An integer, or a fraction with few binary places, is exact as its odd part times 2^shift at exponent 0, or
    // times 5^-shift at exponent shift: the exponents closest to 0 that IEEE 754 5.4.2 prefers
    const auto trailing_zeros {int128::detail::countr_zero(significand)};
    const auto odd {significand >> trailing_zeros};
    const auto shift {binary_exponent + trailing_zeros};
    if (shift >= 0)
    {
        if (shift < 64 && odd < (enough_digits >> shift))
        {
            return pack_in_range<DecimalType>(odd << shift, 0, sign);
        }
    }
    else if (shift >= -19)
    {
        const auto pow5 {pow10(static_cast<std::uint64_t>(-shift)) >> -shift};
        const auto product {dragonbox::umul128(odd, pow5)};
        if (product.high == 0U && product.low < enough_digits)
        {
            return pack_in_range<DecimalType>(product.low, shift, sign);
        }
    }

    // Rounding needs a digit past the precision whenever the value is inexact
    const auto scaled {dragonbox_scaled(significand * 2U, binary_exponent)};
    if (!scaled.inexact || scaled.significand >= enough_digits)
    {
        return to_decimal_type<DecimalType>(scaled.significand, scaled.exponent, scaled.inexact != 0U, sign);
    }

    const auto exact {exact_scaled(significand, binary_exponent)};
    return to_decimal_type<DecimalType>(exact.significand, exact.exponent, exact.inexact != 0U, sign);
}

// A function call, because MSVC warns about a constant condition in an if
template <bool value>
constexpr auto constant_test() noexcept -> bool { return value; }

// Shortest digits one past the precision that end in 5 are a tie, which the exact value of the binary number breaks:
// the 5 alone rounds as above it, and a 4 with the sticky digit as below it
template <typename DecimalType, typename Float, typename UInt>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_shortest(const Float val,
                                                 const dragonbox::floating_decimal<UInt>& components) noexcept
    -> DecimalType
{
    // The test for one digit past the precision wraps to false if the precision holds all the digits of UInt
    constexpr auto pow10_precision {precision_v<DecimalType> < std::numeric_limits<UInt>::digits10 ?
                                    pow10(static_cast<UInt>(precision_v<DecimalType>)) : UInt {}};
    if (components.mantissa - pow10_precision < pow10_precision * 9U && components.mantissa % 10U == 5U)
    {
        // A subnormal number has fewer digits, but Dragonbox takes its shortest digits from the same scale
        using binary = narrow_binary<Float>;
        const auto scaled {scaled_bits<binary>(bit_cast<binary_uint<binary>>(static_cast<binary>(val)))};
        const auto shortest {static_cast<std::uint64_t>(components.mantissa) *
                             pow10(static_cast<std::uint64_t>(components.exponent - scaled.exponent))};
        if (shortest != scaled.significand || scaled.inexact)
        {
            const std::uint64_t digits {scaled.significand >= shortest ? components.mantissa :
                                                                         components.mantissa - 1U};
            return to_decimal_type<DecimalType>(digits, components.exponent, true, components.sign);
        }
    }

    // A float to 7 digits has one or two digits more, and a double to 16 one, so constant divisors round them. A fast
    // type rounds a double in its constructor, which pack_in_range calls. Digits that fit are packed directly.
    constexpr int extra_digits {std::numeric_limits<Float>::max_digits10 - precision_v<DecimalType>};
    constexpr bool rounds {extra_digits == 2 || (extra_digits == 1 && !is_fast_type_v<DecimalType>)};
    if (constant_test<rounds>() && components.mantissa >= pow10_precision)
    {
        const auto two {components.mantissa >= pow10_precision * 10U};
        auto rounded {two ? components.mantissa / 10U : components.mantissa};
        const auto sticky {rounded * (two ? 10U : 1U) != components.mantissa};
        const auto exponent {components.exponent + static_cast<int>(two) +
                             fenv_round<DecimalType>(rounded, components.sign, sticky)};
        return pack_in_range<DecimalType>(rounded, exponent, components.sign);
    }
    if (constant_test<extra_digits <= 0 || rounds || is_fast_type_v<DecimalType>>())
    {
        return pack_in_range<DecimalType>(components.mantissa, components.exponent, components.sign);
    }
    return DecimalType {components.mantissa, components.exponent, components.sign};
}

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t ten_33 {UINT64_C(0x3159FC9E26C9), UINT64_C(0x5B6A4B4BF9480000)};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t ten_34 {UINT64_C(0x1ED09BEAD87C0), UINT64_C(0x378D8E6400000000)};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t ten_35 {UINT64_C(0x13426172C74D82), UINT64_C(0x2B878FE800000000)};

// 2^64 is 6 modulo 10, so the last digit needs no division of 128 bits
constexpr auto last_digit(const int128::uint128_t value) noexcept -> std::uint64_t
{
    return (value.low % 10U + 6U * (value.high % 10U)) % 10U;
}

// Digits below 10^36 split at 10^17 into two words, whose trailing zeros are cheaper to remove than those of 128 bits
constexpr auto remove_trailing_zeros_split(const int128::uint128_t value) noexcept
    -> remove_trailing_zeros_return<int128::uint128_t>
{
    const auto split {impl::divmod_pow10_uint128(value, 17)};
    const auto high {static_cast<std::uint64_t>(split.quotient)};
    const auto low {static_cast<std::uint64_t>(split.remainder)};
    if (low == 0U)
    {
        const auto trimmed {remove_trailing_zeros(high)};
        return {int128::uint128_t {trimmed.trimmed_number}, trimmed.number_of_removed_zeros + 17U};
    }
    const auto trimmed {remove_trailing_zeros(low)};
    const std::uint64_t kept_zeros {17U - trimmed.number_of_removed_zeros};
    const auto scale {pow10(kept_zeros)};
    return {dragonbox::umul128(high, scale) + trimmed.trimmed_number, trimmed.number_of_removed_zeros};
}

// Shortest digits that fit: a fast type keeps all its digits, so they move up to the precision, and the others lose
// their trailing zeros
template <typename DecimalType, typename UInt>
constexpr auto pack_fitting(const UInt digits, const int exponent, const bool sign) noexcept -> DecimalType
{
    if (is_fast_type_v<DecimalType>)
    {
        const auto room {precision_v<DecimalType> - num_digits(digits)};
        return pack_in_range<DecimalType>(int128::uint128_t {digits} * pow10(static_cast<int128::uint128_t>(room)),
                                          exponent - room, sign);
    }
    if (last_digit(int128::uint128_t {digits}) != 0U)
    {
        return pack_in_range<DecimalType>(digits, exponent, sign);
    }
    const auto trimmed {remove_trailing_zeros(digits)};
    return pack_in_range<DecimalType>(trimmed.trimmed_number,
                                      exponent + static_cast<int>(trimmed.number_of_removed_zeros), sign);
}

// The wider types take the side of the exact value from Ryu
template <typename DecimalType>
constexpr auto from_shortest(const ryu::floating_decimal_128& components) noexcept -> DecimalType
{
    // To 34 digits, binary128 has one or two digits more. If x is below them, it is above them less 1, so the digits
    // round with a sticky bit, which also breaks a tie.
    if (precision_v<DecimalType> == 34 && components.mantissa >= ten_34)
    {
        auto digits {components.direction < 0 ? components.mantissa - 1U : components.mantissa};
        auto sticky {components.direction != 0};
        auto exponent {components.exponent};
        if (components.mantissa > ten_35)
        {
            const auto split {impl::divmod10(digits)};
            digits = split.quotient;
            sticky = sticky || split.remainder != 0U;
            ++exponent;
        }
        exponent += fenv_round<DecimalType>(digits, components.sign, sticky);
        return pack_in_range<DecimalType>(digits, exponent, components.sign);
    }

    if (components.direction != 0 && last_digit(components.mantissa) == 5U &&
        num_digits(components.mantissa) == precision_v<DecimalType> + 1)
    {
        const auto digits {components.direction > 0 ? components.mantissa : components.mantissa - 1U};
        return to_decimal_type<DecimalType>(digits, components.exponent, true, components.sign);
    }
    return DecimalType {components.mantissa, components.exponent, components.sign};
}

// The interval of a double holds at most one 7-digit decimal, and of a float two. If no 7-digit neighbor of x is in it,
// the shortest digits are longer and x rounds once; if one alone is in it, it is the shortest. Else Dragonbox decides.
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_double_narrow(const Float val) noexcept -> DecimalType
{
    // A float takes the float table and a double the double table
    using binary = narrow_binary<Float>;
    using format = std::conditional_t<sizeof(binary) == 4, dragonbox::binary32_format, dragonbox::binary64_format>;
    constexpr int fraction_bits {std::numeric_limits<binary>::digits - 1};
    const auto bits {bit_cast<binary_uint<binary>>(static_cast<binary>(val))};
    const auto exponent_bits {static_cast<int>(bits >> fraction_bits) & (2 * std::numeric_limits<binary>::max_exponent - 1)};
    const auto sign {(bits >> (sizeof(binary) * 8U - 1U)) != 0U};
    if (exponent_bits == 0 || (bits & ((binary_uint<binary> {1U} << fraction_bits) - 1U)) == 0U)
    {
        return from_shortest<DecimalType>(val, floating_point_to_fd(val));
    }

    // 8 to 10 digits for a float, 18 or 19 for a double, and the interval width in the same unit
    const auto scaled {scaled_bits<binary>(bits)};
    const auto beta {exponent_bits - (std::numeric_limits<binary>::max_exponent - 2 + std::numeric_limits<binary>::digits) +
                     dragonbox::floor_log2_pow10(-scaled.exponent)};
    const auto delta {format::compute_delta(format::get_cache(scaled.exponent), beta)};

    // Each division names its constant divisor, so that the compiler multiplies instead
    constexpr int fewest_digits {sizeof(binary) == 4 ? 8 : 18};
    constexpr auto short_unit {pow10(static_cast<std::uint64_t>(fewest_digits - precision_v<DecimalType>))};
    const int form {scaled.significand >= pow10(static_cast<std::uint64_t>(fewest_digits + 1)) ? 2 :
                    scaled.significand >= pow10(static_cast<std::uint64_t>(fewest_digits)) ? 1 : 0};
    const auto unit {form == 2 ? short_unit * 100U : form == 1 ? short_unit * 10U : short_unit};
    const auto digits {form == 2 ? scaled.significand / (short_unit * 100U) :
                       form == 1 ? scaled.significand / (short_unit * 10U) : scaled.significand / short_unit};
    const auto remainder {scaled.significand - digits * unit};
    const auto exponent {scaled.exponent + fewest_digits + form - precision_v<DecimalType>};
    if (exponent + precision_v<DecimalType> - 1 < std::numeric_limits<DecimalType>::min_exponent10)
    {
        return from_shortest<DecimalType>(val, floating_point_to_fd(val));
    }

    // The fraction below the last digit is less than 1, and delta is the width rounded down
    const auto below_outside {2U * remainder > delta + 1U};
    const auto above_outside {2U * (unit - remainder - 1U) >= delta + 1U};
    const auto inside {2U * (remainder + 1U) <= delta && above_outside ? digits :
                       2U * (unit - remainder) < delta && below_outside ? digits + 1U : 0U};
    if (inside != 0U)
    {
        // A fast type keeps all its digits
        if (is_fast_type_v<DecimalType>)
        {
            return pack_in_range<DecimalType>(inside, exponent, sign);
        }
        const auto trimmed {remove_trailing_zeros(static_cast<std::uint32_t>(inside))};
        return pack_in_range<DecimalType>(trimmed.trimmed_number,
                                          exponent + static_cast<int>(trimmed.number_of_removed_zeros), sign);
    }
    if (below_outside && above_outside)
    {
        // The digits and the one after them round with the sticky rest, as to_decimal_type does
        const auto tenth {unit / 10U};
        auto rounded {form == 2 ? scaled.significand / (short_unit * 10U) :
                      form == 1 ? scaled.significand / short_unit : scaled.significand / (short_unit / 10U)};
        const auto sticky {scaled.inexact || rounded * tenth != scaled.significand};
        const auto rounded_exponent {exponent - 1 + fenv_round<DecimalType>(rounded, sign, sticky)};
        return pack_in_range<DecimalType>(rounded, rounded_exponent, sign);
    }
    return from_shortest<DecimalType>(val, floating_point_to_fd(val));
}

// The binary exponents whose 10^k is in cache_128, which also keep floor_log10_pow2 in the range of int
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int min_cached_binary_exponent {-1080};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int max_cached_binary_exponent {1130};

static_assert(dragonbox::floor_log10_pow2(-max_cached_binary_exponent) >= fast_float::significand_table_min_power &&
              dragonbox::floor_log10_pow2(-min_cached_binary_exponent) <= fast_float::significand_table_max_power,
              "The cached binary exponents need 10^k in the table");

// The products of a significand and cache_128 below are off by less than this many units of their fraction
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t product_error {4U};

struct next_digit_result
{
    std::uint64_t digit;
    std::uint64_t distance;
    bool rounded_up;
};

// The digit after the integer part of x * 10^k, rounded to nearest, from the fraction in units of 2^-fraction_bits.
// The distance is from 10 times the fraction to the digit, and is off by less than 10 * product_error.
template <unsigned fraction_bits>
BOOST_DECIMAL_CXX20_CONSTEXPR auto next_digit(const std::uint64_t fraction) noexcept -> next_digit_result
{
    constexpr auto one {UINT64_C(1) << fraction_bits};
    const auto tenths {dragonbox::umul128(fraction, 10U)};
    const auto rest {static_cast<std::uint64_t>(tenths) & (one - 1U)};
    const auto rounded_up {rest >= one / 2U};
    return {static_cast<std::uint64_t>(tenths >> fraction_bits) + static_cast<std::uint64_t>(rounded_up),
            rounded_up ? one - rest : rest, rounded_up};
}

// The shortest digits of an x87 or binary128 number from its bits, in the format of Float
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto wide_shortest(const int128::uint128_t& bits) noexcept -> ryu::floating_decimal_128
{
    constexpr bool x87 {std::numeric_limits<Float>::digits == 64};
    return ryu::generic_binary_to_decimal(bits, x87 ? 64U : 112U, 15U, x87);
}

// x * 10^k * 2^63 from the significand of an x87 or binary128 number and 10^k rounded up, within 4 units, for
// 18 to 20 digits. The interval of x holds at most one decimal of 16 digits, and its half width is a shift of 10^k.
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_wide_scaled(const Float& val) noexcept -> DecimalType
{
    constexpr bool x87 {std::numeric_limits<Float>::digits == 64};
    constexpr int bottom_bits {x87 ? 0 : 49};
    const auto bits {wide_bits(val)};
    const auto biased {wide_biased_exponent<Float>(bits)};
    const auto sign {wide_sign<Float>(bits)};
    const auto top {x87 ? bits.low : (UINT64_C(1) << 63U) | (bits.high << 15U) | (bits.low >> 49U)};
    const auto bottom {x87 ? UINT64_C(0) : bits.low & ((UINT64_C(1) << 49U) - 1U)};
    const auto binary_exponent {biased - wide_exponent_bias - 63};
    // The constructors leave the x87 infinity and NaN here, as Clang reads the bits through the stack and that stalls.
    // An unnormal number is a NaN as for the hardware.
    if (x87 && x87_not_finite(biased, top))
    {
        constexpr auto infinity {std::numeric_limits<DecimalType>::infinity()};
        const auto nan {decode_wide<Float>(bits).nan};
        return nan ? std::numeric_limits<DecimalType>::quiet_NaN() : sign ? -infinity : infinity;
    }
    // Ryu takes the bits where it decides, as GCC else keeps an x87 copy of x for it, and the reads of the bits above
    // wait on that store
    if (biased == 0 || top < (UINT64_C(1) << 63U) ||
        binary_exponent < min_cached_binary_exponent || binary_exponent > max_cached_binary_exponent)
    {
        return from_shortest<DecimalType>(wide_shortest<Float>(bits));
    }
    const bool power_of_two {top == (UINT64_C(1) << 63U) && bottom == 0U};
    const auto k {dragonbox::floor_log10_pow2(-binary_exponent)};

    // The 192-bit product (top + bottom / 2^49) * 10^k, without its lowest word, shifted right by 64 to 68 bits
    const auto cache {dragonbox::cache_128(k)};
    const auto shift {127 - dragonbox::floor_log2_pow10(k) - binary_exponent};
    BOOST_DECIMAL_ASSERT(shift >= 127 && shift <= 131);
    const auto high {dragonbox::umul128(top, cache.high)};
    const auto middle {dragonbox::umul128(top, cache.low).high};
    const auto extra {(dragonbox::umul128(bottom, cache.high) << 15U) + (dragonbox::umul128(bottom, cache.low) >> 49U)};
    int128::uint128_t product {high.high, high.low + middle};
    product.high += static_cast<std::uint64_t>(product.low < middle);
    product += extra.high;
    const auto value {product >> (shift - 127)};
    const auto scaled {static_cast<std::uint64_t>(value >> 63U)};
    const auto fraction {static_cast<std::uint64_t>(value) & ((UINT64_C(1) << 63U) - 1U)};
    const std::uint64_t half_width {cache.high >> (shift - 126 + bottom_bits)};

    constexpr bool keep_shortest {precision_v<DecimalType> > 20};
    // scaled has 18 to 20 digits, so constant divisors split them
    constexpr int short_dropped {keep_shortest ? 0 : 18 - precision_v<DecimalType>};
    constexpr auto short_unit {pow10(static_cast<std::uint64_t>(short_dropped))};
    constexpr auto top_digits {pow10(static_cast<std::uint64_t>(keep_shortest ? 0 : precision_v<DecimalType>))};
    const auto truncated {scaled / short_unit};
    const auto more {keep_shortest ? 0 : truncated >= top_digits * 10U ? 2 : truncated >= top_digits ? 1 : 0};
    const auto dropped {short_dropped + more};
    const auto unit {more == 2 ? short_unit * 100U : more == 1 ? short_unit * 10U : short_unit};
    const auto digits {more == 2 ? truncated / 100U : more == 1 ? truncated / 10U : truncated};
    const auto remainder {scaled - digits * unit};
    const auto exponent {dropped - k};
    if (exponent < std::numeric_limits<DecimalType>::min_exponent10 ||
        exponent + precision_v<DecimalType> - 1 > std::numeric_limits<DecimalType>::max_exponent10)
    {
        return from_shortest<DecimalType>(wide_shortest<Float>(bits));
    }

    // Distances below 2^63 to the decimal below, the one above and the midpoint, each within product_error
    constexpr auto one {UINT64_C(1) << 63U};
    // The interval of a power of two is half as wide below x, so the distance below counts double. This x has one
    // truncation and the error of 10^k, each below 1 unit, so the doubled distance stays within product_error.
    const auto below {remainder != 0U ? one : !power_of_two ? fraction : fraction < one / 2U ? fraction * 2U : one};
    const auto above {remainder == unit - 1U ? one - fraction : one};
    const auto midpoint {remainder == unit / 2U ? fraction : remainder == unit / 2U - 1U ? one - fraction : one};
    // One test of the nearer distance, as a short decimal is on either side of x with equal chance
    if ((above < below ? above : below) + product_error < half_width)
    {
        // A fast type keeps all its digits. A carry to 10^precision has one digit more, which its constructor takes.
        const auto inside {digits + static_cast<std::uint64_t>(above < below)};
        if (is_fast_type_v<DecimalType> && !keep_shortest)
        {
            return pack_in_range<DecimalType>(inside, exponent, sign);
        }
        const auto carry {!keep_shortest && inside == top_digits};
        return pack_fitting<DecimalType>(carry ? inside / 10U : inside, exponent + static_cast<int>(carry), sign);
    }
    if (keep_shortest && below > half_width + product_error && above > half_width + product_error)
    {
        // No integer is inside, so the shortest digits have one digit more and are the nearest to x. The half width,
        // rounded down, is off by less than 1 unit, which is 10 units more of error next to 10 times the fraction.
        const auto next {next_digit<63>(fraction)};
        if (next.distance + 10U * product_error < one / 2U &&
            (half_width > one / 8U || next.distance + 10U * (product_error + 1U) <
                                      half_width * (power_of_two && !next.rounded_up ? 5U : 10U)))
        {
            // A fast type counts the digits of scaled, as GCC calls the count of 128 bits out of line
            BOOST_DECIMAL_IF_CONSTEXPR (is_fast_type_v<DecimalType>)
            {
                const auto room {precision_v<DecimalType> - 1 - num_digits(scaled)};
                return pack_in_range<DecimalType>((int128::uint128_t {scaled} * 10U + next.digit) *
                                                  pow10(static_cast<int128::uint128_t>(room)), -k - 1 - room, sign);
            }
            return pack_fitting<DecimalType>(int128::uint128_t {scaled} * 10U + next.digit, -k - 1, sign);
        }
    }
    if (keep_shortest || below <= half_width + product_error || above <= half_width + product_error)
    {
        return from_shortest<DecimalType>(wide_shortest<Float>(bits));
    }
    // The digits and the one after them round with a sticky fraction, as to_decimal_type does. Next to the midpoint,
    // an exact compare gives the side of x: the digit 4 rounds as below it, and the 5 without the sticky bit as on it.
    constexpr auto tenth_unit {keep_shortest ? UINT64_C(1) : short_unit / 10U};
    auto rounded {more == 2 ? truncated / 10U : more == 1 ? truncated : scaled / tenth_unit};
    auto sticky {true};
    if (midpoint <= product_error)
    {
        rounded = digits * 10U + 5U;
        const auto significand {(int128::uint128_t {top} << bottom_bits) | bottom};
        const auto order {fast_float::compare_exact(exponent - 1, rounded, significand, binary_exponent - bottom_bits)};
        rounded -= static_cast<std::uint64_t>(order > 0);
        sticky = order != 0;
    }
    const auto rounded_exponent {exponent - 1 + fenv_round<DecimalType>(rounded, sign, sticky)};
    return pack_in_range<DecimalType>(rounded, rounded_exponent, sign);
}

// x * 10^k * 2^14 from the 113 significand bits of a binary128 number and 10^k rounded up, within 4 units, for 33 to
// 35 digits. As for x87, the shortest digits are the integer in the interval, or else the nearest with one digit more.
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_quad_scaled(const Float val) noexcept -> DecimalType
{
    const auto bits {bit_cast<int128::uint128_t>(val)};
    const auto biased {wide_biased_exponent<Float>(bits)};
    const auto sign {wide_sign<Float>(bits)};
    const auto top {(bits.high & ((UINT64_C(1) << 48U) - 1U)) | (UINT64_C(1) << 48U)};
    const auto binary_exponent {biased - wide_exponent_bias - 112};
    if (biased == 0 || binary_exponent < min_cached_binary_exponent || binary_exponent > max_cached_binary_exponent)
    {
        return from_shortest<DecimalType>(wide_shortest<Float>(bits));
    }
    const bool power_of_two {top == (UINT64_C(1) << 48U) && bits.low == 0U};
    const auto k {dragonbox::floor_log10_pow2(-binary_exponent)};

    // The 241-bit product of the significand and 10^k, without its lowest word, shifted right by 49 to 55 bits
    const auto cache {dragonbox::cache_128(k)};
    const auto shift {127 - dragonbox::floor_log2_pow10(k) - binary_exponent - 78};
    BOOST_DECIMAL_ASSERT(shift >= 49 && shift <= 53);
    const auto cross {dragonbox::umul128(top, cache.low) + dragonbox::umul128(bits.low, cache.low).high};
    const auto middle {dragonbox::umul128(bits.low, cache.high) + cross};
    const auto carry {static_cast<std::uint64_t>(middle < cross)};
    const auto high {dragonbox::umul128(top, cache.high) + int128::uint128_t {carry, 0U}};
    const auto value {(high << (64 - shift)) + (middle >> shift)};
    const auto scaled {value >> 14U};
    const auto fraction {static_cast<std::uint64_t>(value) & ((UINT64_C(1) << 14U) - 1U)};
    const std::uint64_t half_width {cache.high >> (shift + 1)};

    constexpr auto one {UINT64_C(1) << 14U};
    const auto below {!power_of_two ? fraction : fraction < one / 2U ? fraction * 2U : one};
    const auto above {one - fraction};
    if (below + product_error < half_width || above + product_error < half_width)
    {
        // 34 digits fit, and the digits of a fast type need no trailing zeros removed
        const auto integer {scaled + static_cast<std::uint64_t>(above < below)};
        if (integer < ten_34)
        {
            if (is_fast_type_v<DecimalType> || last_digit(integer) != 0U)
            {
                return pack_in_range<DecimalType>(integer, -k, sign);
            }
            const auto trimmed {remove_trailing_zeros_split(integer)};
            return pack_in_range<DecimalType>(trimmed.trimmed_number,
                                              static_cast<int>(trimmed.number_of_removed_zeros) - k, sign);
        }

        // A tie of 35 digits needs the side of x, which the error hides if x is close to the integer
        if ((above < below ? above : below) <= product_error && last_digit(integer) == 5U)
        {
            return from_shortest<DecimalType>(wide_shortest<Float>(bits));
        }
        const auto trimmed {remove_trailing_zeros_split(integer)};
        const auto exponent {static_cast<std::int32_t>(trimmed.number_of_removed_zeros) - k};
        return from_shortest<DecimalType>(ryu::floating_decimal_128 {trimmed.trimmed_number, exponent, sign,
                                                                          above < below ? -1 : 1});
    }
    if (below > half_width + product_error && above > half_width + product_error)
    {
        // The nearest decimal with one digit more, as for x87, and the side of x as above
        const auto next {next_digit<14>(fraction)};
        if (next.distance + 10U * product_error < one / 2U &&
            next.distance + 10U * (product_error + 1U) < half_width * (power_of_two && !next.rounded_up ? 5U : 10U) &&
            (next.distance > 10U * product_error || next.digit != 5U))
        {
            // With 34 digits, x rounds once from the truncated next digit and a sticky bit. fenv_round needs only the
            // parity of the digits, so it rounds a stand-in of two digits; one below 10^34 - 1 cannot carry.
            if (scaled >= ten_33 && scaled < ten_34 - 1U)
            {
                const auto parity {scaled.low & 1U};
                auto stand_in {parity * 10U + next.digit - static_cast<std::uint64_t>(next.rounded_up)};
                fenv_round<DecimalType>(stand_in, sign, true);
                return pack_in_range<DecimalType>(scaled + (stand_in - parity), -k, sign);
            }
            return from_shortest<DecimalType>(ryu::floating_decimal_128 {scaled * 10U + next.digit, -k - 1, sign,
                                                                              next.rounded_up ? -1 : 1});
        }
    }
    return from_shortest<DecimalType>(wide_shortest<Float>(bits));
}

// The significand bits, which numeric_limits does not give for __float128 on the libstdc++ of GCC 15 and before
template <typename Float>
struct significand_digits : std::integral_constant<int, std::numeric_limits<Float>::digits> {};

#ifdef BOOST_DECIMAL_HAS_FLOAT128
template <>
struct significand_digits<__float128> : std::integral_constant<int, 113> {};
#endif

// The tests do not draw -Wfloat-equal, but a signaling NaN can raise the invalid exception of the binary environment.
// The infinity comes through double because numeric_limits<__float128> is empty on old libstdc++.
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto is_finite_value(const Float val, std::false_type) noexcept -> bool
{
    constexpr Float inf {infinity_value<Float>()};
    return !is_nan_value(val) && val < inf && val > -inf;
}

// binary128 compares in software, so its exponent bits decide
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto is_finite_value(const Float val, std::true_type) noexcept -> bool
{
    return wide_biased_exponent<Float>(bit_cast<int128::uint128_t>(val)) != wide_max_biased_exponent;
}

// Whether a constructor gives val to from_float: a finite value, or each x87 number, as from_float picks
// from_wide_scaled for it, which tests the number itself. Clang reads its bits through the stack, and each read waits.
template <typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto takes_from_float(const Float val) noexcept -> bool
{
    if (constant_test<sizeof(Float) == 16 && significand_digits<Float>::value == 64>())
    {
        return true;
    }
    return is_finite_value(val, std::integral_constant<bool, sizeof(Float) == 16 &&
                                                             significand_digits<Float>::value == 113> {});
}

// The decimal infinity or quiet NaN of a binary number that takes_from_float leaves to the constructor
template <typename DecimalType, typename Float>
constexpr auto special_value(const Float val) noexcept -> DecimalType
{
    if (is_nan_value(val))
    {
        return std::numeric_limits<DecimalType>::quiet_NaN();
    }
    // The sign bit, as Clang can compute val > 0 before the test for a NaN, and that compare raises FE_INVALID
    #if defined(__GNUC__) || defined(__clang__)
    const bool negative {__builtin_signbit(val) != 0};
    #else
    const bool negative {val < 0};
    #endif
    return negative ? -std::numeric_limits<DecimalType>::infinity() : std::numeric_limits<DecimalType>::infinity();
}

// The conversion that from_float picks for a binary type and a decimal type
enum class float_path
{
    shortest,
    ryu,
    narrow,
    wide,
    quad,
    fast_shortest
};

template <float_path path>
using float_path_tag = std::integral_constant<float_path, path>;

// Dragonbox gives the shortest digits
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float val, float_path_tag<float_path::shortest>) noexcept -> DecimalType
{
    return from_shortest<DecimalType>(val, floating_point_to_fd(val));
}

// Ryu gives the shortest digits and their side of the exact value
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float val, float_path_tag<float_path::ryu>) noexcept -> DecimalType
{
    return from_shortest<DecimalType>(floating_point_to_fd(val));
}

// A float or a double to 7 digits
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float val, float_path_tag<float_path::narrow>) noexcept -> DecimalType
{
    return from_double_narrow<DecimalType>(val);
}

// An x87 number to 7, 16 or 34 digits, or a binary128 number to 7 or 16 digits
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float& val, float_path_tag<float_path::wide>) noexcept -> DecimalType
{
    return from_wide_scaled<DecimalType>(val);
}

// A binary128 number to 34 digits
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float val, float_path_tag<float_path::quad>) noexcept -> DecimalType
{
    return from_quad_scaled<DecimalType>(val);
}

// A fast type keeps all its digits, so Dragonbox can leave the trailing zeros
template <typename DecimalType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const float val, float_path_tag<float_path::fast_shortest>) noexcept -> DecimalType
{
    return from_shortest<DecimalType>(val, dragonbox::floating_point_to_fd32<false>(val));
}

template <typename DecimalType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const double val, float_path_tag<float_path::fast_shortest>) noexcept -> DecimalType
{
    return from_shortest<DecimalType>(val, dragonbox::floating_point_to_fd64<false>(val));
}

// A float or double to 7 digits, or 64 or 113 significand bits to 7 or 16 digits, rounds once unless near a decimal.
// To 34 digits, an x87 or binary128 number keeps its shortest digits: the integer in its interval, or one digit more.
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float& val) noexcept -> DecimalType
{
    using components = decltype(floating_point_to_fd(val));
    constexpr int digits {significand_digits<Float>::value};
    constexpr bool wide {std::is_same<components, ryu::floating_decimal_128>::value && sizeof(Float) == 16 &&
                         (digits == 64 || digits == 113)};
    constexpr bool narrow {std::is_same<components, dragonbox::floating_decimal<std::uint64_t>>::value ||
                           (std::is_same<components, dragonbox::floating_decimal<std::uint32_t>>::value &&
                            digits == 24)};
    constexpr bool plain {std::is_same<Float, float>::value || std::is_same<Float, double>::value};
    constexpr float_path path {precision_v<DecimalType> == 7 && narrow ? float_path::narrow :
                               is_fast_type_v<DecimalType> && plain ? float_path::fast_shortest :
                               wide && precision_v<DecimalType> <= 16 ? float_path::wide :
                               wide && precision_v<DecimalType> > 20 ?
                                   (digits == 64 ? float_path::wide : float_path::quad) :
                               std::is_same<components, ryu::floating_decimal_128>::value ? float_path::ryu :
                                                                                            float_path::shortest};
    static_assert(path == float_path::wide || sizeof(Float) != 16 || digits != 64,
                  "takes_from_float leaves the test of an x87 number to from_wide_scaled");
    return from_float<DecimalType>(val, float_path_tag<path> {});
}

} // namespace binary_to_decimal
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP
