// Copyright 2026 Pavel P
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP
#define BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/add_impl.hpp>
#include <boost/decimal/detail/attributes.hpp>
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
    bool inexact;
};

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

// x * 10^k for Dragonbox's k, 18 or 19 digits for a normal double; both results are exact, as they come from the
// product and integer check of Dragonbox's compute_mul_parity, whose cache precision is proven for every double
constexpr auto dragonbox_scaled(const std::uint64_t two_fc, const int binary_exponent) noexcept
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
    return {integer_part, minus_k, !is_integer};
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
        return to_decimal_type<DecimalType>(scaled.significand, scaled.exponent, scaled.inexact, sign);
    }

    const auto exact {exact_scaled(significand, binary_exponent)};
    return to_decimal_type<DecimalType>(exact.significand, exact.exponent, exact.inexact, sign);
}

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
        // A subnormal double has fewer digits, but Dragonbox takes its shortest digits from the same scale
        const auto bits {bit_cast<std::uint64_t>(static_cast<double>(val))};
        const auto exponent_bits {static_cast<int>((bits >> 52U) & 0x7FFU)};
        const auto fraction {bits & ((UINT64_C(1) << 52U) - 1U)};
        const auto significand {exponent_bits == 0 ? fraction : fraction | (UINT64_C(1) << 52U)};
        const auto scaled {dragonbox_scaled(significand * 2U, exponent_bits == 0 ? -1074 : exponent_bits - 1075)};
        const auto shortest {static_cast<std::uint64_t>(components.mantissa) *
                             pow10(static_cast<std::uint64_t>(components.exponent - scaled.exponent))};
        if (shortest != scaled.significand || scaled.inexact)
        {
            const std::uint64_t digits {scaled.significand >= shortest ? components.mantissa :
                                                                         components.mantissa - 1U};
            return to_decimal_type<DecimalType>(digits, components.exponent, true, components.sign);
        }
    }

    return DecimalType {components.mantissa, components.exponent, components.sign};
}

// 2^64 is 6 modulo 10, so the last digit needs no division of 128 bits
constexpr auto last_digit(const int128::uint128_t value) noexcept -> std::uint64_t
{
    return (value.low % 10U + 6U * (value.high % 10U)) % 10U;
}

// The wider types take the side of the exact value from Ryu
template <typename DecimalType, typename Float>
constexpr auto from_shortest(const Float, const ryu::floating_decimal_128& components) noexcept -> DecimalType
{
    if (components.direction != 0 && last_digit(components.mantissa) == 5U &&
        num_digits(components.mantissa) == precision_v<DecimalType> + 1)
    {
        const auto digits {components.direction > 0 ? components.mantissa : components.mantissa - 1U};
        return to_decimal_type<DecimalType>(digits, components.exponent, true, components.sign);
    }
    return DecimalType {components.mantissa, components.exponent, components.sign};
}

// Dragonbox or Ryu gives the shortest digits
template <typename DecimalType, typename Float>
BOOST_DECIMAL_CXX20_CONSTEXPR auto from_float(const Float val) noexcept -> DecimalType
{
    return from_shortest<DecimalType>(val, floating_point_to_fd(val));
}

} // namespace binary_to_decimal
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_BINARY_TO_DECIMAL_HPP
