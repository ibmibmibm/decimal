// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt


#ifndef BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP
#define BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/binary_bits.hpp>
#include <boost/decimal/detail/fast_float/compute_float64.hpp>
#include <boost/decimal/detail/fast_float/compute_float80_128.hpp>
#include <boost/decimal/detail/concepts.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cerrno>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

// Duplicated branches when on a machine with 64 bit long doubles
#if defined(__GNUC__) && __GNUC__ >= 6
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wduplicated-branches"
#endif
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace detail {

#if defined(BOOST_DECIMAL_HAS_FLOAT16) || defined(BOOST_DECIMAL_HAS_BRAINFLOAT16)

// A 16-bit type rounds the bits of the double again. Only a double on a midpoint of TargetType can be wrong, so the
// exact value sig * 10^q gives the side of that midpoint. The bits are built here, as a cast can round in other ways.
template <typename TargetType>
constexpr auto round_to_narrow(const double value, const std::int64_t q, const int128::uint128_t sig) noexcept
    -> TargetType
{
    constexpr int binary_digits {std::numeric_limits<TargetType>::digits};
    constexpr int min_exponent {std::numeric_limits<TargetType>::min_exponent - 1};
    constexpr int max_exponent {std::numeric_limits<TargetType>::max_exponent - 1};
    constexpr std::uint64_t infinity_bits {static_cast<std::uint64_t>(2 * max_exponent + 1) << (binary_digits - 1)};

    const auto bits {bit_cast<std::uint64_t>(value)};
    const auto exponent_bits {static_cast<int>((bits >> binary64_fraction_bits) & binary64_max_biased_exponent)};
    const auto exponent {exponent_bits - binary64_exponent_bias};
    constexpr auto leading_bit {UINT64_C(1) << binary64_fraction_bits};
    const auto significand {(bits & (leading_bit - 1U)) | leading_bit};

    // The bits of the double below the last place of TargetType, more for its subnormal values
    const auto low_exponent {exponent > min_exponent ? exponent : min_exponent};
    const auto dropped {low_exponent - (binary_digits - 1) - (exponent - binary64_fraction_bits)};
    std::uint64_t magnitude {};
    if (exponent > max_exponent)
    {
        magnitude = infinity_bits;
    }
    else if (exponent_bits != 0 && dropped <= 53)
    {
        const auto half {UINT64_C(1) << (dropped - 1)};
        const auto rest {significand & ((half << 1U) - 1U)};
        auto kept {significand >> dropped};
        const auto order {rest != half ? (rest > half ? 1 : -1) :
                          fast_float::compare_exact(q, sig, significand, exponent - binary64_fraction_bits)};
        if (order > 0 || (order == 0 && (kept & 1U) != 0U))
        {
            ++kept;
        }
        magnitude = (static_cast<std::uint64_t>(low_exponent - min_exponent) << (binary_digits - 1)) + kept;
        magnitude = magnitude < infinity_bits ? magnitude : infinity_bits;
    }
    return bit_cast<TargetType>(static_cast<std::uint16_t>(((bits >> 48U) & 0x8000U) | magnitude));
}

#endif

// A coefficient with 18 trailing zeros is at most a word times 5^18 * 2^18. The inverse of 5^18 modulo 2^64 gives the
// word, and one product makes sure of it. It gives 0 for any other coefficient.
BOOST_DECIMAL_NO_INLINE constexpr auto word_above_low_zeros(const int128::uint128_t sig) noexcept -> std::uint64_t
{
    const auto odd {sig >> 18U};
    const auto word {odd.low * UINT64_C(13940676669711910505)};
    return fast_float::umul128(word, UINT64_C(3814697265625)) == odd ? word : UINT64_C(0);
}

// The format of a binary target from numeric_limits, which has nothing for __float128 with libc++ or old libstdc++
template <typename T>
struct target_format
{
    static constexpr bool from_limits {true};
    static constexpr int digits {std::numeric_limits<T>::digits};
    static constexpr int min_exponent {std::numeric_limits<T>::min_exponent};
    static constexpr int max_exponent {std::numeric_limits<T>::max_exponent};
};

#if defined(BOOST_DECIMAL_HAS_FLOAT128) && !defined(BOOST_DECIMAL_LDBL_IS_FLOAT128)
template <>
struct target_format<__float128>
{
    static constexpr bool from_limits {false};
    static constexpr int digits {113};
    static constexpr int min_exponent {-16381};
    static constexpr int max_exponent {16384};
};
#endif

// An infinity or a NaN of a binary target, from numeric_limits or else from the bits of binary128
template <typename T>
constexpr auto target_special(const bool nan, const bool signaling, const bool negative, std::true_type) noexcept -> T
{
    return !nan ? (negative ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity()) :
           signaling ? std::numeric_limits<T>::signaling_NaN() : std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto target_special(const bool nan, const bool signaling, const bool negative,
                                                  std::false_type) noexcept -> T
{
    const std::uint64_t high {!nan ? (static_cast<std::uint64_t>(negative) << 63U) | UINT64_C(0x7FFF000000000000) :
                              signaling ? UINT64_C(0x7FFF400000000000) : UINT64_C(0x7FFF800000000000)};
    return bit_cast<T>(int128::uint128_t {high, 0U});
}

// A type wider than double in the format of long double
template <typename TargetType, typename Significand>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_wide_float(const std::int64_t q, const Significand& sig, const bool negative,
                                                 std::true_type) noexcept -> TargetType
{
    return static_cast<TargetType>(fast_float::compute_float80_128(q, sig, negative));
}

// A binary128 type, such as std::float128_t where long double has another format
template <typename TargetType, typename Significand>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_wide_float(const std::int64_t q, const Significand& sig, const bool negative,
                                                 std::false_type) noexcept -> TargetType
{
    using limits = target_format<TargetType>;
    const auto parts {fast_float::compute_long_double_parts<limits::digits, limits::min_exponent, limits::max_exponent>(
        q, static_cast<int128::uint128_t>(sig))};
    return fast_float::to_binary128<TargetType>(parts, negative);
}

} // namespace detail

template <typename Decimal, typename TargetType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_float(const Decimal& val) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, Decimal, detail::is_floating_point_v, TargetType, TargetType)
{
    using limits = detail::target_format<TargetType>;
    if (!isfinite(val))
    {
        return detail::target_special<TargetType>(!isinf(val), issignaling(val), val.isneg(),
                                                  std::integral_constant<bool, limits::from_limits> {});
    }

    // The casts to result are redundant, but in pre C++17 modes MSVC warns about implicit conversions
    TargetType result {};

    auto sig {val.full_significand()};
    auto exp {val.biased_exponent()};
    // A coefficient above the precision, which a raw encoding can hold, is zero as IEEE 754 reads the encoding
    if (sig == 0U || sig > detail::max_significand_v<Decimal>)
    {
        return val.isneg() ? -TargetType{0} : TargetType{0};
    }

    // A short value of a fast type has zeros that hide it from the paths for 64 bits, and few other values have them
    const auto full_sig {static_cast<int128::uint128_t>(sig)};
    if (BOOST_DECIMAL_UNLIKELY(detail::is_fast_type_v<Decimal> && full_sig.high != 0U &&
                               (full_sig.low & ((UINT64_C(1) << 18U) - 1U)) == 0U))
    {
        const auto word {detail::word_above_low_zeros(full_sig)};
        if (word != 0U)
        {
            sig = static_cast<decltype(sig)>(word);
            exp += 18;
        }
    }

    // A wider type than double is rounded in the format of long double, or in binary128
    constexpr bool long_double_format {limits::digits == std::numeric_limits<long double>::digits &&
                                       limits::max_exponent == std::numeric_limits<long double>::max_exponent};
    static_assert(limits::digits != 0 &&
                  (limits::digits <= std::numeric_limits<double>::digits || long_double_format ||
                   (limits::digits == 113 && limits::max_exponent == 16384 && sizeof(TargetType) == 16)),
                  "to_float supports the formats up to double, the format of long double, and binary128");

    BOOST_DECIMAL_IF_CONSTEXPR (limits::digits > std::numeric_limits<double>::digits)
    {
        // The tag also takes the narrow types, which compile this branch in C++14
        constexpr bool in_long_double {long_double_format || limits::digits <= std::numeric_limits<double>::digits};
        result = detail::to_wide_float<TargetType>(exp, sig, val.isneg(),
                                                   std::integral_constant<bool, in_long_double> {});
    }
    else
    {
        // Rounds the exact value once, also to a subnormal value or to infinity
        constexpr bool float_format {limits::digits == std::numeric_limits<float>::digits &&
                                     limits::max_exponent == std::numeric_limits<float>::max_exponent};
        using binary_type = std::conditional_t<float_format, float, double>;
        const auto wide_sig {static_cast<int128::uint128_t>(sig)};
        binary_type value {};
        if (wide_sig.high != 0U)
        {
            value = detail::fast_float::to_binary<binary_type>(detail::fast_float::compute_float_wide<binary_type>(exp, wide_sig), val.isneg());
        }
        else if (!detail::fast_float::clinger_fast_path(exp, wide_sig.low, val.isneg(), value))
        {
            value = detail::fast_float::to_binary<binary_type>(detail::fast_float::compute_float<binary_type>(exp, wide_sig.low), val.isneg());
        }
        #if defined(BOOST_DECIMAL_HAS_FLOAT16) || defined(BOOST_DECIMAL_HAS_BRAINFLOAT16)
        if constexpr (limits::digits < std::numeric_limits<float>::digits)
        {
            result = detail::round_to_narrow<TargetType>(value, exp, wide_sig);
        }
        else
        #endif
        {
            result = static_cast<TargetType>(value);
        }
    }

    return result;
}

#if defined(__GNUC__) && __GNUC__ >= 6
#  pragma GCC diagnostic pop
#endif
#ifdef _MSC_VER
#  pragma warning(pop)
#endif

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP
