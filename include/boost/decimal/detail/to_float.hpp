// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt


#ifndef BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP
#define BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/dragonbox/dragonbox.hpp>
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

// A type narrower than float rounds the double again, which is wrong only for a double on its midpoint that is not
// the exact value sig * 10^q. The double then moves one unit toward the exact value, so the second rounding agrees.
template <typename TargetType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto away_from_midpoint(const double value, const std::int64_t q,
                                                      const int128::uint128_t sig) noexcept -> double
{
    constexpr int binary_digits {std::numeric_limits<TargetType>::digits};
    constexpr int min_exponent {std::numeric_limits<TargetType>::min_exponent - 1};

    const auto bits {bit_cast<std::uint64_t>(value)};
    const auto exponent {static_cast<int>((bits >> 52U) & 0x7FFU) - 1023};
    const auto significand {(bits & ((UINT64_C(1) << 52U) - 1U)) | (UINT64_C(1) << 52U)};

    // The bits of the double below the last place of TargetType, more for its subnormal values
    const auto dropped {(exponent > min_exponent ? exponent : min_exponent) - (binary_digits - 1) - (exponent - 52)};
    if (dropped > 53 || (significand & ((UINT64_C(2) << (dropped - 1)) - 1U)) != (UINT64_C(1) << (dropped - 1)))
    {
        return value;
    }

    const auto order {fast_float::compare_exact(q, sig, significand, exponent - 52)};
    return order == 0 ? value : bit_cast<double>(order > 0 ? bits + 1U : bits - 1U);
}

// A coefficient with 18 trailing zeros is at most a word times 5^18 * 2^18. The inverse of 5^18 modulo 2^64 gives the
// word, and one product makes sure of it. It gives 0 for any other coefficient.
BOOST_DECIMAL_NO_INLINE constexpr auto word_above_low_zeros(const int128::uint128_t sig) noexcept -> std::uint64_t
{
    const auto odd {sig >> 18U};
    const auto word {odd.low * UINT64_C(13940676669711910505)};
    return dragonbox::umul128(word, UINT64_C(3814697265625)) == odd ? word : UINT64_C(0);
}

} // namespace detail

template <typename Decimal, typename TargetType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_float(const Decimal& val) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, Decimal, detail::is_floating_point_v, TargetType, TargetType)
{
    if (!isfinite(val))
    {
        if (isinf(val))
        {
            return val.isneg() ? -std::numeric_limits<TargetType>::infinity() : std::numeric_limits<TargetType>::infinity();
        }
        if (issignaling(val))
        {
            return std::numeric_limits<TargetType>::signaling_NaN();
        }
        return std::numeric_limits<TargetType>::quiet_NaN();
    }

    // The casts to result are redundant, but in pre C++17 modes MSVC warns about implicit conversions
    TargetType result {};

    auto sig {val.full_significand()};
    auto exp {val.biased_exponent()};
    if (sig == 0U)
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

    #if BOOST_DECIMAL_LDBL_BITS == 64
    constexpr bool binary64_long_double {true};
    #else
    constexpr bool binary64_long_double {false};
    #endif

    BOOST_DECIMAL_IF_CONSTEXPR (!std::is_same<TargetType, long double>::value || binary64_long_double)
    {
        // Rounds the exact value once, also to a subnormal value or to infinity
        using binary_type = std::conditional_t<std::is_same<TargetType, float>::value, float, double>;
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
        BOOST_DECIMAL_IF_CONSTEXPR (std::numeric_limits<TargetType>::digits < std::numeric_limits<float>::digits)
        {
            value = static_cast<binary_type>(detail::away_from_midpoint<TargetType>(static_cast<double>(value), exp, wide_sig));
        }
        result = static_cast<TargetType>(value);
    }
    else
    {
        result = static_cast<TargetType>(detail::fast_float::compute_float80_128(exp, sig, val.isneg()));
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
