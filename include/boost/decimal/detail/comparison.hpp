// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_COMPARISON_HPP
#define BOOST_DECIMAL_DETAIL_COMPARISON_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/apply_sign.hpp>
#include <boost/decimal/detail/normalize.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/to_decimal.hpp>
#include <boost/decimal/detail/shrink_significand.hpp>
#include <boost/decimal/detail/cmath/isfinite.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/u256.hpp>
#include "int128.hpp"
#include <boost/decimal/detail/attributes.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <cstdint>
#endif

namespace boost {
namespace decimal {

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4127)
#endif

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto equality_impl(DecimalType lhs, DecimalType rhs) noexcept -> bool
{
    using comp_type = typename DecimalType::significand_type;

    // Step 1: Check for NANs per IEEE 754
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(lhs) || !isfinite(rhs))
    {
        if (isnan(lhs) || isnan(rhs))
        {
            return false;
        }
        // The payload of an infinity has no meaning, so two infinities compare by sign
        if (isinf(lhs) && isinf(rhs))
        {
            return lhs.isneg() == rhs.isneg();
        }
    }
    #endif

    // Step 2: Fast path
    if (lhs.bits_ == rhs.bits_)
    {
        return true;
    }

    const auto lhs_components {lhs.to_components()};
    const auto rhs_components {rhs.to_components()};

    // References, because a copy of a 128-bit significand reloads it just after two 64-bit stores
    const auto& lhs_sig {lhs_components.sig};
    const auto& rhs_sig {rhs_components.sig};

    // Step 4: Check signs
    if (lhs_components.sign != rhs_components.sign)
    {
        return (lhs_sig == 0U && rhs_sig == 0U);
    }

    // Step 4b: Same-sign zeros from any cohort compare equal regardless of their exponents
    // (IEEE 754-2008 3.5.1). Without this, two zeros with delta_exp greater than the type's
    // precision would fall through to the early-return below and wrongly compare unequal.
    if (lhs_sig == 0U && rhs_sig == 0U)
    {
        return true;
    }

    // Step 5: Check the exponents
    // If the difference is greater than we can represent in the significand than we can assume they are different
    const auto lhs_exp {lhs_components.exp};
    const auto rhs_exp {rhs_components.exp};

    const auto delta_exp {lhs_exp - rhs_exp};

    if (BOOST_DECIMAL_UNLIKELY(delta_exp > detail::precision_v<DecimalType> || delta_exp < -detail::precision_v<DecimalType>))
    {
        return false;
    }

    // Step 6: Normalize the significand and compare
    BOOST_DECIMAL_IF_CONSTEXPR (detail::decimal_val_v<DecimalType> >= 128)
    {
        // The significand with the larger exponent is scaled when the product fits, which needs no
        // digit count for a small delta. Else it has more digits than any significand.
        constexpr auto max_delta_diff {std::numeric_limits<comp_type>::digits10 - detail::precision_v<DecimalType>};
        constexpr auto max_digits {std::numeric_limits<comp_type>::digits10};
        if (delta_exp > 0)
        {
            if (delta_exp > max_delta_diff && detail::num_digits(lhs_sig) + delta_exp > max_digits)
            {
                return false;
            }
            return lhs_sig * detail::pow10(static_cast<comp_type>(delta_exp)) == rhs_sig;
        }
        else if (delta_exp < 0)
        {
            if (delta_exp < -max_delta_diff && detail::num_digits(rhs_sig) - delta_exp > max_digits)
            {
                return false;
            }
            return lhs_sig == rhs_sig * detail::pow10(static_cast<comp_type>(-delta_exp));
        }

        return lhs_sig == rhs_sig;
    }
    else
    {
        using promoted_sig_type = std::conditional_t<std::is_same<typename DecimalType::significand_type, std::uint32_t>::value, std::uint64_t, int128::uint128_t>;
        promoted_sig_type promotes_lhs {lhs_sig};
        promoted_sig_type promotes_rhs {rhs_sig};

        if (delta_exp > 0)
        {
            promotes_lhs *= detail::pow10(static_cast<promoted_sig_type>(delta_exp));
        }
        else if (delta_exp < 0)
        {
            promotes_rhs *= detail::pow10(static_cast<promoted_sig_type>(-delta_exp));
        }

        return promotes_lhs == promotes_rhs;
    }
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

template <BOOST_DECIMAL_FAST_DECIMAL_FLOATING_TYPE DecimalType>
BOOST_DECIMAL_FORCE_INLINE constexpr auto fast_equality_impl(const DecimalType& lhs, const DecimalType& rhs) noexcept -> bool
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs) || isnan(rhs))
    {
        return false;
    }
    #endif

    const auto lhs_sig {lhs.significand_};
    const auto rhs_sig {rhs.significand_};

    if (lhs_sig == 0U && rhs_sig == 0U)
    {
        // -0 == +0
        return true;
    }

    return  lhs_sig == rhs_sig &&
            lhs.exponent_ == rhs.exponent_ &&
            lhs.sign_ == rhs.sign_;
}

template <BOOST_DECIMAL_FAST_DECIMAL_FLOATING_TYPE DecimalType>
BOOST_DECIMAL_FORCE_INLINE constexpr auto fast_inequality_impl(const DecimalType& lhs, const DecimalType& rhs) noexcept -> bool
{
    return !fast_equality_impl(lhs, rhs);
}

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4127) // BOOST_DECIMAL_IF_CONSTEXPR prior to C++17 is constant conditional expression
#endif

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType = decimal32_t, BOOST_DECIMAL_INTEGRAL T1,
          BOOST_DECIMAL_INTEGRAL U1, BOOST_DECIMAL_INTEGRAL T2, BOOST_DECIMAL_INTEGRAL U2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto equal_parts_impl(T1 lhs_sig, U1 lhs_exp, bool lhs_sign,
                                T2 rhs_sig, U2 rhs_exp, bool rhs_sign) noexcept -> std::enable_if_t<detail::is_ieee_type_v<DecimalType>, bool>
{
    using comp_type = detail::make_unsigned_t<std::conditional_t<(std::numeric_limits<T1>::digits10 > std::numeric_limits<T2>::digits10), T1, T2>>;

    BOOST_DECIMAL_ASSERT(lhs_sig >= 0U);
    BOOST_DECIMAL_ASSERT(rhs_sig >= 0U);

    if (lhs_sig == 0U && rhs_sig == 0U)
    {
        // +0 == -0
        return true;
    }

    // We con compare signs right away
    if (lhs_sign != rhs_sign)
    {
        return false;
    }

    const auto delta_exp {lhs_exp - rhs_exp};

    // Check the value of delta exp to avoid to large a value for pow10
    // Also if only one of the significands is 0 then we know the values have to be mismatched
    if (BOOST_DECIMAL_UNLIKELY(delta_exp > detail::precision_v<DecimalType> || delta_exp < -detail::precision_v<DecimalType>))
    {
        return false;
    }

    auto new_lhs_sig {static_cast<comp_type>(lhs_sig)};
    auto new_rhs_sig {static_cast<comp_type>(rhs_sig)};

    // Step 5: Normalize the significand and compare
    // Instead of multiplying the larger number, divide the smaller one
    BOOST_DECIMAL_IF_CONSTEXPR (detail::decimal_val_v<DecimalType> >= 128)
    {
        // Instead of multiplying the larger number, divide the smaller one
        //
        // We try for multiplication even though it's a small range
        // Since it's an order of magnitude faster
        if (delta_exp <= 4 && delta_exp >= -4)
        {
            if (delta_exp > 0)
            {
                new_lhs_sig *= detail::pow10(static_cast<comp_type>(delta_exp));
            }
            else if (delta_exp < 0)
            {
                new_rhs_sig *= detail::pow10(static_cast<comp_type>(-delta_exp));
            }
        }
        else if (delta_exp > 0)
        {
            // Check if we can divide rhs_sig safely
            // E.g. 9e0 != 90000000204928e-13 but if we just did division we would falsely get 9 ?= 9
            if (new_rhs_sig % detail::pow10(static_cast<comp_type>(delta_exp)) != 0U)
            {
                return false;
            }
            new_rhs_sig /= detail::pow10(static_cast<comp_type>(delta_exp));
        }
        else if (delta_exp < 0)
        {
            // Check if we can divide lhs_sig safely
            if (new_lhs_sig % detail::pow10(static_cast<comp_type>(-delta_exp)) != 0U)
            {
                return false;
            }
            new_lhs_sig /= detail::pow10(static_cast<comp_type>(-delta_exp));
        }

        return new_lhs_sig == new_rhs_sig;
    }
    else
    {
        using promoted_sig_type = std::conditional_t<std::is_same<comp_type, std::uint32_t>::value, std::uint64_t, int128::uint128_t>;
        promoted_sig_type promotes_lhs {new_lhs_sig};
        promoted_sig_type promotes_rhs {new_rhs_sig};

        if (delta_exp > 0)
        {
            promotes_lhs *= detail::pow10(static_cast<promoted_sig_type>(delta_exp));
        }
        else if (delta_exp < 0)
        {
            promotes_rhs *= detail::pow10(static_cast<promoted_sig_type>(-delta_exp));
        }

        return promotes_lhs == promotes_rhs;
    }
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType = decimal32_t, BOOST_DECIMAL_INTEGRAL T1,
          BOOST_DECIMAL_INTEGRAL U1, BOOST_DECIMAL_INTEGRAL T2, BOOST_DECIMAL_INTEGRAL U2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto equal_parts_impl(T1 lhs_sig, U1 lhs_exp, const bool lhs_sign,
                                T2 rhs_sig, U2 rhs_exp, const bool rhs_sign) noexcept -> std::enable_if_t<detail::is_fast_type_v<DecimalType>, bool>
{
    using comp_type = std::conditional_t<(std::numeric_limits<T1>::digits10 > std::numeric_limits<T2>::digits10), T1, T2>;

    BOOST_DECIMAL_ASSERT(lhs_sig >= 0U);
    BOOST_DECIMAL_ASSERT(rhs_sig >= 0U);

    if (lhs_sig == 0U && rhs_sig == 0U)
    {
        // +0 == -0
        return true;
    }

    auto new_lhs_sig {static_cast<comp_type>(lhs_sig)};
    auto new_rhs_sig {static_cast<comp_type>(rhs_sig)};

    detail::normalize<DecimalType>(new_lhs_sig, lhs_exp);
    detail::normalize<DecimalType>(new_rhs_sig, rhs_exp);

    #ifdef BOOST_DECIMAL_DEBUG_EQUAL
    std::cerr << "Normalized Values"
              << "\nlhs_sig: " << new_lhs_sig
              << "\nlhs_exp: " << lhs_exp
              << "\nrhs_sig: " << new_rhs_sig
              << "\nrhs_exp: " << rhs_exp << std::endl;
    #endif

    return lhs_sign == rhs_sign &&
           lhs_exp == rhs_exp &&
           new_lhs_sig == new_rhs_sig;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal, BOOST_DECIMAL_INTEGRAL Integer>
BOOST_DECIMAL_CUDA_CONSTEXPR auto mixed_equality_impl(Decimal lhs, Integer rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal> && detail::is_integral_v<Integer>), bool>
{
    using exp_type = typename Decimal::biased_exponent_type;

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs) || isinf(lhs))
    {
        return false;
    }
    #endif

    bool rhs_isneg {false};
    BOOST_DECIMAL_IF_CONSTEXPR (detail::is_signed_v<Integer>)
    {
        if (rhs < static_cast<Integer>(0))
        {
            rhs_isneg = true;
        }
    }

    const auto rhs_significand {detail::make_positive_unsigned(rhs)};

    return equal_parts_impl<Decimal>(lhs.full_significand(), lhs.biased_exponent(), lhs.isneg(),
                                     rhs_significand, static_cast<exp_type>(0), rhs_isneg);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto mixed_decimal_equality_impl(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    static_assert(!std::is_same<Decimal1, Decimal2>::value, "Equality of same type exists in simpler form");
    using Bigger_Decimal_Type = std::conditional_t<(sizeof(lhs) > sizeof(rhs)), Decimal1, Decimal2>;

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs) || isnan(rhs))
    {
        return false;
    }
    #endif

    const auto new_lhs = to_decimal<Bigger_Decimal_Type>(lhs);
    const auto new_rhs = to_decimal<Bigger_Decimal_Type>(rhs);

    return new_lhs == new_rhs;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator==(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    return mixed_decimal_equality_impl(lhs, rhs);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator!=(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    return !(mixed_decimal_equality_impl(lhs, rhs));
}

template <BOOST_DECIMAL_FAST_DECIMAL_FLOATING_TYPE DecimalType>
BOOST_DECIMAL_FORCE_INLINE constexpr auto fast_less_impl(const DecimalType& lhs, const DecimalType& rhs) noexcept -> bool
{
    // One test finds a zero or a value that is not finite in either operand, except in
    // two words, where the test costs more than the separate tests.
    using sig_type = decltype(lhs.significand_);
    constexpr auto max_sig {static_cast<sig_type>(detail::max_significand_v<DecimalType>)};
    if (detail::decimal_val_v<DecimalType> >= 128 ||
        BOOST_DECIMAL_UNLIKELY((static_cast<unsigned>(static_cast<sig_type>(lhs.significand_ - 1U) >= max_sig) |
                                static_cast<unsigned>(static_cast<sig_type>(rhs.significand_ - 1U) >= max_sig)) != 0U))
    {
    #ifndef BOOST_DECIMAL_FAST_MATH
        if (not_finite(lhs) || not_finite(rhs))
        {
            if (isnan(lhs) || isnan(rhs) ||
                (!lhs.isneg() && rhs.isneg()))
            {
                return false;
            }
            else if (lhs.isneg() && !rhs.isneg())
            {
                return true;
            }
            else if (isfinite(lhs) && isinf(rhs))
            {
                return !signbit(rhs);
            }
            else if (isinf(lhs) && isfinite(rhs))
            {
                return signbit(rhs);
            }
        }
    #endif

        // Needed to correctly compare signed and unsigned zeros.
        // Every zero compares equal to every other zero, so neither of a pair of them is less
        // (IEEE 754-2008 5.11). Ordering -0 before +0 is the job of totalorder, not of operator<.
        if (lhs.significand_ == 0U || rhs.significand_ == 0U)
        {
            if (lhs.significand_ == 0U && rhs.significand_ == 0U)
            {
                return false;
            }
            return lhs.significand_ == 0U ? !rhs.sign_ : lhs.sign_;
        }
    }

    BOOST_DECIMAL_IF_CONSTEXPR (detail::decimal_val_v<DecimalType> < 128)
    {
        // Normalized values order by exponent and then significand, so one signed key holds the
        // whole order. The key of a nonzero value is positive, so its negation is less than it.
        constexpr auto below_min_sig {detail::pow10(static_cast<std::uint64_t>(detail::precision_v<DecimalType> > 16 ? 0 : detail::precision_v<DecimalType> - 1)) - 1U};
        const auto lhs_mask {-static_cast<std::int64_t>(lhs.sign_)};
        const auto rhs_mask {-static_cast<std::int64_t>(rhs.sign_)};
        const auto lhs_key {(static_cast<std::int64_t>(lhs.exponent_) << 53) + static_cast<std::int64_t>(static_cast<std::uint64_t>(lhs.significand_) - below_min_sig)};
        const auto rhs_key {(static_cast<std::int64_t>(rhs.exponent_) << 53) + static_cast<std::int64_t>(static_cast<std::uint64_t>(rhs.significand_) - below_min_sig)};
        return ((lhs_key ^ lhs_mask) - lhs_mask) < ((rhs_key ^ rhs_mask) - rhs_mask);
    }
    else
    {
        // The same key in two words: the exponent sits above the 49 high bits of the significand
        const auto lhs_high {(static_cast<std::uint64_t>(lhs.exponent_) << 49) + int128::uint128_t{lhs.significand_}.high};
        const auto rhs_high {(static_cast<std::uint64_t>(rhs.exponent_) << 49) + int128::uint128_t{rhs.significand_}.high};
        const auto lhs_low {int128::uint128_t{lhs.significand_}.low};
        const auto rhs_low {int128::uint128_t{rhs.significand_}.low};
        const auto same_high {static_cast<unsigned>(lhs_high == rhs_high)};
        const auto above {static_cast<unsigned>(lhs_high > rhs_high) | (same_high & static_cast<unsigned>(lhs_low > rhs_low))};
        const auto below {static_cast<unsigned>(lhs_high < rhs_high) | (same_high & static_cast<unsigned>(lhs_low < rhs_low))};
        const auto lhs_negative {static_cast<unsigned>(lhs.sign_)};
        const auto signs_differ {lhs_negative ^ static_cast<unsigned>(rhs.sign_)};
        const auto same_sign {(lhs_negative & above) | (~lhs_negative & below & 1U)};
        return ((signs_differ & lhs_negative) | (~signs_differ & same_sign & 1U)) != 0U;
    }
}

template <BOOST_DECIMAL_INTEGRAL T, BOOST_DECIMAL_INTEGRAL U>
BOOST_DECIMAL_FORCE_INLINE constexpr auto fast_type_less_parts_impl(T lhs_sig, U lhs_exp, bool lhs_sign,
                                                                    T rhs_sig, U rhs_exp, bool rhs_sign) noexcept -> bool
{
    // A pair of zeros compares equal no matter what their signs are (IEEE 754-2008 5.11)
    if (lhs_sig == 0U || rhs_sig == 0U)
    {
        if (lhs_sig == 0U && rhs_sig == 0U)
        {
            return false;
        }
        return lhs_sig == 0U ? !rhs_sign : lhs_sign;
    }

    if (lhs_sign != rhs_sign)
    {
        return lhs_sign;
    }

    if (lhs_exp != rhs_exp)
    {
        return lhs_sign ? lhs_exp > rhs_exp : lhs_exp < rhs_exp;
    }

    return lhs_sign ? lhs_sig > rhs_sig : lhs_sig < rhs_sig;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto sequential_less_impl(DecimalType lhs, DecimalType rhs) noexcept -> bool
{
    using comp_type = std::conditional_t<detail::decimal_val_v<DecimalType> < 64, std::uint_fast64_t, int128::uint128_t>;

    // Non-finite handling is performed by the callers. With NaN ruled out,
    // identical bit patterns encode the same finite value, so lhs < rhs is false.
    if (lhs.bits_ == rhs.bits_)
    {
        return false;
    }

    // One combination-field branch per operand instead of three.
    const auto lhs_components {lhs.to_components()};
    const auto rhs_components {rhs.to_components()};

    const auto lhs_sign {lhs_components.sign};
    const auto rhs_sign {rhs_components.sign};

    // References, because a copy of a 128-bit significand reloads it just after two 64-bit stores
    const auto& lhs_sig {lhs_components.sig};
    const auto& rhs_sig {rhs_components.sig};

    // Signed zeros: +0 == -0, neither is less.
    if (lhs_sig == 0U || rhs_sig == 0U)
    {
        if (lhs_sig == rhs_sig)
        {
            return false;
        }
        return lhs_sig == 0U ? !rhs_sign : lhs_sign;
    }

    const auto lhs_exp {lhs_components.exp};
    const auto rhs_exp {rhs_components.exp};
    const auto delta_exp {lhs_exp - rhs_exp};

    BOOST_DECIMAL_IF_CONSTEXPR (detail::decimal_val_v<DecimalType> >= 128)
    {
        if (lhs_sign != rhs_sign)
        {
            return lhs_sign;
        }

        // The exponents of the leading digits decide first. When they are equal, the significand
        // with the larger exponent has fewer digits, and scaled to the digits of the other it fits.
        if (delta_exp != 0)
        {
            const auto lead_delta {delta_exp + detail::num_digits(lhs_sig) - detail::num_digits(rhs_sig)};
            if (lead_delta != 0)
            {
                return lhs_sign ? lead_delta > 0 : lead_delta < 0;
            }

            if (delta_exp > 0)
            {
                const auto scaled {lhs_sig * detail::pow10(static_cast<comp_type>(delta_exp))};
                return lhs_sign ? scaled > rhs_sig : scaled < rhs_sig;
            }

            const auto scaled {rhs_sig * detail::pow10(static_cast<comp_type>(-delta_exp))};
            return lhs_sign ? lhs_sig > scaled : lhs_sig < scaled;
        }

        return lhs_sign ? lhs_sig > rhs_sig : lhs_sig < rhs_sig;
    }

    // A gap of p or more decides by the exponents alone. Smaller gaps scale both sides from the table,
    // so a random gap costs no mispredicted branch, and the products have at most 18 or 35 digits.
    constexpr auto max_delta_diff {detail::decimal_val_v<DecimalType> < 64 ? 11 : 19};

    if (BOOST_DECIMAL_UNLIKELY(delta_exp > max_delta_diff || delta_exp < -max_delta_diff))
    {
        if (lhs_sign != rhs_sign)
        {
            return lhs_sign;
        }
        return lhs_sign ? lhs_exp > rhs_exp : lhs_exp < rhs_exp;
    }

    BOOST_DECIMAL_IF_CONSTEXPR (sizeof(comp_type) == sizeof(std::uint64_t))
    {
        const auto lhs_scaled {static_cast<std::uint64_t>(lhs_sig) * detail::pow10(static_cast<std::uint64_t>(delta_exp > 0 ? delta_exp : 0))};
        const auto rhs_scaled {static_cast<std::uint64_t>(rhs_sig) * detail::pow10(static_cast<std::uint64_t>(delta_exp < 0 ? -delta_exp : 0))};
        const auto lhs_mask {-static_cast<std::int64_t>(lhs_sign)};
        const auto rhs_mask {-static_cast<std::int64_t>(rhs_sign)};
        return ((static_cast<std::int64_t>(lhs_scaled) ^ lhs_mask) - lhs_mask) < ((static_cast<std::int64_t>(rhs_scaled) ^ rhs_mask) - rhs_mask);
    }

    // One product of the significand with the larger exponent. Only the bits decide the order,
    // because random signs and exponents make a branch mispredict.
    const auto swapped {static_cast<unsigned>(delta_exp < 0)};
    const auto larger {static_cast<std::uint64_t>(swapped != 0U ? rhs_sig : lhs_sig)};
    const auto smaller {static_cast<std::uint64_t>(swapped != 0U ? lhs_sig : rhs_sig)};
    const auto product {detail::umul128(larger, detail::pow10(static_cast<std::uint64_t>(swapped != 0U ? -delta_exp : delta_exp)))};
    const auto above {static_cast<unsigned>(product.high != 0U) | static_cast<unsigned>(product.low > smaller)};
    const auto below {static_cast<unsigned>(product.high == 0U) & static_cast<unsigned>(product.low < smaller)};
    const auto flip {swapped ^ static_cast<unsigned>(lhs_sign)};
    const auto same_sign {(flip & above) | (~flip & below & 1U)};
    const auto signs_differ {static_cast<unsigned>(lhs_sign != rhs_sign)};
    return ((signs_differ & static_cast<unsigned>(lhs_sign)) | (~signs_differ & same_sign & 1U)) != 0U;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType = decimal32_t, BOOST_DECIMAL_INTEGRAL T1,
        BOOST_DECIMAL_INTEGRAL U1, BOOST_DECIMAL_INTEGRAL T2, BOOST_DECIMAL_INTEGRAL U2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto less_parts_impl(T1 lhs_sig, U1 lhs_exp, bool lhs_sign,
                               T2 rhs_sig, U2 rhs_exp, bool rhs_sign, bool normalized = false) noexcept -> std::enable_if_t<detail::decimal_val_v<DecimalType> == 32, bool>
{
    using comp_type = std::uint_fast64_t;

    BOOST_DECIMAL_ASSERT(lhs_sig >= 0);
    BOOST_DECIMAL_ASSERT(rhs_sig >= 0);

    auto new_lhs_sig {static_cast<comp_type>(lhs_sig)};
    auto new_rhs_sig {static_cast<comp_type>(rhs_sig)};

    // Zeros are tested before the signs, since -0 is neither less than nor greater than +0
    if (new_lhs_sig == UINT64_C(0) || new_rhs_sig == UINT64_C(0))
    {
        return (new_lhs_sig == new_rhs_sig) ? false : (new_lhs_sig == 0U ? !rhs_sign : lhs_sign);
    }

    if (lhs_sign != rhs_sign)
    {
        return lhs_sign;
    }

    const auto delta_exp {lhs_exp - rhs_exp};
    constexpr auto max_delta_diff {std::numeric_limits<std::uint_fast64_t>::digits10 - detail::precision_v<DecimalType>};

    // If we can't do this correctly without normalization then do it and try again
    if (delta_exp > max_delta_diff || delta_exp < -max_delta_diff)
    {
        if (!normalized)
        {
            detail::normalize(lhs_sig, lhs_exp);
            detail::normalize(rhs_sig, rhs_exp);
            return less_parts_impl(lhs_sig, lhs_exp, lhs_sign,
                                   rhs_sig, rhs_exp, rhs_sign, true);
        }
        else
        {
            return rhs_sign ? rhs_exp < lhs_exp : rhs_exp > lhs_exp;
        }
    }

    if (delta_exp >= 0)
    {
        new_lhs_sig *= detail::pow10(static_cast<comp_type>(delta_exp));
        lhs_exp -= delta_exp;
    }
    else
    {
        new_rhs_sig *= detail::pow10(static_cast<comp_type>(-delta_exp));
        rhs_exp += delta_exp;
    }

    if (lhs_exp != rhs_exp)
    {
        return lhs_sign ? lhs_exp > rhs_exp : lhs_exp < rhs_exp;
    }

    return lhs_sign ? new_lhs_sig > new_rhs_sig : new_lhs_sig < new_rhs_sig;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType = decimal32_t, BOOST_DECIMAL_INTEGRAL T1,
        BOOST_DECIMAL_INTEGRAL U1, BOOST_DECIMAL_INTEGRAL T2, BOOST_DECIMAL_INTEGRAL U2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto less_parts_impl(T1 lhs_sig, U1 lhs_exp, bool lhs_sign,
                               T2 rhs_sig, U2 rhs_exp, bool rhs_sign) noexcept -> std::enable_if_t<detail::decimal_val_v<DecimalType> == 64 || detail::decimal_val_v<DecimalType> == 128, bool>
{
    using comp_type = typename DecimalType::significand_type;

    BOOST_DECIMAL_ASSERT(lhs_sig >= 0U);
    BOOST_DECIMAL_ASSERT(rhs_sig >= 0U);

    auto new_lhs_sig {static_cast<comp_type>(lhs_sig)};
    auto new_rhs_sig {static_cast<comp_type>(rhs_sig)};

    // Zeros are tested before the signs, since -0 is neither less than nor greater than +0
    if (new_lhs_sig == UINT64_C(0) || new_rhs_sig == UINT64_C(0))
    {
        return (new_lhs_sig == new_rhs_sig) ? false : (new_lhs_sig == 0U ? !rhs_sign : lhs_sign);
    }

    if (lhs_sign != rhs_sign)
    {
        return lhs_sign;
    }

    detail::normalize<DecimalType>(new_lhs_sig, lhs_exp);
    detail::normalize<DecimalType>(new_rhs_sig, rhs_exp);

    if (lhs_exp != rhs_exp)
    {
        return lhs_sign ? lhs_exp > rhs_exp : lhs_exp < rhs_exp;
    }

    return lhs_sign ? new_lhs_sig > new_rhs_sig : new_lhs_sig < new_rhs_sig;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType = decimal32_t, BOOST_DECIMAL_INTEGRAL T1,
          BOOST_DECIMAL_INTEGRAL U1, BOOST_DECIMAL_INTEGRAL T2, BOOST_DECIMAL_INTEGRAL U2>
constexpr auto less_parts_impl(T1 lhs_sig, U1 lhs_exp, bool lhs_sign,
                               T2 rhs_sig, U2 rhs_exp, bool rhs_sign) noexcept -> std::enable_if_t<detail::is_fast_type_v<DecimalType>, bool>
{
    using comp_type = std::conditional_t<(std::numeric_limits<T1>::digits10 > std::numeric_limits<T2>::digits10), T1, T2>;

    BOOST_DECIMAL_ASSERT(lhs_sig >= 0U);
    BOOST_DECIMAL_ASSERT(rhs_sig >= 0U);

    auto new_lhs_sig {static_cast<comp_type>(lhs_sig)};
    auto new_rhs_sig {static_cast<comp_type>(rhs_sig)};

    detail::normalize<DecimalType>(new_lhs_sig, lhs_exp);
    detail::normalize<DecimalType>(new_rhs_sig, rhs_exp);

    if (new_lhs_sig == 0U || new_rhs_sig == 0U)
    {
        if (new_lhs_sig == 0U && new_rhs_sig == 0U)
        {
            return false;
        }
        return new_lhs_sig == 0U ? !rhs_sign : lhs_sign;
    }

    if (lhs_sign != rhs_sign)
    {
        return lhs_sign;
    }

    if (lhs_exp != rhs_exp)
    {
        return lhs_sign ? lhs_exp > rhs_exp : lhs_exp < rhs_exp;
    }

    return lhs_sign ? new_lhs_sig > new_rhs_sig : new_lhs_sig < new_rhs_sig;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal, BOOST_DECIMAL_INTEGRAL Integer>
BOOST_DECIMAL_CUDA_CONSTEXPR auto less_impl(Decimal lhs, Integer rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal> && detail::is_integral_v<Integer>), bool>
{
    using exp_type = typename Decimal::biased_exponent_type;

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs))
    {
        return false;
    }
    else if (isinf(lhs))
    {
        return lhs.isneg();
    }
    #endif

    const auto lhs_significand {lhs.full_significand()};

    // A negative zero is not less than any integer that a positive zero is not less than,
    // so the sign shortcuts below have to stand aside when the decimal is a zero
    const bool lhs_is_zero {lhs_significand == 0U};

    bool lhs_sign {lhs.isneg()};
    bool rhs_sign {false};

    BOOST_DECIMAL_IF_CONSTEXPR (detail::is_signed_v<Integer>)
    {
        if (rhs < static_cast<Integer>(0))
        {
            rhs_sign = true;
        }

        if (!lhs_is_zero)
        {
            if (lhs_sign && !rhs_sign)
            {
                return true;
            }
            else if (!lhs_sign && rhs_sign)
            {
                return false;
            }
        }
    }
    else
    {
        if (lhs_sign && !lhs_is_zero)
        {
            return true;
        }
    }

    const auto rhs_significand {detail::make_positive_unsigned(rhs)};

    return less_parts_impl<Decimal>(lhs_significand, lhs.biased_exponent(), lhs_sign,
                                    rhs_significand, static_cast<exp_type>(0), rhs_sign);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto mixed_decimal_less_impl(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    using Bigger_Decimal_Type = std::conditional_t<(sizeof(lhs) > sizeof(rhs)), Decimal1, Decimal2>;

    // Two zeros are equal whatever their signs and cohorts are, so neither is less
    // (IEEE 754-2008 5.11). This precedes the sign shortcut below, which would otherwise
    // report -0 as less than +0.
    if (
            #ifndef BOOST_DECIMAL_FAST_MATH
            boost::decimal::isfinite(lhs) && boost::decimal::isfinite(rhs) &&
            #endif
            lhs.full_significand() == 0U && rhs.full_significand() == 0U
        )
    {
        return false;
    }

    if (
            #ifndef BOOST_DECIMAL_FAST_MATH
            isnan(lhs) || isnan(rhs) ||
            #endif
            (!lhs.isneg() && rhs.isneg())
        )
    {
        return false;
    }
    else if (lhs.isneg() && !rhs.isneg())
    {
        return true;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (isinf(lhs) || isinf(rhs))
    {
        // The signs are equal here, and the types decode an infinity as different values
        return isinf(lhs) ? lhs.isneg() && !isinf(rhs) : !rhs.isneg();
    }
    #endif

    return less_parts_impl<Bigger_Decimal_Type>(lhs.full_significand(), lhs.biased_exponent(), lhs.isneg(),
                                                rhs.full_significand(), rhs.biased_exponent(), rhs.isneg());
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator<(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    return mixed_decimal_less_impl(lhs, rhs);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator<=(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs) || isnan(rhs))
    {
        return false;
    }
    #endif

    return !(rhs < lhs);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator>(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    return rhs < lhs;
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator>=(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), bool>
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(lhs) || isnan(rhs))
    {
        return false;
    }
    #endif

    return !(lhs < rhs);
}

#ifdef BOOST_DECIMAL_HAS_SPACESHIP_OPERATOR

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal1, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal2>
BOOST_DECIMAL_CUDA_CONSTEXPR auto operator<=>(Decimal1 lhs, Decimal2 rhs) noexcept
    -> std::enable_if_t<(detail::is_decimal_floating_point_v<Decimal1> &&
                         detail::is_decimal_floating_point_v<Decimal2>), std::partial_ordering>
{
    if (lhs < rhs)
    {
        return std::partial_ordering::less;
    }
    else if (lhs > rhs)
    {
        return std::partial_ordering::greater;
    }
    else if (lhs == rhs)
    {
        return std::partial_ordering::equivalent;
    }

    return std::partial_ordering::unordered;
}

#endif

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_COMPARISON_HPP
