//  (C) Copyright John Maddock 2008 - 2023.
//  (C) Copyright Matt Borland 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_DECIMAL_DETAIL_CMATH_NEXT_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_NEXT_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/add_impl.hpp>
#include <boost/decimal/detail/cmath/modf.hpp>
#include <boost/decimal/detail/cmath/abs.hpp>
#include <boost/decimal/detail/cmath/round.hpp>
#include <boost/decimal/detail/cmath/ilogb.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/fast_float/compute_float64.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/utilities.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The neighbor of sig * 10^exp, where frexp10 gave sig with all p digits. For a subnormal
// value, exp can be below etiny, so remove the added zeros to step by one unit of 10^etiny.
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType, typename SigType>
constexpr auto next_value(SigType sig, int exp, const bool is_neg, const bool grow) noexcept -> DecimalType
{
    if (exp < detail::etiny_v<DecimalType>)
    {
        sig /= detail::pow10(static_cast<SigType>(detail::etiny_v<DecimalType> - exp));
        exp = detail::etiny_v<DecimalType>;
    }

    constexpr auto min_full_sig {static_cast<SigType>(detail::max_significand_v<DecimalType> / 10U + 1U)};

    if (grow)
    {
        ++sig;
        if (sig > detail::max_significand_v<DecimalType>)
        {
            sig = min_full_sig;
            ++exp;
        }
    }
    else if (sig == min_full_sig && exp > detail::etiny_v<DecimalType>)
    {
        sig = detail::max_significand_v<DecimalType>;
        --exp;
    }
    else
    {
        // The smallest subnormal value steps to a zero with the exponent etiny
        --sig;
    }

    return detail::pack_in_range<DecimalType>(sig, exp, is_neg);
}

// frexp10 gives each value as a significand with all p digits and an exponent, so the
// two magnitudes compare as (exponent, significand) pairs. Zero has the significand 0.
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
constexpr auto nextafter_impl(const DecimalType val, const DecimalType direction) noexcept -> DecimalType
{
    const bool is_neg {signbit(val)};

    int exp {};
    const auto sig {frexp10(val, &exp)};
    int dir_exp {};
    const auto dir_sig {frexp10(direction, &dir_exp)};

    if (sig == dir_sig && exp == dir_exp && (is_neg == signbit(direction) || sig == 0U))
    {
        return direction;
    }

    if (sig == 0U)
    {
        return signbit(direction) ? -std::numeric_limits<DecimalType>::denorm_min() :
                                    std::numeric_limits<DecimalType>::denorm_min();
    }

    #ifndef BOOST_DECIMAL_FAST_MATH
    const bool dir_inf {isinf(direction)};
    #else
    constexpr bool dir_inf {false};
    #endif

    // The magnitude grows only when direction is on the same side of zero and farther out
    const bool grow {is_neg == signbit(direction) &&
                     (dir_inf || (dir_sig != 0U && (exp != dir_exp ? exp < dir_exp : sig < dir_sig)))};

    return next_value<DecimalType>(sig, exp, is_neg, grow);
}

// A finite long double as mantissa * 2^exponent
struct long_double_parts
{
    int128::uint128_t mantissa;
    int exponent;
    bool negative;
};

BOOST_DECIMAL_CXX20_CONSTEXPR auto split_long_double(const long double val) noexcept -> long_double_parts
{
    #if BOOST_DECIMAL_LDBL_BITS == 64
    const int128::uint128_t bits {bit_cast<std::uint64_t>(val)};
    constexpr int fraction_bits {52};
    constexpr int exponent_bits {11};
    constexpr bool explicit_bit {false};
    #elif BOOST_DECIMAL_LDBL_BITS == 80
    const auto bits {bit_cast<int128::uint128_t>(val)};
    constexpr int fraction_bits {64};
    constexpr int exponent_bits {15};
    constexpr bool explicit_bit {true};
    #else
    const auto bits {bit_cast<int128::uint128_t>(val)};
    constexpr int fraction_bits {112};
    constexpr int exponent_bits {15};
    constexpr bool explicit_bit {false};
    #endif

    constexpr int bias {(1 << (exponent_bits - 1)) - 1};
    const int128::uint128_t one {1U};
    const auto fraction {bits & ((one << fraction_bits) - 1U)};
    const auto biased {static_cast<int>((bits >> fraction_bits) & ((one << exponent_bits) - 1U))};
    const bool negative {((bits >> (fraction_bits + exponent_bits)) & 1U) != 0U};

    // The 80-bit format keeps the leading bit, and the other formats imply it for a normal value
    const auto mantissa {explicit_bit || biased == 0 ? fraction : fraction | (one << fraction_bits)};
    const auto exponent {(biased == 0 ? 1 : biased) - bias - (explicit_bit ? fraction_bits - 1 : fraction_bits)};

    return {mantissa, exponent, negative};
}

// The sign of |y| - sig * 10^exp, from big integers of Capacity words
template <std::size_t Capacity>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compare_magnitude_exact(const long_double_parts& y, const int128::uint128_t sig, const int exp) noexcept -> int
{
    fast_float::compare_bigint<Capacity> value {sig};
    fast_float::compare_bigint<Capacity> limit {y.mantissa};
    if (exp >= 0)
    {
        value.mul_pow5(exp);
    }
    else
    {
        limit.mul_pow5(-exp);
    }

    const std::int64_t two_shift {static_cast<std::int64_t>(exp) - y.exponent};
    if (two_shift >= 0)
    {
        value.shl(two_shift);
    }
    else
    {
        limit.shl(-two_shift);
    }

    return limit.compare(value);
}

// The sign of y - sig * 10^exp, from an exact compare. The magnitudes compare first by their
// bit lengths, and only values within a few binades of each other need the big integers.
template <typename SigType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compare_long_double(const long_double_parts& y, const SigType sig, const int exp, const bool is_neg) noexcept -> int
{
    if (y.mantissa == 0U || sig == 0U)
    {
        if (y.mantissa == 0U && sig == 0U)
        {
            return 0;
        }
        return y.mantissa != 0U ? (y.negative ? -1 : 1) : (is_neg ? 1 : -1);
    }

    if (y.negative != is_neg)
    {
        return y.negative ? -1 : 1;
    }

    const int128::uint128_t wide_sig {sig};
    const int y_bits {y.exponent + 128 - countl_zero(y.mantissa)};

    // 1701 / 512 is log2(10) to within 1e-3, so the estimate of the bit length of the
    // decimal value is off by less than 8 bits over the long double range
    const int x_bits {128 - countl_zero(wide_sig) + exp * 1701 / 512};

    int order {};
    if (y_bits > x_bits + 8)
    {
        order = 1;
    }
    else if (y_bits < x_bits - 8)
    {
        order = -1;
    }
    else if (exp >= -550 && exp <= 550)
    {
        // Both values are near 2^y_bits, so each big integer has fewer than 1450 bits
        order = compare_magnitude_exact<24>(y, wide_sig, exp);
    }
    else
    {
        // |exp| is at most about 5000 in the long double range
        order = compare_magnitude_exact<200>(y, wide_sig, exp);
    }

    return is_neg ? -order : order;
}

} // namespace detail

BOOST_DECIMAL_EXPORT
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T1,
          BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T2>
constexpr auto nextafter(const T1 val, const T2 direction) noexcept
    BOOST_DECIMAL_REQUIRES_TWO(detail::is_decimal_floating_point_v, T1, detail::is_decimal_floating_point_v, T2)
{
    // The step is one unit in the common type, as for the other functions of two arguments
    using ret_type = detail::promote_args_t<T1, T2>;

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(val))
    {
        return static_cast<ret_type>(val);
    }
    else if (isnan(direction))
    {
        return static_cast<ret_type>(direction);
    }
    else if (isinf(val))
    {
        if (val == direction)
        {
            return static_cast<ret_type>(direction);
        }

        // Stepping inward from an infinity lands on the largest finite value
        return signbit(val) ? -(std::numeric_limits<ret_type>::max)()
                            : (std::numeric_limits<ret_type>::max)();
    }
    #endif

    return detail::nextafter_impl(static_cast<ret_type>(val), static_cast<ret_type>(direction));
}

BOOST_DECIMAL_EXPORT template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto nexttoward(const T val, const long double direction) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // direction can round to val in T, so the two compare exactly as numbers
    #ifndef BOOST_DECIMAL_FAST_MATH
    constexpr long double inf {detail::infinity_value<long double>()};

    if (isnan(val))
    {
        return val;
    }
    else if (detail::is_nan_value(direction))
    {
        return static_cast<T>(direction);
    }
    else if (isinf(val))
    {
        if (signbit(val) ? !(direction > -inf) : !(direction < inf))
        {
            return val;
        }

        // Stepping inward from an infinity lands on the largest finite value
        return signbit(val) ? -(std::numeric_limits<T>::max)() : (std::numeric_limits<T>::max)();
    }
    #endif

    const bool is_neg {signbit(val)};
    int exp {};
    const auto sig {frexp10(val, &exp)};

    int order {};
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!(direction < inf) || !(direction > -inf))
    {
        order = direction > 0 ? 1 : -1;
    }
    else
    #endif
    {
        const auto parts {detail::split_long_double(direction)};
        order = detail::compare_long_double(parts, sig, exp, is_neg);

        // Equal values differ only in the sign of a zero, and the result has the sign of direction
        if (order == 0)
        {
            return parts.negative != is_neg ? -val : val;
        }
    }

    if (sig == 0U)
    {
        return order > 0 ? std::numeric_limits<T>::denorm_min() : -std::numeric_limits<T>::denorm_min();
    }

    return detail::next_value<T>(sig, exp, is_neg, (order > 0) != is_neg);
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_NEXT_HPP
