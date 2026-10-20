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
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/impl/compare_long_double.hpp>
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
BOOST_DECIMAL_FORCE_INLINE constexpr auto next_value(SigType sig, int exp, const bool is_neg, const bool grow) noexcept -> DecimalType
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

#ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

BOOST_DECIMAL_EXPORT template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto nexttoward(const T val, const long double direction) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // direction can round to val in T, so the two compare exactly as numbers. Its bits also give a NaN or an
    // infinity, as an x87 compare loads it again.
    const auto parts {detail::split_long_double(direction)};

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(val))
    {
        return val;
    }
    else if (parts.nan)
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    else if (isinf(val))
    {
        if (parts.infinite && parts.negative == signbit(val))
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
    if (parts.infinite)
    {
        order = parts.negative ? -1 : 1;
    }
    else
    #endif
    {
        order = detail::compare_long_double<T>(parts, sig, exp, is_neg);

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

#endif // BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_NEXT_HPP
