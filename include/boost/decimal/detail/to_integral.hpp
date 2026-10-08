// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_TO_INTEGRAL_HPP
#define BOOST_DECIMAL_DETAIL_TO_INTEGRAL_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/apply_sign.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/cmath/decompose.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cerrno>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

// MSVC 14.1 warns of unary minus being applied to unsigned type from numeric_limits::min
// 14.2 and on get it right
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4146)
#endif

namespace detail {

// Multiplies value by power, and returns false when the product overflows
BOOST_DECIMAL_CUDA_CONSTEXPR auto checked_mul(std::uint64_t& value, const std::uint64_t power) noexcept -> bool
{
    std::uint64_t high {};
    value = int128::detail::umul(value, power, high);
    return high == 0U;
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto checked_mul(int128::uint128_t& value, const int128::uint128_t& power) noexcept -> bool
{
    const auto product {umul256(value, power)};
    value = int128::uint128_t {product.bytes[1], product.bytes[0]};
    return (product.bytes[2] | product.bytes[3]) == 0U;
}

// With BOOST_DECIMAL_DETAIL_INT128_NO_BUILTIN_INT128 the int128 classes have no conversion to a built-in
// 128-bit integer, so a magnitude for a type wider than 64 bits is made from both words
template <typename T, typename U, std::enable_if_t<(sizeof(T) > sizeof(std::uint64_t) && sizeof(U) > sizeof(std::uint64_t)), bool> = true>
BOOST_DECIMAL_CUDA_CONSTEXPR auto cast_magnitude(const U& value) noexcept -> T
{
    return (static_cast<T>(value.high) << 64U) | static_cast<T>(value.low);
}

template <typename T, typename U, std::enable_if_t<!(sizeof(T) > sizeof(std::uint64_t) && sizeof(U) > sizeof(std::uint64_t)), bool> = true>
BOOST_DECIMAL_CUDA_CONSTEXPR auto cast_magnitude(const U& value) noexcept -> T
{
    return static_cast<T>(value);
}

} // namespace detail

// Attempts conversion to integral type:
// If this is nan sets errno to EINVAL and returns the largest value of the type
// If this is not representable sets errno to ERANGE and returns the largest value of the type
template <typename Decimal, typename TargetType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_integral(Decimal val) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, Decimal, detail::is_integral_v, TargetType, TargetType)
{
    using Conversion_Type = std::conditional_t<(sizeof(TargetType) > sizeof(std::uint64_t) || sizeof(decltype(decompose(val).sig)) > sizeof(std::uint64_t)),
                                               int128::uint128_t, std::uint64_t>;

    // std::numeric_limits has no __int128 in a strict mode of an old libstdc++, so the size gives the digits there
    constexpr int target_digits {std::numeric_limits<TargetType>::is_specialized ? std::numeric_limits<TargetType>::digits :
                                 static_cast<int>(sizeof(TargetType) * 8U) - (detail::is_signed_v<TargetType> ? 1 : 0)};
    static_assert(target_digits <= std::numeric_limits<Conversion_Type>::digits, "The target type has more than 128 bits");
    constexpr int unused_bits {std::numeric_limits<Conversion_Type>::digits - target_digits};
    constexpr Conversion_Type max_magnitude {(std::numeric_limits<Conversion_Type>::max)() >> unused_bits};

    if (isnan(val))
    {
        #if defined(__clang__) && __clang_major__ >= 20
        if (!BOOST_DECIMAL_IS_CONSTANT_EVALUATED(val))
        #endif
        {
            errno = EINVAL;
        }
        return detail::cast_magnitude<TargetType>(max_magnitude);
    }

    // The range check is on the truncated magnitude, because the limits can round in the decimal type.
    // A negative value has the limit -min, which is 0 for an unsigned type.
    const auto components {decompose(val)};
    const auto limit {!components.sign ? max_magnitude :
                      detail::is_signed_v<TargetType> ? static_cast<Conversion_Type>(max_magnitude + 1U) : Conversion_Type {}};

    auto result {static_cast<Conversion_Type>(components.sig)};
    auto expval {components.exp};
    const auto abs_exp_val {detail::make_positive_unsigned(expval)};
    bool in_range {!isinf(val)};

    if (expval > 0 && result != 0U)
    {
        const auto power {expval > std::numeric_limits<Conversion_Type>::digits10 ? Conversion_Type {} :
                                                                                    detail::pow10<Conversion_Type>(static_cast<Conversion_Type>(expval))};
        in_range = in_range && power != 0U && detail::checked_mul(result, power);
    }
    else if (expval < 0)
    {
        result = abs_exp_val > std::numeric_limits<Conversion_Type>::digits10 ? Conversion_Type {} :
                                                                                result / detail::pow10<Conversion_Type>(static_cast<Conversion_Type>(abs_exp_val));
    }

    if (!in_range || result > limit)
    {
        #if defined(__clang__) && __clang_major__ >= 20
        if (!BOOST_DECIMAL_IS_CONSTANT_EVALUATED(val))
        #endif
        {
            errno = ERANGE;
        }

        return detail::cast_magnitude<TargetType>(max_magnitude);
    }

    if (!components.sign || result == 0U)
    {
        return detail::cast_magnitude<TargetType>(result);
    }

    // -result is -(result - 1) - 1, and result - 1 fits in the signed type. An unsigned type has result 0 here.
    using signed_type = std::conditional_t<detail::is_signed_v<TargetType>, TargetType, Conversion_Type>;
    const auto below {detail::cast_magnitude<signed_type>(static_cast<Conversion_Type>(result - 1U))};
    return static_cast<TargetType>(-below - 1);
}

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_TO_INTEGRAL_HPP
