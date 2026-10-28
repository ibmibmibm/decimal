// Copyright 2025 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_WRITE_PAYLOAD_HPP
#define BOOST_DECIMAL_DETAIL_WRITE_PAYLOAD_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/power_tables.hpp>

namespace boost {
namespace decimal {

BOOST_DECIMAL_CUDA_CONSTEXPR auto from_bits(const std::uint32_t rhs) noexcept -> decimal32_t;
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_bits(const std::uint64_t rhs) noexcept -> decimal64_t;
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_bits(const int128::uint128_t rhs) noexcept -> decimal128_t;

namespace detail {

// The largest canonical payload of an IEEE type is 10^(p - 1) - 1 (IEEE 754 3.5.2)
template <typename T>
BOOST_DECIMAL_CUDA_CONSTEXPR auto max_ieee_payload() noexcept -> typename T::significand_type
{
    return static_cast<typename T::significand_type>(pow10(static_cast<typename T::significand_type>(precision_v<T> - 1)) - 1U);
}

template <typename TargetDecimalType, bool is_snan>
BOOST_DECIMAL_CUDA_CONSTEXPR auto write_payload(typename TargetDecimalType::significand_type payload_value)
    BOOST_DECIMAL_REQUIRES(detail::is_fast_type_v, TargetDecimalType)
{
    using sig_type = typename TargetDecimalType::significand_type;

    constexpr TargetDecimalType nan_type {is_snan ? std::numeric_limits<TargetDecimalType>::signaling_NaN() :
                                                    std::numeric_limits<TargetDecimalType>::quiet_NaN()};

    constexpr std::uint32_t significand_field_bits {decimal_val_v<TargetDecimalType> < 64 ? 23U :
                                                    decimal_val_v<TargetDecimalType> < 128 ? 53U : 110U};

    constexpr sig_type max_payload_value {(static_cast<sig_type>(1) << significand_field_bits) - 1U};

    TargetDecimalType return_value {nan_type};
    if (payload_value < max_payload_value)
    {
        return_value.significand_ |= payload_value;
    }

    return return_value;
}

template <typename TargetDecimalType, bool is_snan>
BOOST_DECIMAL_CUDA_CONSTEXPR auto write_payload(typename TargetDecimalType::significand_type payload_value)
    BOOST_DECIMAL_REQUIRES(detail::is_ieee_type_v, TargetDecimalType)
{
    constexpr TargetDecimalType nan_type {is_snan ? std::numeric_limits<TargetDecimalType>::signaling_NaN() :
                                                    std::numeric_limits<TargetDecimalType>::quiet_NaN()};

    auto return_value {nan_type.bits_};
    if (payload_value <= max_ieee_payload<TargetDecimalType>())
    {
        return_value = payload_value | return_value;
    }

    return from_bits(return_value);
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto read_payload(const T value) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_ieee_type_v, T, typename T::significand_type)
{
    using sig_type = typename T::significand_type;

    // The payload is the trailing significand field. A larger value than the canonical limit gives 0.
    constexpr sig_type field_mask {(static_cast<sig_type>(1) << detail::trailing_significand_field_width_v<T>) - 1U};
    const auto payload {static_cast<sig_type>(value.bits_ & field_mask)};

    return isnan(value) && payload <= detail::max_ieee_payload<T>() ? payload : static_cast<sig_type>(0U);
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto read_payload(const T value) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_fast_type_v, T, typename T::significand_type)
{
    if (!isnan(value))
    {
        return 0U;
    }
    else if (issignaling(value))
    {
        return value.significand_ ^ std::numeric_limits<T>::signaling_NaN().significand_;
    }
    else
    {
        return value.significand_ ^ std::numeric_limits<T>::quiet_NaN().significand_;
    }
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_WRITE_PAYLOAD_HPP
