// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP
#define BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/write_payload.hpp>
#include <boost/decimal/detail/cmath/decompose.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cerrno>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The payload stays if the target holds it (IEEE 754 6.2.3), and only a change
// of precision is an operation, which gives a quiet NaN (IEEE 754 6.2)
template <typename TargetType, bool keep_signaling = false, typename Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto nan_to_decimal(const Decimal val) noexcept -> TargetType
{
    using target_sig = typename TargetType::significand_type;
    const auto payload {read_payload(val)};
    const auto target_payload {static_cast<target_sig>(payload) == payload ? static_cast<target_sig>(payload) : target_sig {0}};
    const auto nan_val {(keep_signaling || std::numeric_limits<TargetType>::digits10 == std::numeric_limits<Decimal>::digits10) && issignaling(val) ?
                        detail::write_payload<TargetType, true>(target_payload) :
                        detail::write_payload<TargetType, false>(target_payload)};

    return signbit(val) ? -nan_val : nan_val;
}

template <typename TargetType, bool keep_signaling, typename Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto convert_decimal(const Decimal val, std::false_type) noexcept -> TargetType
{
    if (isinf(val))
    {
        return signbit(val) ? -std::numeric_limits<TargetType>::infinity() :
                               std::numeric_limits<TargetType>::infinity();
    }
    else if (BOOST_DECIMAL_UNLIKELY(isnan(val)))
    {
        return nan_to_decimal<TargetType, keep_signaling>(val);
    }

    const auto c {decompose(val)};
    return TargetType{c.sig, c.exp, c.sign};
}

// A conversion to the same type returns the value unchanged
template <typename TargetType, bool keep_signaling, typename Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto convert_decimal(const Decimal val, std::true_type) noexcept -> TargetType
{
    return val;
}

// Mixed arithmetic keeps a signaling NaN operand, so the operation selects the NaN as for one type
template <typename TargetType, typename Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto promote_operand(const Decimal val) noexcept -> TargetType
{
    return convert_decimal<TargetType, true>(val, std::is_same<TargetType, Decimal> {});
}

} // namespace detail

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE TargetType, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_decimal(Decimal val) noexcept -> TargetType
{
    return detail::convert_decimal<TargetType, false>(val, std::is_same<TargetType, Decimal> {});
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP
