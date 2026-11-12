// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP
#define BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/add_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cerrno>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// True if the target type holds each coefficient of the source and needs no search for its digits.
// A fast type has no coefficient above its precision, so it also packs to the IEEE type of its size.
template <typename TargetType, typename Decimal>
using is_exact_repack = std::integral_constant<bool, (precision_v<TargetType> > precision_v<Decimal> &&
                                                      (is_fast_type_v<Decimal> || !is_fast_type_v<TargetType>)) ||
                                                     (precision_v<TargetType> == precision_v<Decimal> &&
                                                      is_fast_type_v<Decimal> && !is_fast_type_v<TargetType>)>;

// The value keeps its coefficient. Two fast types differ by a constant count of zeros,
// and a zero keeps the exponent that the constructor gives it.
template <typename TargetType, typename Decimal, typename SigType, typename ExpType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto finite_to_decimal(const SigType sig, const ExpType exp, const bool sign,
                                                    std::true_type) noexcept -> TargetType
{
    using target_sig = typename TargetType::significand_type;
    constexpr int zeros {is_fast_type_v<TargetType> ? precision_v<TargetType> - precision_v<Decimal> : 0};
    if (is_true(zeros != 0) && sig == 0U)
    {
        return TargetType{sig, exp, sign};
    }

    const auto shift {pow10(static_cast<target_sig>(zeros))};
    const auto wide_sig {static_cast<target_sig>(static_cast<target_sig>(sig) * shift)};
    return pack_in_range<TargetType>(wide_sig, static_cast<int>(exp) - zeros, sign);
}

template <typename TargetType, typename Decimal, typename SigType, typename ExpType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto finite_to_decimal(const SigType sig, const ExpType exp, const bool sign,
                                                    std::false_type) noexcept -> TargetType
{
    return TargetType{sig, exp, sign};
}

} // namespace detail

#ifndef BOOST_DECIMAL_NO_CXX17_IF_CONSTEXPR

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE TargetType, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_decimal(Decimal val) noexcept -> TargetType
{
    if (isinf(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::infinity() :
                              std::numeric_limits<TargetType>::infinity();
    }
    else if (issignaling(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::signaling_NaN() :
                              std::numeric_limits<TargetType>::signaling_NaN();
    }
    else if (isnan(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::quiet_NaN() :
                              std::numeric_limits<TargetType>::quiet_NaN();
    }

    BOOST_DECIMAL_IF_CONSTEXPR (std::is_same<TargetType, Decimal>::value)
    {
        return val;
    }
    else
    {
        return detail::finite_to_decimal<TargetType, Decimal>(
            val.full_significand(), val.biased_exponent(), val.isneg(), detail::is_exact_repack<TargetType, Decimal>{});
    }
}

#else

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE TargetType, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE Decimal>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_decimal(Decimal val) noexcept -> TargetType
{
    if (isinf(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::infinity() :
                              std::numeric_limits<TargetType>::infinity();
    }
    else if (issignaling(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::signaling_NaN() :
                              std::numeric_limits<TargetType>::signaling_NaN();
    }
    else if (isnan(val))
    {
        return val.isneg() ? -std::numeric_limits<TargetType>::quiet_NaN() :
                              std::numeric_limits<TargetType>::quiet_NaN();
    }

    return detail::finite_to_decimal<TargetType, Decimal>(
        val.full_significand(), val.biased_exponent(), val.isneg(), detail::is_exact_repack<TargetType, Decimal>{});
}

#endif

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_TO_DECIMAL_HPP
