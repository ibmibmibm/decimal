// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt


#ifndef BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP
#define BOOST_DECIMAL_DETAIL_TO_FLOAT_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/fast_float/compute_float32.hpp>
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

template <typename Decimal, typename TargetType>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_float(Decimal val) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, Decimal, detail::is_floating_point_v, TargetType, TargetType)
{
    auto fp_class = fpclassify(val);

    switch (fp_class)
    {
        case FP_NAN:
            if (issignaling(val))
            {
                return std::numeric_limits<TargetType>::signaling_NaN();
            }
            return std::numeric_limits<TargetType>::quiet_NaN();
        case FP_INFINITE:
            return val.isneg() ? -std::numeric_limits<TargetType>::infinity() : std::numeric_limits<TargetType>::infinity();
        case FP_ZERO:
            return val.isneg() ? -TargetType{0} : TargetType{0};
        default:
            break;
    }

    // The casts to result are redundant, but in pre C++17 modes MSVC warns about implicit conversions
    TargetType result {};

    auto sig {val.full_significand()};
    auto exp {val.biased_exponent()};

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
        if (wide_sig.high != 0U || !detail::fast_float::clinger_fast_path(exp, wide_sig.low, val.isneg(), value))
        {
            value = detail::fast_float::to_binary<binary_type>(detail::fast_float::compute_float_wide<binary_type>(exp, wide_sig), val.isneg());
        }
        result = static_cast<TargetType>(value);
    }
    else
    {
        bool success {};
        result = static_cast<TargetType>(detail::fast_float::compute_float80_128(exp, sig, val.isneg(), success));
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
