// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_FREXP10_HPP
#define BOOST_DECIMAL_DETAIL_FREXP10_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>
#include <boost/decimal/detail/countl.hpp>
#include "../int128.hpp"
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The encoding keeps a non-canonical significand at p + 1 digits and below 10^(p + 1) - 5,
// so fenv_round removes one digit. The helper also writes the exponent, so frexp10 tail calls it.
template <typename TargetType, typename T1, typename T2>
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto frexp10_round_non_canonical(T1 sig, T2 exp, int* expptr) noexcept -> T1
{
    static_cast<void>(fenv_round<TargetType>(sig));
    *expptr = static_cast<int>(exp) + 1;
    return sig;
}

// Hopefully the compiler makes this NOOP for the fast case
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE TargetType, BOOST_DECIMAL_INTEGRAL T1, BOOST_DECIMAL_INTEGRAL T2>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto frexp10_normalize(T1&, T2&) noexcept -> std::enable_if_t<is_fast_type_v<TargetType>, bool> { return false; }

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE TargetType, BOOST_DECIMAL_INTEGRAL T1, BOOST_DECIMAL_INTEGRAL T2>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto frexp10_normalize(T1& sig, T2& exp) noexcept -> std::enable_if_t<!is_fast_type_v<TargetType>, bool>
{
    // A significand with p digits is already normal. Returns true only for a non-canonical one,
    // which still needs rounding.
    constexpr int p {precision_v<TargetType>};
    constexpr auto ten_pm1 {pow10(static_cast<T1>(p - 1))};
    constexpr auto ten_p {pow10(static_cast<T1>(p))};
    if (sig < ten_pm1)
    {
        int t {};
        BOOST_DECIMAL_IF_CONSTEXPR (sizeof(T1) <= 4U)
        {
            // Count the digits with compares, because x86 bsr waits for the old value of its output
            for (int k {1}; k < p - 1; ++k)
            {
                t += static_cast<int>(sig >= pow10(static_cast<T1>(k)));
            }
        }
        else BOOST_DECIMAL_IF_CONSTEXPR (sizeof(T1) <= 8U)
        {
            // 1233 / 4096 is close to log10(2), so t is floor(log10(sig)) or one more. A factor of
            // ten fixes a short product.
            const int bits {static_cast<int>(sizeof(T1) * 8U) - countl_zero(sig)};
            t = (bits * 1233) >> 12;
            sig *= pow10(static_cast<T1>(p - 1 - t));
            exp -= p - 1 - t;
            // Add 9 * sig through a mask. Clang makes a branch of an if, which mispredicts
            // when the digit counts are mixed.
            const bool short_product {sig < ten_pm1};
            sig += (sig * 9U) & (T1{0U} - static_cast<T1>(short_product));
            exp -= static_cast<T2>(short_product);
            return false;
        }
        else
        {
            // floor(log10(sig)) is t or t - 1, and one compare picks it without a branch
            const int bits {static_cast<int>(sizeof(T1) * 8U) - countl_zero(sig)};
            t = (bits * 1233) >> 12;
            t -= static_cast<int>(sig < pow10(static_cast<T1>(t)));
        }
        const int zeros_needed {p - 1 - t};
        sig *= pow10(static_cast<T1>(zeros_needed));
        exp -= zeros_needed;
        return false;
    }
    return sig >= ten_p;
}

} // namespace detail

// Returns the normalized significand and exponent to be cohort agnostic
// Returns num in the range
//   [1e06, 1e06 - 1] for decimal32_t
//   [1e15, 1e15 - 1] for decimal64_t
// If the conversion can not be performed returns UINT32_MAX and exp = 0
BOOST_DECIMAL_EXPORT template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
BOOST_DECIMAL_CUDA_CONSTEXPR auto frexp10(const T num, int* expptr) noexcept -> typename T::significand_type
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(num))
    {
        *expptr = 0;
        return (std::numeric_limits<typename T::significand_type>::max)();
    }
    #endif

    auto num_exp {num.biased_exponent()};
    auto num_sig {num.full_significand()};

    // Normalize the handling of zeros
    if (num_sig == 0U)
    {
        *expptr = 0;
        return 0U;
    }

    if (detail::frexp10_normalize<T>(num_sig, num_exp))
    {
        return detail::frexp10_round_non_canonical<T>(num_sig, num_exp, expptr);
    }

    *expptr = static_cast<int>(num_exp);

    return num_sig;
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_FREXP10_HPP
