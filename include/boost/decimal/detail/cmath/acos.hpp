// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ACOS_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ACOS_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/asin.hpp>


#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <cstdint>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace atan_detail {

// acos of x not NaN, rounded once in the current mode
template <typename T>
constexpr auto acos_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // |x| >= 1 and infinity have e + p > 0: acos(1) = 0, acos(-1) = pi, and NaN above
    if (sig != 0U && e + p > 0)
    {
        if (asin_is_one(sig, e, p))
        {
            return neg ? pi_quarters<T>(4, false) : T {0};
        }
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    // For |x| < 1/256 and x = 0, acos(x) = pi/2 - asin(x)
    if (sig == 0U || decimal_tiny(sig, e))
    {
        const auto u {fx_from_decimal<words>(sig, -e)};
        const auto half_pi {fx_top<words>(fx_table::half_pi)};
        const auto s {fx_mul(u, asin_series(fx_mul(u, u)))};
        return fx_wide_round<T, words>(fx_widen(neg ? fx_add(half_pi, s) : fx_sub(half_pi, s)), false);
    }

    // For x = 1 - w with w < 2^-15, w is exact and acos(1 - w) = sqrt(2w) * asin_series(w / 2).
    // 2w = q * 10^(k + nd) with q in [0.01, 1), and an even k and nd.
    const auto one {pow10(static_cast<boost::int128::uint128_t>(-e))};
    const auto w {one - sig};
    if (!neg && e == -p && w * 32768U < one)
    {
        auto two_w {w * 2U};
        int k {e};
        if ((k & 1) != 0)
        {
            two_w *= 10U;
            --k;
        }
        const int nd {(num_digits(two_w) + 1) & ~1};
        const auto g {fx_mul(fx_sqrt(fx_from_decimal<words>(two_w, nd)), asin_series(fx_shr(fx_from_decimal<words>(w, -e), 1U)))};
        return fx_wide_round<T, words>(fx_widen(g), false, (k + nd) / 2);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // The same for x < 0.99, where the result is above 0.1: pi/2 and pi add 1 unit each
        T r {};
        if (neg || !decimal_at_least(sig, e, p, 99U, 0))
        {
            // acos(x) = pi - acos(|x|) for x < 0
            constexpr int n {words == 3 ? 2 : 1};
            const auto theta {asin_value<n>(sig, e, p, true)};
            if (fx_wide_try_round(fx_widen(neg ? fx_sub(fx_top<n>(fx_table::pi), theta) : theta), false, r, words == 3 ? 6U : 3U))
            {
                return r;
            }
        }
    }
    const auto theta {asin_value<words>(sig, e, p, true)};
    return fx_wide_round<T, words>(fx_widen(neg ? fx_sub(fx_top<words>(fx_table::pi), theta) : theta), false);
}

} // namespace atan_detail

template <typename T>
constexpr auto acos_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(x))
    {
        return issignaling(x) ? nan_conversion(x) : x;
    }
    #endif

    return atan_detail::acos_fixed(x);
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto acos(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::acos_impl(static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ACOS_HPP
