// Copyright 2024 Matt Borland
// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_RIEMANN_ZETA_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_RIEMANN_ZETA_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/impl/riemann_zeta_impl.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <cstdint>
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

template <typename T>
constexpr auto riemann_zeta_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    constexpr int p {precision_v<T>};
    using sig_type = boost::int128::uint128_t;

    const auto fpc = fpclassify(x);

    constexpr T one  { 1, 0 };

    const bool is_neg { signbit(x) };

    if (fpc == FP_ZERO)
    {
        // The value of riemann_zeta(0) is -1/2.

        return T { 5U, -1, true };
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (fpc != FP_NORMAL && fpc != FP_SUBNORMAL)
    {
        // zeta has no limit at -infinity, where it changes sign at each even integer
        return ((fpc == FP_INFINITE) ? (is_neg ? std::numeric_limits<T>::quiet_NaN() : one) : (issignaling(x) ? nan_conversion(x) : x));
    }
    #endif
    if (x == one)
    {
        // The argument is exactly one. The result is complex-infinity.
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T{0};
        #endif
    }

    int e {};
    const auto sig {static_cast<sig_type>(frexp10(x, &e))};

    // zeta(x) = -1/2 - x ln(2 pi) / 2 + O(x^2) for |x| < 10^-30, or 10^-20 and 10^-11 for fewer digits
    if (e + p <= (p <= 7 ? -11 : p <= 16 ? -20 : -30))
    {
        const auto c {fx_top<fx_words(p)>(fx_table::ln_sqrt_2pi)};
        const auto t {e + 37 >= 0 ? fx_scale(sig * pow10(static_cast<sig_type>(e + 37)), c)
                                  : (-e - 37 < 39 ? fx_scale(sig, c) / pow10(static_cast<sig_type>(-e - 37)) : sig_type {0U})};
        const auto half {pow10(static_cast<sig_type>(36)) * 5U};
        return fx_round<T>(is_neg ? half - t - 1U : half + t, -37, true);
    }

    // zeta(x) = 1 + 2^-x + ... rounds as a value just above 1 where 2^-x is below half a unit in the last place
    if (x > T {p <= 7 ? 21 : p <= 16 ? 51 : 111})
    {
        return fx_round<T>(pow10(static_cast<sig_type>(37)), -37, false);
    }
    return is_neg ? riemann_zeta_detail::zeta_negative<T>(sig, e) : riemann_zeta_detail::zeta_positive(x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto riemann_zeta(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::riemann_zeta_impl(static_cast<evaluation_type>(x)));
}

BOOST_DECIMAL_EXPORT template <typename T, typename IntegralType>
constexpr auto riemann_zeta(const IntegralType n) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, T, detail::is_integral_v, IntegralType, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    // TODO(ckormanyos) Consider making an integral-argument specialization.
    // Some exact values are known. Some simplifications for small-n are possible.

    return static_cast<T>(detail::riemann_zeta_impl(static_cast<evaluation_type>(n)));
}

} //namespace decimal
} //namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_RIEMANN_ZETA_HPP
