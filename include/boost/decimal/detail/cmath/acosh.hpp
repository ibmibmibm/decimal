// Copyright 2023 Matt Borland
// Copyright 2023 - 2026 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ACOSH_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ACOSH_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/asinh.hpp>
#include <boost/decimal/detail/concepts.hpp>


#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace log_detail {

// u + sqrt(u^2 - 10^(-2 n10)) for x = u * 10^n10 >= 1 and u in [0.1, 1)
template <int N>
constexpr auto acosh_v(const boost::int128::uint128_t& sig, const int p, const int n10) noexcept -> fx<N>
{
    const auto u {fx_from_decimal<N>(sig, p)};
    return fx_add(u, fx_sqrt(fx_sub(fx_mul(u, u), fx_from_decimal<N>(1U, 2 * n10))));
}

// acosh of a finite x != 0, rounded once in the current mode
template <typename T>
constexpr auto acosh_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};

    // x < 1 is negative or has e + p <= 0, and acosh(x) is -NaN
    if (signbit(x) || e + p <= 0)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return -std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }

    // For x = 1 + w with w < 2^-15, w is exact and acosh(1 + w) = sqrt(2w) * asinh_series(w / 2).
    // 2w = q * 10^(k + nd) with q in [0.01, 1), and an even k and nd.
    if (e == 1 - p)
    {
        const auto one {pow10(static_cast<boost::int128::uint128_t>(p - 1))};
        const auto w {sig - one};

        // A non-canonical x such as 10000001 * 10^-7 can give sig = 10^(p - 1), and acosh(1) = 0
        if (w == 0U)
        {
            return T {0};
        }
        if (w * 32768U < one)
        {
            auto two_w {w * 2U};
            int k {e};
            if ((k & 1) != 0)
            {
                two_w *= 10U;
                --k;
            }
            const int nd {(num_digits(two_w) + 1) & ~1};
            const auto g {fx_mul(fx_sqrt(fx_from_decimal<words>(two_w, nd)), asinh_series(fx_shr(fx_from_decimal<words>(w, -e), 1U)))};
            return fx_wide_round<T, words>(fx_widen(g), false, (k + nd) / 2);
        }
    }

    // acosh(x) = ln(x + sqrt(x^2 - 1)) with x = u * 10^n10 and u in [0.1, 1):
    // x + sqrt(x^2 - 1) = (u + sqrt(u^2 - 10^(-2 n10))) * 10^n10.
    const int n10 {e + p};
    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer for x >= 2, with an error below 256 units of its last bit
        T r {};
        if (decimal_at_least(sig, e, p, 20U, 1)
            && fx_wide_try_round(log_wide_of(acosh_v<words == 3 ? 2 : 1>(sig, p, n10), n10), false, r, 8U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(log_wide_of(acosh_v<words>(sig, p, n10), n10), false);
}

} // namespace log_detail

template <typename T>
constexpr auto acosh_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpc = fpclassify(x);

    T result { };

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (fpc != FP_NORMAL)
    {
        if ((fpc == FP_INFINITE) && (!signbit(x)))
        {
            result = x;
        }
        else if (fpc == FP_ZERO)
        {
            result = -std::numeric_limits<T>::quiet_NaN();
        }
        else
        {
            result = x;
        }
    }
    #else
    if (fpc == FP_ZERO)
    {
        result = T{0, 0};
    }
    #endif
    else
    {
        result = log_detail::acosh_fixed(x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto acosh(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::acosh_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_ACOSH_HPP
