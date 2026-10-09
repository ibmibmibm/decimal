// Copyright 2023 Matt Borland
// Copyright 2023 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_LOG1P_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_LOG1P_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace log_detail {

// ln(1 + x) for |x| >= 1/256 and x > -1: 1 + x is exact with the exponent e for e < 0, and with the exponent 0
// where it has at most 38 digits. Else ln(1 + x) = ln(x) + 1/x - ... with 1/x = (10^38 / sig) * 10^-(38 + e) below 10^-38.
template <int N>
constexpr auto log1p_wide(const boost::int128::uint128_t& sig, const int e, const bool neg, const int p) noexcept -> fx_wide<N>
{
    if (e < 0)
    {
        const auto one {pow10(static_cast<boost::int128::uint128_t>(-e))};
        const auto y {neg ? one - sig : one + sig};
        return log_wide_of<N>(y, num_digits(y), e, log_base::e);
    }
    if (e + p <= 38)
    {
        const auto y {sig * pow10(static_cast<boost::int128::uint128_t>(e)) + 1U};
        return log_wide_of<N>(y, num_digits(y), 0, log_base::e);
    }
    const auto inv {fx_from_decimal<N>(pow10(static_cast<boost::int128::uint128_t>(38)) / sig, 38 + e)};
    return fx_wide_add(log_wide_of<N>(sig, p, e, log_base::e), fx_widen(inv));
}

// log1p of a finite x != 0, rounded once in the current mode
template <typename T>
constexpr auto log1p_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // x <= -1 has e + p > 0: log1p(-1) = -infinity, and NaN below, also for a non-canonical x such as
    // -10000001 * 10^-7 that frexp10 rounds to -1
    if (neg && e + p > 0)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return e + p == 1 && sig == pow10(static_cast<boost::int128::uint128_t>(p - 1)) && x == -T {1} ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    if (decimal_tiny(sig, e))
    {
        return log_near<T>(sig, e, neg, log_base::e);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer for |x| >= 0.1, and in one word for x >= 0.11 or x <= -0.1, with an
        // error below 64 units of its last bit, or 16 units of 2^-62 in one word
        T r {};
        if ((words == 3 ? e + p >= 0 : decimal_at_least(sig, e, p, neg ? 10U : 11U, 0))
            && fx_wide_try_round(log1p_wide<words == 3 ? 2 : 1>(sig, e, neg, p), false, r, words == 3 ? 6U : 4U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(log1p_wide<words>(sig, e, neg, p), false);
}

} // namespace log_detail

template <typename T>
constexpr auto log1p_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    T result { };

    const auto fpc = fpclassify(x);

    if (fpc == FP_ZERO)
    {
        result = x;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpc != FP_NORMAL && fpc != FP_SUBNORMAL)
    {
        result =
        (
            ((fpc == FP_INFINITE) && signbit(x)) ? std::numeric_limits<T>::quiet_NaN() : x
        );
    }
    #endif
    else
    {
        result = log_detail::log1p_fixed(x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto log1p(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::log1p_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_LOG1P_HPP
