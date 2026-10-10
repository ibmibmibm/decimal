// Copyright 2024 - 2025 Matt Borland
// Copyright 2024 - 2025 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ASIN_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ASIN_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/fabs.hpp>
#include <boost/decimal/detail/cmath/atan.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace atan_detail {

// asin(x) / x = 1 + t/6 + 3t^2/40 + ... for t = x^2 < 2^-16, with the coefficients of asinh_series
template <int N>
constexpr auto asin_series(const fx<N>& t) noexcept -> fx<N>
{
    return fx_series<false, N == 1 ? 4U : N == 2 ? 8U : 12U, 4U, 8U>(t, fx_table::asinh_coeff);
}

// sqrt(1 - x^2) for x = sig * 10^e with 1/256 <= |x| < 1. 1 - x^2 = a * b * 10^(2e) with a = 1 + |x| and
// b = 1 - |x| exact, and = ua * ub * 10^n with ua and ub in [0.1, 1), which keeps the digits of b.
template <int N, bool Short = false>
constexpr auto asin_cos(const boost::int128::uint128_t& sig, const int e) noexcept -> fx<N>
{
    // One word with e >= -9, as 10^18 fits in 64 bits: 1 - x^2 = d * 10^(2e) with d = 10^(-2e) - sig^2 exact, and d / 10^n in [0.01, 1) for an even n
    BOOST_DECIMAL_IF_CONSTEXPR (N == 1)
    {
        if (e >= -9)
        {
            const auto d {pow10(static_cast<std::uint64_t>(-2 * e)) - sig.low * sig.low};
            const int n {(num_digits(d) + 1) & ~1};
            const auto r {fx_sqrt<N, Short ? 1 : 2>(fx_from_decimal<N>(d, n))};
            return n == -2 * e ? r : fx_mul(r, fx_from_decimal<N>(1U, -e - n / 2));
        }
    }
    // The same with e >= -19 in 128 bits, and d * 4^s / 10^n in [0.25, 1), so that the sqrt is of full precision
    if (e >= -19)
    {
        const auto d {pow10(static_cast<boost::int128::uint128_t>(-2 * e)) - sig * sig};
        const int n {(num_digits(d) + 1) & ~1};
        const auto top {pow10(static_cast<boost::int128::uint128_t>(n))};
        const unsigned s {static_cast<unsigned>(d * 4U < top) + static_cast<unsigned>(d * 16U < top) + static_cast<unsigned>(d * 64U < top)};
        const auto r {fx_sqrt(fx_from_decimal<N>(d << (2U * s), n))};
        return fx_shr(n == -2 * e ? r : fx_mul(r, fx_from_decimal<N>(1U, -e - n / 2)), s);
    }
    const auto one {pow10(static_cast<boost::int128::uint128_t>(-e))};
    const auto a {one + sig};
    const auto b {one - sig};
    const int na {num_digits(a)};
    const int nb {num_digits(b)};
    auto m {fx_mul(fx_from_decimal<N>(a, na), fx_from_decimal<N>(b, nb))};
    int n {na + nb + 2 * e};

    // a has at most 1 - e digits and b at most -e, thus n <= 1, and m * 10^n <= 1; then an even n <= 0
    if (n == 1)
    {
        m = fx_mul_small(m, 10U);
        n = 0;
    }
    if ((n & 1) != 0)
    {
        m = fx_mul(m, fx_from_decimal<N>(1U, 1));
        ++n;
    }
    const auto r {fx_sqrt<N, Short ? 1 : 2>(m)};
    return n == 0 ? r : fx_mul(r, fx_from_decimal<N>(1U, -n / 2));
}

// asin(r) / r = 1 + t/6 + 3t^2/40 + ... for t = r^2 <= 2^-12, with the term counts of atan_series, or 3 with Short
template <int N, bool Short = false>
constexpr auto asin_series_long(const fx<N>& t) noexcept -> fx<N>
{
    return fx_series<false, Short ? 3U : N == 1 ? 6U : N == 2 ? 11U : 16U, 6U, 11U>(t, fx_table::asinh_coeff);
}

// asin(z) for 0 <= z, w = sqrt(1 - z^2) and z / w < 1 + 1/64. With c = i/32 next to z / w and
// t = atan(c), r = sin(asin(z) - t) = cos(t) * (z - c w), |r| <= 1/64, and asin(z) = t + asin(r).
template <bool Short = false, int N>
constexpr auto asin_of(const fx<N>& z, const fx<N>& w) noexcept -> fx<N>
{
    const auto i {atan_index(z, w)};
    const auto cw {fx_mul_shr(w, i, 5U)};
    const bool below {fx_less(z, cw)};
    const auto r {fx_mul(below ? fx_sub(cw, z) : fx_sub(z, cw), fx_top<N>(fx_table::atan_cos[i]))};
    const auto s {fx_mul(r, asin_series_long<N, Short>(fx_mul(r, r)))};
    const auto c {fx_top<N>(fx_table::atan_c[i])};
    return below ? fx_sub(c, s) : fx_add(c, s);
}

// |x| < 0.71 for x = sig * 10^e with p digits, and then sqrt(1 - x^2) > 0.704
constexpr auto asin_is_small(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> bool
{
    return !decimal_at_least(sig, e, p, 71U, 0);
}

// asin(u) for 0 <= u < 1/2 from the Taylor series at the center y0 of one of eight intervals, in t = 32 |u - y0|.
// The terms are all positive above y0 and alternate below.
constexpr auto asin_poly(const fx<1>& u) noexcept -> fx<1>
{
    const auto i {static_cast<int>(u.w[0] >> 58U)};
    const fx<1> y0 {{static_cast<std::uint64_t>(2 * i + 1) << 57U}};
    const bool below {fx_less(u, y0)};
    const auto t {fx_shl(below ? fx_sub(y0, u) : fx_sub(u, y0), 5U)};
    const auto a {fx_make_span<9>(fx_table::asin_d32, i * 9)};
    fx<1> s {{a[a.size() - 1U]}};
    for (std::size_t k {a.size() - 1U}; k-- > 0U;)
    {
        const auto m {fx_mul(t, s)};
        s = below ? fx_sub(fx<1> {{a[k]}}, m) : fx_add(fx<1> {{a[k]}}, m);
    }
    return s;
}

// asin(|x|) for 1/256 <= |x| < 1 from sqrt(1 - x^2) and pi/2 - asin(sqrt(1 - x^2)) for |x| >= 0.71, or with comp
// acos(|x|) = pi/2 - asin(|x|) = asin(sqrt(1 - x^2)). Short takes a sqrt with one Newton step fewer and 3 series terms.
template <int N, bool Short = false>
constexpr auto asin_value(const boost::int128::uint128_t& sig, const int e, const int p, const bool comp) noexcept -> fx<N>
{
    const auto u {fx_from_decimal<N>(sig, -e)};
    const bool small {asin_is_small(sig, e, p)};
    const auto a {small ? asin_of<Short>(u, fx_sqrt<N, Short ? 1 : 2>(fx_sub(fx_one<N>(), fx_mul(u, u))))
                         : asin_of<Short>(asin_cos<N, Short>(sig, e), u)};
    return small != comp ? a : fx_sub(fx_top<N>(fx_table::half_pi), a);
}

// Whether |x| = 1 for x = sig * 10^e with p digits and e + p > 0
constexpr auto asin_is_one(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> bool
{
    return e + p == 1 && sig == pow10(static_cast<boost::int128::uint128_t>(p - 1));
}

// asin of x != 0 and not NaN, rounded once in the current mode
template <typename T>
constexpr auto asin_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // |x| >= 1 and infinity have e + p > 0: asin(+-1) = +-pi/2, and NaN above
    if (e + p > 0)
    {
        if (asin_is_one(sig, e, p))
        {
            return pi_quarters<T>(2, neg);
        }
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    // For |x| < 1/256, asin(x) = x * asin(x) / x, and asin(x) / x >= 1 gives digits >= m
    if (decimal_tiny(sig, e))
    {
        const auto u {fx_from_decimal<words>(sig, -e)};
        const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(scale10 - p))};
        return fx_round<T>(static_cast<digits_type>(fx_scale(m, asin_series(fx_mul(u, u)))), e - (scale10 - p), neg);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer for |x| >= 0.1, with an error below 64 units of its last bit, and in one word
        // in (-6.1, 5.1) units of 2^-62 for |x| < 0.71 and in (-5, 5.6) above
        T r {};
        if (decimal_at_least(sig, e, p, 10U, 0)
            && fx_wide_try_round(fx_widen(asin_value<words == 3 ? 2 : 1>(sig, e, p, false)), neg, r, words == 3 ? 6U : 3U))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
    {
        // A first try for |x| >= 0.1, where the result is above 0.1. For every decimal32_t x, asin_poly for |x| < 1/2
        // is within 2^20 units of 2^-62 of the full value, and Short above within 2^27.
        T r {};
        if (decimal_at_least(sig, e, p, 10U, 0))
        {
            const bool poly {!decimal_at_least(sig, e, p, 50U, 0)};
            const auto a {poly ? asin_poly(fx_from_decimal<1>(sig, -e)) : asin_value<1, true>(sig, e, p, false)};
            if (fx_wide_try_round(fx_widen(a), neg, r, poly ? 21U : 28U))
            {
                return r;
            }
        }
    }
    return fx_wide_round<T, words>(fx_widen(asin_value<words>(sig, e, p, false)), neg);
}

} // namespace atan_detail

template <typename T>
constexpr auto asin_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpc {fpclassify(x)};

    if (fpc == FP_ZERO
        #ifndef BOOST_DECIMAL_FAST_MATH
        || fpc == FP_NAN
        #endif
        )
    {
        return issignaling(x) ? nan_conversion(x) : x;
    }

    return atan_detail::asin_fixed(x);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto asin(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::asin_impl(static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ASIN_HPP
