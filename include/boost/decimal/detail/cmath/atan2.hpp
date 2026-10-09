// Copyright 2024 Matt Borland
// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP

#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/cmath/atan.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/check_non_finite.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <cstdint>
#include <cmath>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace atan_detail {

// theta = atan(q) for q = |y / x| = (sy / sx) * 10^k <= 1, else pi/2 - atan(1 / q); pi - theta for x < 0
template <int N, bool Scale = false>
constexpr auto atan2_theta(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const int p, const bool xneg) noexcept -> fx<N>
{
    const bool flat {k < 0 || (k == 0 && sy <= sx)};
    const auto& num {flat ? sy : sx};
    const auto& den {flat ? sx : sy};
    // With Scale both get m = 2^s, so that den * 2^s * 10^-p is in [1/2, 1): each truncation then costs below 2 units
    std::uint64_t m {1U};
    BOOST_DECIMAL_IF_CONSTEXPR (Scale)
    {
        const auto top {pow10(static_cast<boost::int128::uint128_t>(p))};
        m <<= static_cast<unsigned>(boost::int128::countl_zero(den) - boost::int128::countl_zero(top));
        m >>= den * m >= top ? 1U : 0U;
    }
    const auto r {atan_quotient(fx_from_decimal<N>(num * m, flat ? p - k : p + k), fx_from_decimal<N>(den * m, p))};
    const auto theta {flat ? r : fx_sub(fx_top<N>(fx_table::half_pi), r)};
    return xneg ? fx_sub(fx_top<N>(fx_table::pi), theta) : theta;
}

// floor(atan(q) * 10^(s10 - k)) for q = |y / x| = (sy / sx) * 10^k < 1/64 in N words, with s10 in out_s10. With wide
// q < 0.11 in one word, from 10 terms of the series with a tail below 2^-68
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto atan2_small_digits(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const int p, const bool wide, int& out_s10) noexcept -> boost::int128::uint128_t
{
    // atan(q) * 10^(s10 - k) = A - B with A = sy * 10^s10 / sx = d + rem / sx and
    // B = A * (1 - atan(q) / q). Its digits are d - floor(B), less one where frac(B) > rem / sx.
    const int s10 {fx_scale10(N) + (sy < sx ? 1 : 0)};
    boost::int128::uint128_t d {};
    boost::int128::uint128_t rem {};
    BOOST_DECIMAL_IF_CONSTEXPR (N == 1)
    {
        const auto n {sy * pow10(static_cast<boost::int128::uint128_t>(s10))};
        d = n / sx;
        rem = n - d * sx;
    }
    else
    {
        const auto dm {impl::div_mod(umul256(sy, pow10(static_cast<boost::int128::uint128_t>(s10))), sx)};
        d = boost::int128::uint128_t {dm.quotient.bytes[1], dm.quotient.bytes[0]};
        rem = boost::int128::uint128_t {dm.remainder.bytes[1], dm.remainder.bytes[0]};
    }
    // For N = 2, q = d * 10^(k - s10) costs less than a division: the part rem / sx * 10^(k - s10) that it drops
    // is below 1 unit
    const auto q {N == 2 ? fx_from_decimal<N>(d, s10 - k) : decimal_ratio<N>(sy, p - k, sx, p)};
    const auto t {fx_mul(q, q)};
    const auto om {fx_sub(fx_one<N>(), wide ? fx_series<true, 10U, 6U, 11U>(t, fx_table::inv_odd) : atan_series(t))};
    const auto f {decimal_ratio<N>(rem, p, sx, p)};

    // B = d * om + f * om in N + 2 words
    std::uint64_t b[static_cast<std::size_t>(N + 2)] {};
    const auto carry_low {fx_mul_add(fx_make_span<N>(b), fx_make_span(om.w), d.low)};
    b[N] = carry_low;
    const auto carry_high {fx_mul_add(fx_make_span<N>(b, 1), fx_make_span(om.w), d.high)};
    b[N + 1] = carry_high;
    const auto fo {fx_mul(f, om)};
    std::uint64_t carry {fx_add_to(fx_make_span<N>(b), fx_make_span(fo.w))};
    for (int i {N}; i < N + 2; ++i)
    {
        b[i] += carry;
        carry = static_cast<std::uint64_t>(b[i] < carry);
    }
    const boost::int128::uint128_t whole {(boost::int128::uint128_t {b[N + 1], b[N]} << 2U) | boost::int128::uint128_t {b[N - 1] >> 62U}};
    const auto frac {fx_frac<N>(b)};
    out_s10 = s10;
    return d - whole - ((rem == 0U || fx_less(f, frac)) ? 1U : 0U);
}

// atan2(y, x) for finite y and x, both not zero, rounded once in the current mode
template <typename T>
constexpr auto atan2_fixed(const T y, const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    using digits_type = fx_digits_t<words>;

    int ey {};
    int ex {};
    const auto sy {static_cast<boost::int128::uint128_t>(frexp10(y, &ey))};
    const auto sx {static_cast<boost::int128::uint128_t>(frexp10(x, &ex))};
    const bool yneg {signbit(y)};
    const bool xneg {signbit(x)};

    // q = |y / x| = (sy / sx) * 10^k
    const int k {ey - ex};
    const bool below {!xneg && (k <= -3 || (k <= 0 && sy * 64U < sx * pow10(static_cast<boost::int128::uint128_t>(-k))))};
    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        // A first try in one word for 19 digits and q < 0.11: om is 2 units of 2^-62 off and d is below 2.2 * 2^62, so
        // with f om and the compare the digits are 7 units off at most
        if (!xneg && (k < -1 || (k <= 0 && sy * pow10(static_cast<boost::int128::uint128_t>(k + 2)) < sx * 11U)))
        {
            int s1 {};
            const auto d1 {atan2_small_digits<1>(sy, sx, k, p, !below, s1)};
            T r {};
            if (d1 >= pow10(static_cast<boost::int128::uint128_t>(18)) && fx_try_round(static_cast<std::uint64_t>(d1), k - s1, yneg, 8U, r))
            {
                return r;
            }
        }
    }
    if (below)
    {
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // The same in two words, where a window of 16 units covers an error of 7 units in 10^38 digits
            int s2 {};
            const auto d2 {atan2_small_digits<2>(sy, sx, k, p, false, s2)};
            if (k - s2 + 5 >= etiny_v<T> && !fx_near_boundary(d2, 16U))
            {
                return fx_round<T>(static_cast<fx_digits_t<2>>(d2), k - s2, yneg);
            }
        }
        int s10 {};
        const auto digits {atan2_small_digits<words>(sy, sx, k, p, false, s10)};
        return fx_round<T>(static_cast<digits_type>(digits), k - s10, yneg);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // The same, where x < 0 or |y / x| >= 1/64 in two words and |y / x| >= 0.11 in one word. With Scale, 16 units
        // cover an error below 9 units in one word and below 12 units in two words.
        T r {};
        const bool mid {!xneg && (k == -1 || (k == 0 && sy <= sx)) && sy * pow10(static_cast<boost::int128::uint128_t>(k + 2)) >= sx * 11U};
        const bool steep {words == 3 || k > 0 || (k == 0 && sy > sx)};
        if ((xneg || mid || steep) && fx_wide_try_round(fx_widen(atan2_theta<words == 3 ? 2 : 1, true>(sy, sx, k, p, xneg)), yneg, r, 4U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(fx_widen(atan2_theta<words>(sy, sx, k, p, xneg)), yneg);
}

} // namespace atan_detail

template <typename T>
constexpr auto atan2_impl(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpcx {fpclassify(x)};
    const auto fpcy {fpclassify(y)};
    const auto signx {signbit(x)}; // True if neg
    const auto signy {signbit(y)};
    const auto isfinitex {fpcx != FP_INFINITE && fpcx != FP_NAN};
    const auto isfinitey {fpcy != FP_INFINITE && fpcy != FP_NAN};

    T result { };

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (fpcx == FP_NAN)
    {
        result = detail::check_non_finite(x);
    }
    else if (fpcy == FP_NAN)
    {
        result = detail::check_non_finite(y);
    }
    else
    #endif
    if (fpcy == FP_ZERO && signx)
    {
        result = atan_detail::pi_quarters<T>(4, signy);
    }
    else if (fpcy == FP_ZERO && !signx)
    {
        result = y;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpcy == FP_INFINITE && isfinitex)
    {
        result = atan_detail::pi_quarters<T>(2, signy);
    }
    else if (fpcy == FP_INFINITE && fpcx == FP_INFINITE && signx)
    {
        result = atan_detail::pi_quarters<T>(3, signy);
    }
    else if (fpcy == FP_INFINITE && fpcx == FP_INFINITE && !signx)
    {
        result = atan_detail::pi_quarters<T>(1, signy);
    }
    #endif
    else if (fpcx == FP_ZERO)
    {
        result = atan_detail::pi_quarters<T>(2, signy);
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpcx == FP_INFINITE && signx && isfinitey)
    {
        result = atan_detail::pi_quarters<T>(4, signy);
    }
    else if (fpcx == FP_INFINITE && !signx && isfinitey)
    {
        constexpr T zero { 0, 0 };

        result = signy ? -zero : zero;
    }
    #endif
    else
    {
        result = atan_detail::atan2_fixed(y, x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto atan2(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::atan2_impl(static_cast<evaluation_type>(y), static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP
