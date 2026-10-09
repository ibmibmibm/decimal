// Copyright 2024 Matt Borland
// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ATAN2_HPP

#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/cmath/atan.hpp>
#include <boost/decimal/detail/type_traits.hpp>
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
template <int N>
constexpr auto atan2_theta(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const int p, const bool xneg) noexcept -> fx<N>
{
    const auto theta {k < 0 || (k == 0 && sy <= sx) ? atan_quotient(fx_from_decimal<N>(sy, p - k), fx_from_decimal<N>(sx, p))
                                                   : fx_sub(fx_top<N>(fx_table::half_pi), atan_quotient(fx_from_decimal<N>(sx, p + k), fx_from_decimal<N>(sy, p)))};
    return xneg ? fx_sub(fx_top<N>(fx_table::pi), theta) : theta;
}

// atan2(y, x) for finite y and x, both not zero, rounded once in the current mode
template <typename T>
constexpr auto atan2_fixed(const T y, const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int ey {};
    int ex {};
    const auto sy {static_cast<boost::int128::uint128_t>(frexp10(y, &ey))};
    const auto sx {static_cast<boost::int128::uint128_t>(frexp10(x, &ex))};
    const bool yneg {signbit(y)};
    const bool xneg {signbit(x)};

    // q = |y / x| = (sy / sx) * 10^k
    const int k {ey - ex};
    if (!xneg && (k <= -4 || (k <= 0 && sy * 256U < sx * pow10(static_cast<boost::int128::uint128_t>(-k)))))
    {
        // For q < 1/256, atan(q) * 10^(s10 - k) = A - B with A = sy * 10^s10 / sx = d + rem / sx and
        // B = A * (1 - atan(q) / q). Its digits are d - floor(B), less one where frac(B) > rem / sx.
        const int s10 {scale10 + (sy < sx ? 1 : 0)};
        boost::int128::uint128_t d {};
        boost::int128::uint128_t rem {};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
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
        const auto q {decimal_ratio<words>(sy, p - k, sx, p)};
        const auto om {fx_sub(fx_one<words>(), atan_series(fx_mul(q, q)))};
        const auto f {decimal_ratio<words>(rem, p, sx, p)};

        // B = d * om + f * om in words + 2 words
        std::uint64_t b[static_cast<std::size_t>(words + 2)] {};
        b[words] = fx_mul_add(fx_make_span<words>(b), fx_make_span(om.w), d.low);
        b[words + 1] = fx_mul_add(fx_make_span<words>(b, 1), fx_make_span(om.w), d.high);
        const auto fo {fx_mul(f, om)};
        std::uint64_t carry {fx_add_to(fx_make_span<words>(b), fx_make_span(fo.w))};
        for (int i {words}; i < words + 2; ++i)
        {
            b[i] += carry;
            carry = static_cast<std::uint64_t>(b[i] < carry);
        }
        const boost::int128::uint128_t whole {(boost::int128::uint128_t {b[words + 1], b[words]} << 2U) | boost::int128::uint128_t {b[words - 1] >> 62U}};
        const auto frac {fx_frac<words>(b)};
        const auto digits {d - whole - ((rem == 0U || fx_less(f, frac)) ? 1U : 0U)};
        return fx_round<T>(static_cast<digits_type>(digits), k - s10, yneg);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // The same, where x < 0 or |y / x| > 0.01 in two words and |y / x| >= 1.6 in one word. One word has 32 units
        // for x > 0 and 0.11 <= |y / x| <= 1.
        T r {};
        const bool mid {words == 2 && !xneg && (k == -1 || (k == 0 && sy <= sx)) && sy * pow10(static_cast<boost::int128::uint128_t>(k + 2)) >= sx * 11U};
        if ((xneg || mid || (words == 3 ? k >= -1 : k > 1 || (k >= 0 && sy * pow10(static_cast<boost::int128::uint128_t>(k + 1)) >= sx * 16U)))
            && fx_wide_try_round(fx_widen(atan2_theta<words == 3 ? 2 : 1>(sy, sx, k, p, xneg)), yneg, r, mid ? 5U : 6U))
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
        result = issignaling(x) ? nan_conversion(x) : x;
    }
    else if (fpcy == FP_NAN)
    {
        result = issignaling(y) ? nan_conversion(y) : y;
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
