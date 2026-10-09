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

// B = d g + fg in N + 2 words for the fractions g and fg: the whole part and the fraction of B
template <int N>
struct atan2_parts
{
    boost::int128::uint128_t whole;
    fx<N> frac;
};

template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto atan2_split(const boost::int128::uint128_t& d, const fx<N>& g, const fx<N>& fg) noexcept -> atan2_parts<N>
{
    std::uint64_t b[static_cast<std::size_t>(N + 2)] {};
    const auto carry_low {fx_mul_add(fx_make_span<N>(b), fx_make_span(g.w), d.low)};
    b[N] = carry_low;
    const auto carry_high {fx_mul_add(fx_make_span<N>(b, 1), fx_make_span(g.w), d.high)};
    b[N + 1] = carry_high;
    std::uint64_t carry {fx_add_to(fx_make_span<N>(b), fx_make_span(fg.w))};
    for (int i {N}; i < N + 2; ++i)
    {
        b[i] += carry;
        carry = static_cast<std::uint64_t>(b[i] < carry);
    }
    return {(boost::int128::uint128_t {b[N + 1], b[N]} << 2U) | boost::int128::uint128_t {b[N - 1] >> 62U}, fx_frac<N>(b)};
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
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(N == 1))
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

    // B = d * om + f * om
    const auto b {atan2_split(d, om, fx_mul(f, om))};
    out_s10 = s10;
    return d - b.whole - ((rem == 0U || fx_less(f, b.frac)) ? 1U : 0U);
}

// atan2_small_digits where B is next to an integer or to rem / sx: B = A t^2 R 10^-2a for q = t 10^-a with t in
// [0.1, 1) and R = 1/3 - q^2/5 + ..., so that B has its own precision. A tie goes down, as B is below its lead term.
template <int N>
BOOST_DECIMAL_NO_INLINE constexpr auto atan2_small_lead_digits(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const int p) noexcept -> boost::int128::uint128_t
{
    using u128 = boost::int128::uint128_t;
    const int s10 {fx_scale10(N) + (sy < sx ? 1 : 0)};
    const auto dm {impl::div_mod(umul256(sy, pow10(static_cast<u128>(s10))), sx)};
    const u128 d {dm.quotient.bytes[1], dm.quotient.bytes[0]};
    const u128 rem {dm.remainder.bytes[1], dm.remainder.bytes[0]};
    const int a {-k - (sy < sx ? 0 : 1)};
    const auto t {decimal_ratio<N>(sy, p + (sy < sx ? 0 : 1), sx, p)};
    const auto q {decimal_ratio<N>(sy, p - k, sx, p)};
    const auto g {fx_mul(fx_mul(t, t), fx_horner_range<true>(fx_mul(q, q), fx_table::inv_odd, 1U, 16U, fx<N> {}))};

    // B 10^2a = d g + (rem / sx) g
    const auto b {atan2_split(d, g, fx_mul(decimal_ratio<N>(rem, p, sx, p), g))};
    auto w {b.whole};
    const auto& frac {b.frac};

    // Past 38 digits the units are 10^(2a - 38) times larger, with the digits below them in low, and B is below 1
    const bool coarse {2 * a > 38};
    const auto cut {pow10(static_cast<u128>(coarse ? 2 * a - 38 : 0))};
    const auto low {w % cut};
    w /= cut;
    const auto unit {pow10(static_cast<u128>(coarse ? 38 : 2 * a))};
    auto whole {w / unit};
    const auto part {w - whole * unit};
    // rem / sx in the same units, as an integer and a fraction
    const auto fm {impl::div_mod(umul256(rem, unit), sx)};
    u128 f_rem {fm.remainder.bytes[1], fm.remainder.bytes[0]};
    bool above {part > u128 {fm.quotient.bytes[1], fm.quotient.bytes[0]}};
    bool equal {part == u128 {fm.quotient.bytes[1], fm.quotient.bytes[0]}};
    if (coarse && equal)
    {
        const auto fl {impl::div_mod(umul256(f_rem, cut), sx)};
        above = low > u128 {fl.quotient.bytes[1], fl.quotient.bytes[0]};
        equal = low == u128 {fl.quotient.bytes[1], fl.quotient.bytes[0]};
        f_rem = u128 {fl.remainder.bytes[1], fl.remainder.bytes[0]};
    }
    fx<N> f_frac {decimal_ratio<N>(f_rem, p, sx, p)};
    // The error of the sum is below 2^-56 units
    f_frac.w[N - 1] += 64U;
    bool borrow {rem == 0U || above || (equal && fx_less(f_frac, frac))};
    if (!coarse && whole != 0U && part == 0U && frac.w[N - 1] < 64U)
    {
        --whole;
        borrow = true;
    }
    return d - whole - (borrow ? 1U : 0U);
}

// atan2(y, x) for x > 0 and a tiny q = |y / x| = (sy / sx) * 10^k: atan(q) is below A = q * 10^(s10 - k) by less than
// 1 / sx, thus its digits are those of A, less one unit where A is an integer
template <typename T>
constexpr auto atan2_tiny(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const bool yneg) noexcept -> T
{
    constexpr int words {fx_words(std::numeric_limits<T>::digits10)};
    const int s10 {fx_scale10(words == 3 ? 2 : 1) + (sy < sx ? 1 : 0)};
    boost::int128::uint128_t d {};
    bool exact {};
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 3))
    {
        const auto dm {impl::div_mod(umul256(sy, pow10(static_cast<boost::int128::uint128_t>(s10))), sx)};
        d = boost::int128::uint128_t {dm.quotient.bytes[1], dm.quotient.bytes[0]};
        exact = (dm.remainder.bytes[0] | dm.remainder.bytes[1]) == 0U;
    }
    else
    {
        const auto n {sy * pow10(static_cast<boost::int128::uint128_t>(s10))};
        d = n / sx;
        exact = d * sx == n;
    }
    // Two words take 18 digits more
    constexpr int more {words == 2 ? 18 : 0};
    const auto digits {d * pow10(static_cast<boost::int128::uint128_t>(more)) - (exact ? 1U : 0U)};
    return fx_round<T>(static_cast<fx_digits_t<words>>(digits), k - s10 - more, yneg);
}

// Whether 37 or 38 digits with an error below w units are next to a rounding boundary of p digits
constexpr auto atan2_near_boundary(const boost::int128::uint128_t& digits, const int p, const std::uint64_t w = 16U) noexcept -> bool
{
    const auto short_by {digits < pow10(static_cast<boost::int128::uint128_t>(37)) ? 1 : 0};
    const auto m {pow10(static_cast<boost::int128::uint128_t>(37 - short_by - p)) * 5U};
    const auto r {digits % m};
    return r < w || r > m - w;
}

// A table value in six words from its three words top and the three words low below them
constexpr auto atan2_join(const fx<3>& top, const fx<3>& low) noexcept -> fx<6>
{
    return {{low.w[0], low.w[1], low.w[2], top.w[0], top.w[1], top.w[2]}};
}

// atan2 in six words, rounded once, for decimal128_t where three words leave the digits next to a rounding boundary:
// two free significands put y / x that close to it. q = num / (den * 10^j) <= 1 with j <= 38, as in atan_quotient.
template <typename T>
BOOST_DECIMAL_NO_INLINE constexpr auto atan2_wide(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const bool xneg, const bool yneg) noexcept -> T
{
    const bool flat {k < 0 || (k == 0 && sy <= sx)};
    const auto& num {flat ? sy : sx};
    const auto wide_den {umul256(flat ? sx : sy, pow10(static_cast<boost::int128::uint128_t>(flat ? -k : k)))};
    fx<6> b {{wide_den.bytes[0], wide_den.bytes[1], wide_den.bytes[2], wide_den.bytes[3], 0U, 0U}};
    const int s {fx_countl_zero(b) - 2};
    b = fx_shift(b, s);
    const auto a {fx_shift(fx<6> {{num.low, num.high, 0U, 0U, 0U, 0U}}, s)};

    const auto i {atan_index(a, b)};
    const auto cb {fx_mul_shr(b, i, 5U)};
    const bool under {fx_less(a, cb)};
    const auto n {under ? fx_sub(cb, a) : fx_sub(a, cb)};
    const auto den {fx_add(b, fx_mul_shr(a, i, 5U))};
    const auto d {static_cast<unsigned>(countl_zero(den.w[5]) - 1)};
    const auto q {d == 0U ? fx_div(n, den) : fx_div(fx_shl(n, d), fx_shl(den, d))};

    // atan(q) / q from 32 terms for q^2 <= 2^-12
    const auto t {fx_mul(q, q)};
    fx<6> series {};
    for (std::uint64_t j {32U}; j-- > 0U;)
    {
        series = fx_sub(fx_div_small(fx_one<6>(), 2U * j + 1U), fx_mul(t, series));
    }
    const auto r {fx_mul(q, series)};
    const auto c {atan2_join(fx_table::atan_c[i], fx_table::atan_c_low[i])};
    auto theta {under ? fx_sub(c, r) : fx_add(c, r)};
    if (!flat)
    {
        theta = fx_sub(atan2_join(fx_table::half_pi, fx_table::half_pi_low), theta);
    }
    if (xneg)
    {
        theta = fx_sub(atan2_join(fx_table::pi, fx_table::pi_low), theta);
    }
    return fx_wide_round<T, 6>(fx_widen(theta), yneg);
}

// atan2(y, x) for finite y and x, both not zero, rounded once in the current mode
template <typename T>
constexpr auto atan2_fixed(const T y, const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int ey {};
    int ex {};
    const auto sy {static_cast<boost::int128::uint128_t>(frexp10(y, &ey))};
    const auto sx {static_cast<boost::int128::uint128_t>(frexp10(x, &ex))};
    const bool yneg {signbit(y)};
    const bool xneg {signbit(x)};
    if (!xneg && ex == 1 - p && sx == pow10(static_cast<boost::int128::uint128_t>(p - 1)))
    {
        return atan(y);
    }

    // q = |y / x| = (sy / sx) * 10^k
    const int k {ey - ex};
    if (!xneg && k <= -(fx_scale10(words == 3 ? 2 : 1) + p + 4) / 2)
    {
        return atan2_tiny<T>(sy, sx, k, yneg);
    }
    const bool below {!xneg && (k <= -3 || (k <= 0 && sy * 64U < sx * pow10(static_cast<boost::int128::uint128_t>(-k))))};
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words < 3))
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
            // For q < 10^-10 and an integer A = d1 + 1, atan(q) * 10^(s1 - k) is below A by less than 1 unit
            if (k <= -11 && (d1 + 1U) * sx == sy * pow10(static_cast<boost::int128::uint128_t>(s1)))
            {
                return fx_round<T>((d1 + 1U) * pow10(static_cast<boost::int128::uint128_t>(18)) - 1U, k - s1 - 18, yneg);
            }
        }
    }
    if (below)
    {
        BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 3))
        {
            // The same in two words, where a window of 16 units covers an error of 7 units in 10^38 digits
            int s2 {};
            const auto d2 {atan2_small_digits<2>(sy, sx, k, p, false, s2)};
            // The same for q < 10^-19 and an integer A = d2 + 1: atan(q) * 10^(s2 - k) is below A by less than 1 unit
            if ((k - s2 + 5 >= etiny_v<T> && !fx_near_boundary(d2, 16U))
                || (k <= -20 && umul256(d2 + 1U, sx) == umul256(sy, pow10(static_cast<boost::int128::uint128_t>(s2)))))
            {
                return fx_round<T>(static_cast<fx_digits_t<2>>(d2), k - s2, yneg);
            }
        }
        // In three words, or in two after the first try in one. Next to a rounding boundary, A can be next to it by
        // B, or a short q puts B next to an integer or to rem / sx: B then takes its own precision.
        constexpr int last {words == 3 ? 3 : 2};
        int s10 {};
        auto digits {atan2_small_digits<last>(sy, sx, k, p, false, s10)};
        BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 3))
        {
            // For q >= 10^-15, y / x can put B closer to rem / sx than three words tell
            if (k >= -14 && fx_near_boundary(digits, 2U))
            {
                return atan2_wide<T>(sy, sx, k, false, yneg);
            }
        }
        if (words == 3 ? fx_near_boundary(digits, 2U) : atan2_near_boundary(digits, p))
        {
            digits = atan2_small_lead_digits<3>(sy, sx, k, p);
        }
        return fx_round<T>(digits, k - s10, yneg);
    }

    // The same, where x < 0 or |y / x| >= 1/64 in two words and |y / x| >= 0.11 in one word. With Scale, 16 units
    // cover an error below 9 units in one word and below 12 units in two words.
    T r {};
    const bool mid {!xneg && (k == -1 || (k == 0 && sy <= sx)) && sy * pow10(static_cast<boost::int128::uint128_t>(k + 2)) >= sx * 11U};
    const bool steep {words == 3 || k > 0 || (k == 0 && sy > sx)};
    if ((xneg || mid || steep) && fx_unsigned_try_round(atan2_theta<words == 3 ? 2 : 1, true>(sy, sx, k, p, xneg), yneg, r, 4U))
    {
        return r;
    }
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 3))
    {
        // For |y / x| >= 10^12, the result is pi/2 -+ atan(|x / y|), and for x < 0 and |y / x| <= 10^-12 it is
        // pi - atan(|y / x|): the small quotient q = t * 10^-a has t in [0.1, 1)
        if (k >= 13 || (xneg && k <= -13))
        {
            const bool flat {k < 0};
            const auto& num {flat ? sy : sx};
            const auto& den {flat ? sx : sy};
            const int a {(flat ? -k : k) - (num < den ? 0 : 1)};
            const auto g {atan_small_factor(decimal_ratio<3>(num, p + (num < den ? 0 : 1), den, p), a)};
            const auto digits {pi_digits(flat, g, a, !flat && xneg)};
            // A q below 10^-38 leaves the digits of pi/2 and pi, which are not next to a rounding boundary
            return fx_near_boundary(digits, 2U) && a <= 37 ? atan2_wide<T>(sy, sx, k, xneg, yneg) : fx_round<T>(digits, -37, yneg);
        }
    }
    // Three words also for decimal64_t: next to pi/2 or pi, y / x puts the result closer to a tie than two words tell.
    // Two words first, which decide but for a result next to a rounding boundary: their error is far below 2^20 units.
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 2))
    {
        bool neg2 {yneg};
        int e2 {};
        const auto digits2 {fx_wide_digits(fx_widen(atan2_theta<2>(sy, sx, k, p, xneg)), neg2, e2)};
        if (!atan2_near_boundary(digits2, p, UINT64_C(1) << 20U))
        {
            return fx_round<T>(digits2, e2, neg2);
        }
    }
    constexpr int last {words == 1 ? 2 : 3};
    bool neg {yneg};
    int exp10 {};
    const auto digits {fx_wide_digits(fx_widen(atan2_theta<last>(sy, sx, k, p, xneg)), neg, exp10)};
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(words == 3))
    {
        if (fx_near_boundary(digits, 2U))
        {
            return atan2_wide<T>(sy, sx, k, xneg, yneg);
        }
    }
    return fx_round<T>(digits, exp10, neg);
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
    #ifndef BOOST_DECIMAL_FAST_MATH
    const auto isfinitex {fpcx != FP_INFINITE && fpcx != FP_NAN};
    const auto isfinitey {fpcy != FP_INFINITE && fpcy != FP_NAN};
    #endif

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
