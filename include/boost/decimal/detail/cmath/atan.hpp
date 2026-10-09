// Copyright 2023 - 2024 Matt Borland
// Copyright 2023 - 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ATAN_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ATAN_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/u256.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace atan_detail {

// atan(q) / q = 1 - t/3 + t^2/5 - ... for t = q^2 <= 2^-12, to below 2^-75, 2^-136 and 2^-197 in one, two
// and three words. The terms from t^6 and t^11 on are below 2^-72 and 2^-132.
template <int N>
constexpr auto atan_series(const fx<N>& t) noexcept -> fx<N>
{
    return fx_series<true, N == 1 ? 6U : N == 2 ? 11U : 16U, 6U, 11U>(t, fx_table::inv_odd);
}

// atan(z) for 0 <= z <= 1. With c = i/32 next to z, atan(z) = atan(c) + atan(q) for
// q = (z - c) / (1 + z c), and |q| <= 1/64.
template <int N>
constexpr auto atan_of(const fx<N>& z) noexcept -> fx<N>
{
    const auto i {static_cast<std::size_t>((z.w[N - 1] + (UINT64_C(1) << 56U)) >> 57U)};
    fx<N> c {};
    c.w[N - 1] = static_cast<std::uint64_t>(i) << 57U;
    const bool below {fx_less(z, c)};
    const auto q {fx_div(below ? fx_sub(c, z) : fx_sub(z, c), fx_add(fx_one<N>(), fx_mul_shr(z, i, 5U)))};
    const auto s {fx_mul(q, atan_series(fx_mul(q, q)))};
    const auto a {fx_top<N>(fx_table::atan_c[i])};
    return below ? fx_sub(a, s) : fx_add(a, s);
}

// i with i/32 next to a / b <= 1 + 1/64, by a short division of the top words
template <int N>
constexpr auto atan_index(const fx<N>& a, const fx<N>& b) noexcept -> std::size_t
{
    const auto ta {a.w[N - 1] >> 8U};
    const auto tb {b.w[N - 1] >> 8U};
    return static_cast<std::size_t>((64U * ta + tb) / (2U * tb));
}

// atan(a / b) for 0 <= a <= b < 1 with b > 2^-54, with one reciprocal: c = i/32 next to a / b comes from a short
// division of the top words, and atan(a / b) = atan(c) + atan(q) for q = (a - c b) / (b + c a), with |q| <= 1/64
template <int N>
constexpr auto atan_quotient(fx<N> a, fx<N> b) noexcept -> fx<N>
{
    // b in [1/2, 1), and b + c a in [1/2, 2)
    const auto s {static_cast<unsigned>(countl_zero(b.w[N - 1]) - 2)};
    if (s != 0U)
    {
        a = fx_shl(a, s);
        b = fx_shl(b, s);
    }
    const auto i {atan_index(a, b)};
    const auto cb {fx_mul_shr(b, i, 5U)};
    const bool below {fx_less(a, cb)};
    const auto num {below ? fx_sub(cb, a) : fx_sub(a, cb)};
    const auto den {fx_add(b, fx_mul_shr(a, i, 5U))};
    const auto d {static_cast<unsigned>(countl_zero(den.w[N - 1]) - 1)};
    const auto q {d == 0U ? fx_div(num, den) : fx_div(fx_shl(num, d), fx_shl(den, d))};
    const auto r {fx_mul(q, atan_series(fx_mul(q, q)))};
    const auto c {fx_top<N>(fx_table::atan_c[i])};
    return below ? fx_sub(c, r) : fx_add(c, r);
}

// a / b for a = sa * 10^-ka < b = sb * 10^-p, with sb of p digits
template <int N>
constexpr auto atan_ratio(const boost::int128::uint128_t& sa, const int ka, const boost::int128::uint128_t& sb, const int p) noexcept -> fx<N>
{
    unsigned s {};
    const auto v {fx_normalize(fx_from_decimal<N>(sb, p), s)};
    return fx_div(fx_shl(fx_from_decimal<N>(sa, ka), s), v);
}

// k * pi / 4 for k = 1 to 4, rounded once in the current mode with the sign neg
template <typename T>
constexpr auto pi_quarters(const int k, const bool neg) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    const auto half_pi {fx_top<words>(fx_table::pi_2)};
    const auto pi {fx_top<words>(fx_table::pi)};
    const auto v {k == 1 ? fx_shr(half_pi, 1U) : k == 2 ? half_pi : k == 3 ? fx_sub(pi, fx_shr(half_pi, 1U)) : pi};
    return fx_wide_round<T, words>(fx_widen(v), neg);
}

// atan(|x|) for x = sig * 10^e with |x| >= 1/256: atan_of below 1, else pi/2 - atan(1 / |x|) with 1 / |x| =
// 10^-(e + p) / u and u = sig * 10^-p
template <int N>
constexpr auto atan_value(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> fx<N>
{
    return e + p <= 0 ? atan_of(fx_from_decimal<N>(sig, -e))
                      : fx_sub(fx_top<N>(fx_table::pi_2), atan_quotient(fx_from_decimal<N>(1U, e + p), fx_from_decimal<N>(sig, p)));
}

// theta = atan(q) for q = |y / x| = (sy / sx) * 10^k <= 1, else pi/2 - atan(1 / q); pi - theta for x < 0
template <int N>
constexpr auto atan2_theta(const boost::int128::uint128_t& sy, const boost::int128::uint128_t& sx, const int k, const int p, const bool xneg) noexcept -> fx<N>
{
    const auto theta {k < 0 || (k == 0 && sy <= sx) ? atan_quotient(fx_from_decimal<N>(sy, p - k), fx_from_decimal<N>(sx, p))
                                                   : fx_sub(fx_top<N>(fx_table::pi_2), atan_quotient(fx_from_decimal<N>(sx, p + k), fx_from_decimal<N>(sy, p)))};
    return xneg ? fx_sub(fx_top<N>(fx_table::pi), theta) : theta;
}

// atan of a finite x != 0, rounded once in the current mode
template <typename T>
constexpr auto atan_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // For |x| < 1/256, atan(x) = x * atan(x) / x < |x|, also where the series rounds to 1
    if (decimal_tiny(sig, e))
    {
        const auto u {fx_from_decimal<words>(sig, -e)};
        return fx_round_below_x<T>(sig, e, atan_series(fx_mul(u, u)), neg);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer, with an error below 64 units of its last bit, for |x| >= 0.1 in two words
        // and |x| >= 1.6 in one word, where the result is above 1. One word has 8 units for 0.11 <= |x| < 1.
        T r {};
        const bool below_one {words == 2 && e + p == 0 && decimal_at_least(sig, e, p, 11U, 0)};
        if ((words == 3 ? e + p >= 0 : below_one || decimal_at_least(sig, e, p, 16U, 1))
            && fx_wide_try_round(fx_widen(atan_value<words == 3 ? 2 : 1>(sig, e, p)), neg, r, below_one ? 3U : 6U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(fx_widen(atan_value<words>(sig, e, p)), neg);
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
        const auto q {atan_ratio<words>(sy, p - k, sx, p)};
        const auto om {fx_sub(fx_one<words>(), atan_series(fx_mul(q, q)))};
        const auto f {atan_ratio<words>(rem, p, sx, p)};

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
constexpr auto atan_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    const auto fpc { fpclassify(x) };

    T result { };

    if (fpc == FP_ZERO
        #ifndef BOOST_DECIMAL_FAST_MATH
        || fpc == FP_NAN
        #endif
        )
    {
        result = x;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpc == FP_INFINITE)
    {
        result = atan_detail::pi_quarters<T>(2, signbit(x));
    }
    #endif
    else
    {
        result = atan_detail::atan_fixed(x);
    }

    return result;
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto atan(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::atan_impl(static_cast<evaluation_type>(x)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ATAN_HPP
