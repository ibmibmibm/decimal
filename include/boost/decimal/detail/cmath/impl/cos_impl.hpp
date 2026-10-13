// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/impl/trig_reduce.hpp>
#include <boost/decimal/detail/cmath/impl/sin_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace trig {

constexpr auto cos_poly(const fx<1>& z) noexcept -> fx<1>
{
    return fx_alternating(z, fx_table::cos_d32);
}

constexpr auto cos_poly(const fx<2>& z) noexcept -> fx<2>
{
    return fx_alternating_split<fx_table::cos_d64_shift>(z, fx_table::cos_d64_head, fx_table::cos_d64_tail);
}

constexpr auto cos_poly(const fx<3>& z) noexcept -> fx<3>
{
    return fx_alternating_split<fx_table::cos_d128_mid_shift, fx_table::cos_d128_tail_shift>(z, fx_table::cos_d128_head, fx_table::cos_d128_mid, fx_table::cos_d128_tail);
}

// The polynomial of decimal64_t in one word
constexpr auto cos_poly_top(const fx<1>& z) noexcept -> fx<1>
{
    return fx_alternating_top<fx_table::cos_d64_shift>(z, fx_table::cos_d64_head, fx_table::cos_d64_tail);
}

// The polynomial of decimal128_t in two words
constexpr auto cos_poly_top(const fx<2>& z) noexcept -> fx<2>
{
    return fx_alternating_top<fx_table::cos_d128_mid_shift, fx_table::cos_d128_tail_shift>(z, fx_table::cos_d128_head, fx_table::cos_d128_mid, fx_table::cos_d128_tail);
}

// cos(r) for z = r^2 and the polynomial p, which is below 1 for any r != 0.
template <int N>
constexpr auto cos_fx(const fx<N>& z, const fx<N>& p) noexcept -> fx<N>
{
    fx<N> half {};
    half.w[N - 1] = UINT64_C(1) << 61U;
    const auto inner {fx_sub(half, fx_mul(z, p))};
    return fx_clamp_below_one(fx_sub(fx_one<N>(), fx_mul(z, inner)));
}

template <int N>
constexpr auto cos_fx(const fx<N>& z) noexcept -> fx<N>
{
    return cos_fx(z, cos_poly(z));
}

// sin(r) / r = 1 - z/(2*3) * (1 - z/(4*5) * (1 - ...)), or with use_cos cos(r) = 1 - z/(1*2) * (1 - ...), to 22 terms
// with an error below 2^-186 for z <= (pi/4)^2 and no table. With last = 1: the sum that the first z/(...) multiplies.
template <int N>
constexpr auto trig_taylor(const fx<N>& z, const bool use_cos, const std::uint64_t last = 0U) noexcept -> fx<N>
{
    const std::uint64_t c {use_cos ? 1U : 0U};
    auto s {fx_one<N>()};
    for (std::uint64_t k {22}; k > last; --k)
    {
        s = fx_sub(fx_one<N>(), fx_div_small(fx_mul(z, s), (2U * k - c) * (2U * k + 1U - c)));
    }
    return fx_clamp_below_one(s);
}

// sin(r) / r for z = r^2 from the Taylor series with taylor, or else from the polynomial
template <int N>
constexpr auto sinc_fx_of(const fx<N>& z, const bool taylor) noexcept -> fx<N>
{
    return taylor ? trig_taylor(z, false) : sinc_fx(z);
}

// cos(r) for z = r^2 from the Taylor series with taylor, or else from the polynomial
template <int N>
constexpr auto cos_fx_of(const fx<N>& z, const bool taylor) noexcept -> fx<N>
{
    return taylor ? trig_taylor(z, true) : cos_fx(z);
}

// The signs of sin(x) and cos(x) for x = n*pi/2 + r, where the sign of sin(r) is r.neg and cos(r) > 0
template <int N>
constexpr auto sin_neg(const trig_arg<N>& r) noexcept -> bool
{
    return (r.n & 1U) != 0U ? r.n == 3U : r.neg != (r.n == 2U);
}

template <int N>
constexpr auto cos_neg(const trig_arg<N>& r) noexcept -> bool
{
    return (r.n & 1U) != 0U ? r.neg != (r.n == 1U) : r.n == 2U;
}

// The digits of sin_cos_of from the Taylor series in one word, out of line as they are rare
template <int N>
BOOST_DECIMAL_NO_INLINE constexpr auto trig_taylor_digits(const trig_arg<N>& r, const fx<N>& z, const bool use_cos) noexcept -> boost::int128::uint128_t
{
    const auto f {trig_taylor(z, use_cos)};
    return use_cos ? fx_scale(pow10(static_cast<boost::int128::uint128_t>(sig_digits(N))), f) : fx_scale(r.sig, r.tail, f);
}

// floor(sin(r) * 10^k), or floor(cos(r) * 10^digits) with use_cos, for an exact r = sig * 10^-k: sin(r) = r (1 - z p)
// and cos(r) = 1 - z q, with p and q from the Taylor series, and the rest past the lead term takes its own precision
template <int N>
BOOST_DECIMAL_NO_INLINE constexpr auto sin_cos_lead_digits(const trig_arg<N>& r, const fx<N>& z, const bool use_cos) noexcept -> boost::int128::uint128_t
{
    int a {};
    const auto t {fx_lead_t<N>(r.sig, r.k, a)};
    const auto g {fx_mul(fx_mul(t, t), fx_div_small(trig_taylor(z, use_cos, 1U), use_cos ? 2U : 6U))};
    return fx_lead_digits(use_cos ? pow10(static_cast<boost::int128::uint128_t>(sig_digits(N))) : r.sig, g, 2 * a, true);
}

// sin(|r|), or with use_cos cos(r), with the sign neg, for r != 0. The polynomials are good to 2 units of 38 digits for
// three words and to a relative 3 * 10^-23 for two, thus near a rounding boundary the Taylor series decides.
template <typename T>
constexpr auto sin_cos_of(const trig_arg<trig_traits<T>::words>& r, const bool use_cos, const bool neg) noexcept -> T
{
    constexpr int digits {sig_digits(trig_traits<T>::words)};
    constexpr auto one {pow10(static_cast<boost::int128::uint128_t>(digits))};
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    auto v {use_cos ? fx_scale(one, cos_fx(z)) : fx_scale(r.sig, r.tail, sinc_fx(z))};

    // Below z = 2^-62 the polynomials are as good as the series for an r that has an error. For an exact r the
    // series gives the rest past the lead term.
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(trig_traits<T>::words > 1))
    {
        const bool exact {!r.fixed && r.tail == 0U};
        if ((z.w[trig_traits<T>::words - 1] != 0U || exact) && trig_near_boundary<trig_traits<T>::words>(v))
        {
            if (exact)
            {
                v = sin_cos_lead_digits(r, z, use_cos);
            }
            else
            {
                const auto rt {fixed_r_with_tail(r)};
                const auto f {trig_taylor(fx_mul(rt, rt), use_cos)};
                v = use_cos ? fx_scale(one, f) : fx_scale(r.sig, r.tail, f);
            }
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (fx_is(trig_traits<T>::words == 1))
    {
        // In one word the polynomials are good to a relative 5 * 10^-15: in a window of 2^-40 v from a rounding
        // boundary the Taylor series decides
        if (trig_near_boundary_word(static_cast<std::uint64_t>(v), 40U))
        {
            v = trig_taylor_digits(r, z, use_cos);
        }
    }
    return trig_round<T>(v, use_cos ? -digits : -r.k, neg);
}

// sin(r) from the digits of r, with the arguments of trig_reduce_word_digits, and s = sin(r) / r with the error of
// sin_cos_fixed_try. It is out of line with its own trig_arg: this keeps the one of the caller out of memory.
template <typename T>
BOOST_DECIMAL_NO_INLINE constexpr auto sin_digits_try(const boost::int128::uint128_t& high, const std::uint64_t low, const unsigned sh, const fx<1>& s, const bool neg, T& out) noexcept -> bool
{
    trig_arg<1> r {};
    trig_reduce_word_digits(high, low, sh, r);
    return trig_try_round(static_cast<std::uint64_t>(fx_scale(r.sig, r.tail, s)), -r.k, neg, UINT64_C(1) << 16U, out);
}

// A first try of sin(x), or with cosine of cos(x), for decimal32_t and 10^-4 <= |x| < 10^19 from r in fixed point, or
// from its digits for sin(r) < 0.1, rounded once in out where it cannot be next to a rounding boundary.
template <typename T>
constexpr auto sin_cos_fixed_try(const T x, const bool cosine, T& out) noexcept -> bool
{
    constexpr int prec {precision_v<T>};
    constexpr int up {sig_digits(1) - prec};
    const bool neg {signbit(x)};
    int e {};
    const auto m {frexp10(neg ? -x : x, &e)};
    if (m == 0U || prec + e > 19 || prec + e <= -fx_tiny_digits(prec))
    {
        return false;
    }

    // The reduction gives only the quadrant, the sign and rf, with an error below 2 units.
    trig_arg<1> r {};
    r.neg = neg;
    const bool reduced {e > -prec || (e == -prec && m > quarter_pi_digits(sig_digits(1)).low / pow10(static_cast<std::uint64_t>(up)))};
    std::uint64_t big[4] {};
    unsigned sh {};
    if (reduced)
    {
        sh = trig_reduce_word_fixed(m, e, neg, r, big);
    }
    else
    {
        r.rf = fx_from_decimal<1>(m, -e);
    }

    // The polynomials are good to 5 * 10^-15: the errors are below 2^16 units of 2^-62, and of the 19 digits of sin.
    const auto z {fx_mul(r.rf, r.rf)};
    const bool sign {cosine ? cos_neg(r) : sin_neg(r)};
    if (cosine == ((r.n & 1U) == 0U))
    {
        return fx_unsigned_try_round(cos_fx(z), sign, out, 16U);
    }
    if (!reduced)
    {
        const boost::int128::uint128_t sig {static_cast<std::uint64_t>(m) * pow10(static_cast<std::uint64_t>(up))};
        return trig_try_round(static_cast<std::uint64_t>(fx_scale(sig, sinc_fx(z))), e - up, sign, UINT64_C(1) << 16U, out);
    }

    // sin(r) >= 0.1 has 19 digits. Below it, the digits of r give them.
    if (r.rf.w[0] >= UINT64_C(0x0669AD42C3C9EECC))
    {
        return fx_unsigned_try_round(fx_mul(r.rf, sinc_fx(z)), sign, out, 16U);
    }
    return sin_digits_try(boost::int128::uint128_t {big[3], big[2]}, big[1], sh, sinc_fx(z), sign, out);
}

// sin_cos_of in two words for decimal128_t, with the polynomials of the type: the error of v is in (-7.8, 11.5) units
// for sin(r) and in (-6.9, 12.1) for cos(r), thus w = 13. It rounds once in out where v cannot be next to a boundary.
template <typename T, int N>
constexpr auto sin_cos_try(const trig_arg<N>& r, const bool use_cos, const bool neg, T& out) noexcept -> bool
{
    constexpr int digits {sig_digits(N)};
    constexpr auto one {pow10(static_cast<boost::int128::uint128_t>(digits))};
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    const auto v {use_cos ? fx_scale(one, cos_fx(z, cos_poly_top(z))) : fx_scale(r.sig, r.tail, sinc_fx(z, sin_poly_top(z)))};
    return trig_try_round(static_cast<fx_digits_t<N>>(v), use_cos ? -digits : -r.k, neg, 13U, out);
}

// 8 (1 - cos(r)) as sin_rest_top at 8 times the scale, where the inner sum 4 - z p is below 4 and wraps at 4: the
// error is in (-12, 1.4) units.
constexpr auto cos_rest_top(const fx<1>& z) noexcept -> fx<1>
{
    fx<1> s {{fx_tail_sum(z.w[0], fx_table::cos_d64_tail) >> (fx_table::cos_d64_shift - 3)}};
    s = fx_sub(fx_top_shl(fx_table::cos_d64_head[1], 3U), fx_mul(z, s));
    s = fx_sub(fx_top_shl(fx_table::cos_d64_head[0], 3U), fx_mul(z, s));
    return fx_mul(z, fx_sub(fx<1> {}, fx_mul(z, s)));
}

// sin_cos_try in one word for decimal64_t, from sin_rest_top and cos_rest_top: the error of v is in (-1.7, 1.8) units
// for sin(r), which leaves out the tail of r, thus w = 2, and in (-2.4, 2.3) for cos(r), thus w = 3.
template <typename T>
constexpr auto sin_cos_try(const trig_arg<1>& r, const bool use_cos, const bool neg, T& out) noexcept -> bool
{
    using wide = boost::int128::uint128_t;
    constexpr auto one {pow10(static_cast<std::uint64_t>(19))};
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    if (use_cos)
    {
        return trig_try_round(one - 2U - ((static_cast<wide>(one) * cos_rest_top(z).w[0]) >> 65U).low, -19, neg, 3U, out);
    }
    // sig (1 - sin(r) / r), rounded to the nearest
    const auto rest {(static_cast<wide>(r.sig.low) * sin_rest_top(z).w[0] + wide {UINT64_C(2), UINT64_C(0)}) >> 66U};
    return trig_try_round(r.sig.low - rest.low, -r.k, neg, 2U, out);
}

} // namespace trig
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP
