// Copyright 2024 Matt Borland
// Copyright 2002 - 2011, 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_ELLINT_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_ELLINT_IMPL_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/impl/trig_reduce.hpp>
#include <boost/decimal/detail/cmath/impl/cos_impl.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/cmath/atan.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#endif

namespace boost {
namespace decimal {
namespace detail {

namespace ellint_detail {

// The first kind gives F and K, and the second kind gives E
enum class ellint_kind : bool { first, second };

// 1 / a = r * 2^s with r in (1/2, 1], for a > 0, with one more Newton step for Precise
template <int N, bool Precise = false>
constexpr auto ellint_recip(const fx<N>& a, int& s) noexcept -> fx<N>
{
    s = fx_countl_zero(a) - 1;
    const auto v {fx_shift(a, s)};
    const auto r {fx_recip(v)};
    return Precise ? fx_mul(r, fx_sub(fx_int<N>(2U), fx_mul(v, r))) : r;
}

// sqrt(d) as fx_sqrt, with one more Newton step of 1 / sqrt for Precise
template <int N, bool Precise = false>
constexpr auto ellint_sqrt(const fx<N>& d) noexcept -> fx<N>
{
    if (!Precise)
    {
        return fx_sqrt(d);
    }
    const auto k {static_cast<unsigned>(countl_zero(d.w[N - 1]) / 2)};
    const auto m {k == 0U ? d : fx_shl(d, 2U * k)};
    return fx_shr(fx_mul(m, fx_rsqrt_step(m, fx_rsqrt(m))), k);
}

// k = sig * 10^-q < 1 in fixed point, and 1 - k = u * 10^-q where k > 0.99, else u = 0
template <int N>
struct ellint_modulus
{
    boost::int128::uint128_t u;
    fx<N> k;
    int q;
};

template <typename T, int N>
constexpr auto ellint_modulus_of(const T k) noexcept -> ellint_modulus<N>
{
    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(k, &e))};
    ellint_modulus<N> m {0U, fx_from_decimal<N>(sig, -e), -e};
    if (k > T {99, -2})
    {
        m.u = pow10(static_cast<boost::int128::uint128_t>(m.q)) - sig;
    }
    return m;
}

// k' = sqrt((1 - k) (1 + k)) to a few units: for k' below 1/16 (d2 below 2^-8), the exact 1 - k = u * 10^-q scales
// by 10^2j
template <int N, bool Precise = false>
constexpr auto ellint_kp(const ellint_modulus<N>& m) noexcept -> fx<N>
{
    const auto one {fx_one<N>()};
    const auto d2 {fx_mul(fx_sub(one, m.k), fx_add(one, m.k))};
    if (d2.w[N - 1] >= (UINT64_C(1) << 54U))
    {
        return ellint_sqrt<N, Precise>(d2);
    }
    const int j {(m.q - num_digits(m.u)) / 2};
    const auto p {pow10(static_cast<std::uint64_t>(j))};
    return fx_div_small(fx_sqrt(fx_mul(fx_from_decimal<N>(m.u, m.q - 2 * j), fx_add(one, m.k))), p);
}

// atan(a / b) for a, b >= 0 and max(a, b) in [1/2, 1)
template <int N>
constexpr auto ellint_atan(const fx<N>& a, const fx<N>& b) noexcept -> fx<N>
{
    return fx_less(b, a) ? fx_sub(fx_top<N>(fx_table::half_pi), atan_detail::atan_quotient(b, a)) : atan_detail::atan_quotient(a, b);
}

// The quadrant 0 to 3 of x + iy != 0, from x > 0 and y >= 0 counterclockwise
template <int N>
constexpr auto ellint_quadrant(const fx_signed<N>& x, const fx_signed<N>& y) noexcept -> unsigned
{
    const bool xz {fx_countl_zero(x.m) == 64 * N};
    const bool yz {fx_countl_zero(y.m) == 64 * N};
    const bool xp {!x.neg && !xz};
    const bool yp {!y.neg && !yz};
    if (xp && (yp || yz))
    {
        return 0U;
    }
    if (!xp && yp)
    {
        return 1U;
    }
    return (!xp && !xz && !yp) ? 2U : 3U;
}

// The AGM from a = 1 and b = k' with the phase of A&S 17.6.8 from psi, as z = x + iy turned by a x + i b y: F(psi) =
// phi_t / (2^t a), K = pi / (2a), E/K = 1 - sum 2^(n-1) c_n^2 and E(psi) = E/K F(psi) + sum c_n sin phi_n (A&S 17.4.38)
template <int N>
struct ellint_parts
{
    fx<N> u;                // phi_t / 2^t
    fx<N> r;                // 1 / a = r * 2^s
    fx<N> ratio;            // E/K
    fx_signed<N> zeta;      // sum c_n sin phi_n
    int s;
};

template <int N, bool Precise = false>
constexpr auto ellint_phase(const fx<N>& c, const fx<N>& sn, const fx<N>& kp, const fx<N>& k, const bool phase, const ellint_kind kind) noexcept -> ellint_parts<N>
{
    const bool second {kind == ellint_kind::second};
    auto a {fx_one<N>()};
    auto b {kp};
    fx_signed<N> x {c, false};
    fx_signed<N> y {sn, false};
    std::uint64_t q {};
    auto sum {fx_shr(fx_mul(k, k), 1U)};
    ellint_parts<N> out {};
    int t {};
    for (;;)
    {
        if (phase)
        {
            // (a x^2 - b y^2) + i (a + b) x y, scaled so that the larger part is in [1/2, 1)
            const auto nx {fx_signed_add(fx_signed<N> {fx_mul(a, fx_mul(x.m, x.m)), false}, fx_signed<N> {fx_mul(b, fx_mul(y.m, y.m)), true})};
            const fx_signed<N> ny {fx_mul(fx_add(a, b), fx_mul(x.m, y.m)), x.neg != y.neg};
            const auto lx {fx_countl_zero(nx.m)};
            const auto ly {fx_countl_zero(ny.m)};
            const int l {(lx < ly ? lx : ly) - 2};
            x = {fx_shift(nx.m, l), nx.neg};
            y = {fx_shift(ny.m, l), ny.neg};
            // phi_(n+1) = 2 phi_n + d with |d| < pi/2: the quadrant count is 2q - 1, 2q, 2q + 1 or 2q + 2
            const auto d {(ellint_quadrant(x, y) - 2U * static_cast<unsigned>(q)) & 3U};
            q = 2U * q + d - (d == 3U ? 4U : 0U);
        }

        // The truncation can leave b above a by a unit near the end
        const bool over {fx_less(a, b)};
        const auto cn {fx_shr(over ? fx_sub(b, a) : fx_sub(a, b), 1U)};
        if (second)
        {
            sum = fx_add(sum, fx_shift(fx_mul(cn, cn), t));
        }
        if (phase && second)
        {
            // c sin(phi) with sin(phi) = y / |z|
            const auto h {fx_add(fx_mul(x.m, x.m), fx_mul(y.m, y.m))};
            const int e {fx_countl_zero(h) / 2};
            const auto rs {fx_rsqrt(fx_shift(h, 2 * e))};
            out.zeta = fx_signed_add(out.zeta, fx_signed<N> {fx_shift(fx_mul(cn, fx_mul(y.m, rs)), e), y.neg});
        }
        const auto next {over ? fx_add(a, cn) : fx_sub(a, cn)};
        ++t;
        // c falls as its square: below a * 2^-32N, the rest is below 2^-64N
        if (fx_countl_zero(cn) - fx_countl_zero(next) >= 32 * N)
        {
            a = next;
            break;
        }
        b = ellint_sqrt<N, Precise>(fx_mul(a, b));
        a = next;
    }

    if (phase)
    {
        const bool even {(q & 1U) == 0U};
        const auto alpha {ellint_atan(even ? y.m : x.m, even ? x.m : y.m)};
        out.u = fx_add(fx_mul_shr(fx_top<N>(fx_table::half_pi), q, static_cast<unsigned>(t)), fx_shift(alpha, -t));
    }
    out.r = ellint_recip<N, Precise>(a, out.s);
    if (second)
    {
        out.ratio = fx_sub(fx_one<N>(), sum);
    }
    return out;
}

// F(psi) / 32 for k = 1 and psi = r or pi/2 - r with s = sin(psi), sinc of r and cr = cos(r)
template <int N, int W>
constexpr auto ellint_f_k_one(const trig::trig_arg<W>& r, const fx<N>& s, const fx<N>& sinc, const fx<N>& cr) noexcept -> fx<N>
{
    using sig_type = boost::int128::uint128_t;
    const bool odd {(r.n & 1U) != 0U};

    // ln(1 + s) - ln(c), with c = |r| sinc and ln(sig + tail / 2^64) = ln(sig) + tail / (2^64 sig) for odd
    auto f {log_detail::log_wide_of(fx_add(fx_one<N>(), s), 0, 0, log_detail::log_base::e)};
    if (odd)
    {
        f = fx_wide_add(f, log_detail::log_wide_of<N>(r.sig, num_digits(r.sig), -r.k, log_detail::log_base::e), true);
        f = fx_wide_add(f, log_detail::log_wide_of(fx_shl(sinc, 1U), -1, 0, log_detail::log_base::e), true);
        const auto tail {static_cast<sig_type>(r.tail) << 62U};
        const auto t {N == 1 ? sig_type {0U} : tail / (N == 2 ? r.sig : sig_type {r.sig.high})};
        fx<N> d {};
        for (int i {}; i < (N < 2 ? N : 2); ++i)
        {
            d.w[i] = i == 0 ? t.low : t.high;
        }
        f = fx_wide_add(f, fx_widen(d), true);
    }
    else
    {
        f = fx_wide_add(f, log_detail::log_wide_of(fx_shl(cr, 1U), -1, 0, log_detail::log_base::e), true);
    }
    fx<N> value {};
    fx_shr_words(fx_make_span(value.w), fx_make_span<N>(f.w), f.w[N], 5U);
    return value;
}

// The bits below 1 of a phi where F and E take two terms of their series: the first left out, below phi^4 / 24, is
// below 2^-48, 2^-100 and 2^-132
constexpr auto ellint_series_bits(const int p) noexcept -> int
{
    return p <= 7 ? 12 : p <= 16 ? 25 : 33;
}

// F or E for 0 < k <= 1 and phi > 0: phi = m pi + psi gives F(phi) = 2K/pi phi + (F(psi) - 2K/pi psi), where the first
// term keeps the precision of phi. For k = 1, F(psi) = ln((1 + s) / c) and E(psi) = s.
template <typename T>
constexpr auto ellint_fixed(const T k, const T phi, const ellint_kind kind, const bool neg, const bool taylor = false) noexcept -> T
{
    const bool second {kind == ellint_kind::second};
    constexpr int p {precision_v<T>};
    constexpr int words {fx_words(p)};
    using sig_type = boost::int128::uint128_t;

    const bool k_one {k == T {1}};
    const auto m {ellint_modulus_of<T, words>(k_one ? T {0} : k)};
    const auto k2 {k_one ? fx_one<words>() : fx_mul(m.k, m.k)};

    // phi = sig * 10^e, with 37 digits for F and 38 for E
    const int digits {second ? 38 : 37};
    const int scale {second ? 0 : 5};
    int ep {};
    const auto sp {static_cast<sig_type>(frexp10(phi, &ep))};
    const auto sig {sp * pow10(static_cast<sig_type>(digits - num_digits(sp)))};
    const int e {ep + num_digits(sp) - digits};
    // 10^-e for phi > pi/4, where e > -39
    const auto up {(e < 0 ? pow10(static_cast<sig_type>(e > -39 ? -e : 0)) : sig_type {1U}) << scale};

    const auto r {trig::trig_prepare(phi)};
    const auto rf {trig::decimal_r_with_tail<words>(r, 38 - trig::sig_digits(trig::trig_traits<T>::words))};
    const auto z {fx_mul(rf, rf)};
    const bool unreduced {r.n == 0U && phi < T {1}};
    const bool near {unreduced || (r.n == 1U && r.neg && phi < T {2})};
    const int lead {fx_countl_zero(rf)};

    sig_type digits_of {};
    if (unreduced && lead - 2 >= ellint_series_bits(p))
    {
        // F / phi = 1 + k^2 phi^2 / 6 and E / phi = 1 - k^2 phi^2 / 6
        const auto t {fx_div_small(fx_mul(k2, z), 6U)};
        const auto one {fx_one<words>()};
        digits_of = fx_scale(sig << scale, fx_shift(second ? fx_sub(one, t) : fx_add(one, t), -scale));
    }
    else
    {
        const auto sinc {trig::sinc_fx_of(z, taylor)};
        const auto cr {trig::cos_fx_of(z, taylor)};
        const auto sr {fx_mul(rf, sinc)};
        const bool odd {(r.n & 1U) != 0U};
        const auto c {odd ? sr : cr};
        const auto s {odd ? cr : sr};

        // F(psi) / 32 or E(psi), and 2K/pi / 32 or 2E/pi
        fx<words> value {};
        fx<words> g {};
        if (k_one)
        {
            if (!second)
            {
                if (!near)
                {
                    return neg ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
                }
                value = ellint_f_k_one(r, s, sinc, cr);
            }
            else
            {
                value = s;
                g = fx_recip(fx_top<words>(fx_table::half_pi));
            }
        }
        else
        {
            const auto v {ellint_phase(c, s, ellint_kp(m), m.k, true, kind)};
            const auto ur {fx_mul(v.u, v.r)};
            if (second)
            {
                value = fx_signed_add(fx_signed<words> {fx_shift(fx_mul(v.ratio, ur), v.s), false}, v.zeta).m;
                g = fx_shift(fx_mul(v.ratio, v.r), v.s);
            }
            else
            {
                value = fx_shift(ur, v.s - scale);
                g = fx_shift(v.r, v.s - scale);
            }
        }

        if (unreduced)
        {
            // value / phi, with phi in [1, 2) after the shift
            const int l {lead - 1};
            digits_of = fx_scale(sig << scale, fx_div(fx_shift(value, l), fx_shift(rf, l)));
        }
        else if (near)
        {
            digits_of = fx_scale(up, value);
        }
        else
        {
            const auto psi {odd ? fx_sub(fx_top<words>(fx_table::half_pi), rf) : rf};
            const bool psi_neg {odd ? !r.neg : r.neg};
            const auto gp {fx_mul(g, psi)};
            const bool below {fx_less(value, gp)};
            const auto b {below ? fx_sub(gp, value) : fx_sub(value, gp)};

            // a + b with the fractions of both, where b scales by 10^-e
            std::uint64_t fa {};
            std::uint64_t fb {};
            const auto a {fx_scale_frac(sig << scale, g, fa)};
            auto bd {fx_scale_frac(up, b, fb)};
            if (e > 0)
            {
                const auto t {e < 39 ? ((bd << 64U) | sig_type {fb}) / pow10(static_cast<sig_type>(e)) : sig_type {0U}};
                bd = t >> 64U;
                fb = t.low;
            }
            if (below != psi_neg)
            {
                digits_of = a - bd - (fa < fb ? 1U : 0U);
            }
            else
            {
                digits_of = a + bd + (fa + fb < fa ? 1U : 0U);
            }
        }
        // In three words the polynomials of sin and cos are good to a unit of the 38 digits: near a boundary the series decides
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            if (!taylor && fx_near_boundary(digits_of))
            {
                return ellint_fixed(k, phi, kind, neg, true);
            }
        }
    }

    // F > phi and E < phi
    if (second ? digits_of >= sig : digits_of < sig)
    {
        digits_of = second ? sig - 1U : sig;
    }
    return fx_round<T>(digits_of, e, neg);
}

// K(k) or E(k) for 0 <= k < 1 in N words as digits * 10^e, truncated, where the value / 2^scale is below 2
template <typename T, int N, bool Precise = false>
constexpr auto comp_ellint_digits(const T k, const ellint_kind kind, const int scale, const int e) noexcept -> boost::int128::uint128_t
{
    const bool second {kind == ellint_kind::second};
    const auto m {ellint_modulus_of<T, N>(k)};
    const auto v {ellint_phase<N, Precise>(fx_one<N>(), fx<N> {}, ellint_kp<N, Precise>(m), m.k, false, kind)};
    const auto f {fx_mul(fx_top<N>(fx_table::half_pi), second ? fx_mul(v.ratio, v.r) : v.r)};
    return fx_scale(pow10(static_cast<boost::int128::uint128_t>(-e)) << scale, fx_shift(f, v.s - scale));
}

// K(k) or E(k) for 0 <= k < 1, rounded once in the current mode
template <typename T>
constexpr auto comp_ellint_fixed(const T k, const ellint_kind kind) noexcept -> T
{
    const bool second {kind == ellint_kind::second};
    constexpr int words {fx_words(precision_v<T>)};

    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        // For k <= 0.99, a first try in one word: its error is below 70 units of 2^-62 for k <= 0.7, 100 for k <= 0.9,
        // else 220
        T r {};
        if (k <= T {99, -2} && fx_try_round_digits(comp_ellint_digits<T, 1, true>(k, kind, second ? 0 : 2, -37), -37, false, k <= T {7, -1} ? 70U : k <= T {9, -1} ? 100U : 220U, r))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // For k <= 0.99, where K < 4 and k' needs no u, a first try in two words: its error is below 120 units of the
        // 38 digits for k <= 0.7, else 558
        if (k <= T {99, -2})
        {
            const auto d2 {comp_ellint_digits<T, 2, true>(k, kind, second ? 0 : 2, -37)};
            if (!fx_near_boundary(d2, k <= T {7, -1} ? 120U : 558U))
            {
                return fx_round<T>(d2, -37, false);
            }
        }
    }
    // K stays below 48 and E below 2: K / 32 keeps K in fx, and K takes one digit more before the point
    const int e {second ? -37 : -36};
    return fx_round<T>(comp_ellint_digits<T, words>(k, kind, second ? 0 : 5, e), e, false);
}

} //namespace ellint_detail

} //namespace detail
} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_IMPL_ELLINT_IMPL_HPP
