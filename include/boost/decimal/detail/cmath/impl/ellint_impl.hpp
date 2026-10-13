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

// 1 / a = r * 2^s with r in (1/2, 1], for a > 0, with one more Newton step for Precise in two words: one word takes
// one division
template <int N, bool Precise = false>
constexpr auto ellint_recip(const fx<N>& a, int& s) noexcept -> fx<N>
{
    s = fx_countl_zero(a) - 1;
    const auto v {fx_shift(a, s)};
    const auto r {fx_recip(v)};
    return Precise && N > 1 ? fx_mul(r, fx_sub(fx_int<N>(2U), fx_mul(v, r))) : r;
}

// sqrt(d) as fx_sqrt, with one more Newton step of 1 / sqrt for Precise in two words, where in one word its
// truncation is larger than the error of fx_sqrt
template <int N, bool Precise = false>
constexpr auto ellint_sqrt(const fx<N>& d) noexcept -> fx<N>
{
    if (!Precise || N == 1)
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

// The class of k = sig * 10^e with p digits, or k = 0: up to 0.7, 0.9, 0.97 or 0.99, then below 1, then 1
constexpr auto ellint_band(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> unsigned
{
    const auto u {pow10(static_cast<boost::int128::uint128_t>(p - 2))};
    if (e + p != 0)
    {
        return e + p > 0 && sig != 0U ? 5U : 0U;
    }
    return sig <= u * 70U ? 0U : sig <= u * 90U ? 1U : sig <= u * 97U ? 2U : sig <= u * 99U ? 3U : 4U;
}

// For k in the band 5 of ellint_band: true for k > 1, false for k = 1
constexpr auto ellint_above_one(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> bool
{
    return e + p != 1 || sig != pow10(static_cast<boost::int128::uint128_t>(p - 1));
}

template <int N>
constexpr auto ellint_modulus_of(const boost::int128::uint128_t& sig, const int e, const unsigned band) noexcept -> ellint_modulus<N>
{
    ellint_modulus<N> m {0U, fx_from_decimal<N>(sig, -e), -e};
    if (band == 4U)
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
    fx<N> cl {};
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
        // For Precise, the shifts come before the products, which keeps the truncation at one unit
        if (second)
        {
            sum = fx_add(sum, Precise ? fx_mul(fx_shift(cn, t), cn) : fx_shift(fx_mul(cn, cn), t));
        }
        if (phase && second)
        {
            // c sin(phi) with sin(phi) = y / |z|
            const auto h {fx_add(fx_mul(x.m, x.m), fx_mul(y.m, y.m))};
            const int e {fx_countl_zero(h) / 2};
            const auto rs {fx_rsqrt(fx_shift(h, 2 * e))};
            const auto term {Precise ? fx_mul(cn, fx_mul(y.m, fx_shift(rs, e))) : fx_shift(fx_mul(cn, fx_mul(y.m, rs)), e)};
            out.zeta = fx_signed_add(out.zeta, fx_signed<N> {term, y.neg});
        }
        const auto next {over ? fx_add(a, cn) : fx_sub(a, cn)};
        ++t;
        // c falls as its square: below a * 2^-16N, the last step below takes c' = c^2 / 4a without a square root
        if (fx_countl_zero(cn) - fx_countl_zero(next) >= 16 * N)
        {
            a = next;
            cl = cn;
            break;
        }
        b = ellint_sqrt<N, Precise>(fx_mul(a, b));
        a = next;
    }

    // The last step: b = a - 2c' gives M = a - c', 2^t c'^2 in the sum, and phi gains 2 phi - c' / a sin(2 phi) with
    // z turned by a z - 2i c' y. The terms left out are below 2^-64N. With c below 2^(1 - 16N), cr = c' 2^2g < a / 4.
    constexpr int g {16 * N - 1};
    const auto r {ellint_recip<N, Precise>(a, out.s)};
    const auto cg {fx_shift(cl, g)};
    const auto cr {fx_shift(fx_mul(fx_mul(cg, cg), r), out.s - 2)};
    out.r = fx_add(r, fx_shift(fx_mul(r, fx_mul(cr, r)), out.s - 2 * g));
    if (second)
    {
        out.ratio = fx_sub(fx_one<N>(), fx_add(sum, fx_shift(fx_mul(cr, cr), t - 4 * g)));
    }
    if (phase)
    {
        const bool even {(q & 1U) == 0U};
        const auto alpha {ellint_atan(even ? y.m : x.m, even ? x.m : y.m)};
        const auto h {fx_add(fx_mul(x.m, x.m), fx_mul(y.m, y.m))};
        const int e {fx_countl_zero(h) - 1};
        const fx_signed<N> s2 {fx_div(fx_shift(fx_shl(fx_mul(x.m, y.m), 1U), e), fx_shift(h, e)), x.neg != y.neg};
        const auto d {fx_shift(fx_mul(fx_mul(cr, r), s2.m), out.s - t - 1 - 2 * g)};
        const auto u {fx_add(fx_mul_shr(fx_top<N>(fx_table::half_pi), q, static_cast<unsigned>(t)), fx_shift(alpha, -t))};
        out.u = s2.neg ? fx_add(u, d) : fx_sub(u, d);
        if (second)
        {
            out.zeta = fx_signed_add(out.zeta, fx_signed<N> {fx_shift(fx_mul(cr, s2.m), -2 * g), s2.neg});
        }
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
        f = fx_wide_add(f, log_detail::log_wide_of<N>(r.sig, -r.k, log_detail::log_base::e), true);
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

// sin(r) / r and cos(r) for z = r^2, with the polynomials of the next wider type in one or two words for top
template <int N>
constexpr auto ellint_sinc(const fx<N>& z, const bool top, const bool taylor) noexcept -> std::enable_if_t<(N < 3), fx<N>>
{
    return top ? trig::sinc_fx(z, trig::sin_poly_top(z)) : trig::sinc_fx_of(z, taylor);
}

template <int N>
constexpr auto ellint_sinc(const fx<N>& z, const bool, const bool taylor) noexcept -> std::enable_if_t<(N >= 3), fx<N>>
{
    return trig::sinc_fx_of(z, taylor);
}

template <int N>
constexpr auto ellint_cos(const fx<N>& z, const bool top, const bool taylor) noexcept -> std::enable_if_t<(N < 3), fx<N>>
{
    return top ? trig::cos_fx(z, trig::cos_poly_top(z)) : trig::cos_fx_of(z, taylor);
}

template <int N>
constexpr auto ellint_cos(const fx<N>& z, const bool, const bool taylor) noexcept -> std::enable_if_t<(N >= 3), fx<N>>
{
    return trig::cos_fx_of(z, taylor);
}

// k and phi > 0 as sig * 10^e with p digits, and the class of k from ellint_band
template <typename T>
struct ellint_args
{
    T phi;
    boost::int128::uint128_t ks;
    boost::int128::uint128_t ps;
    int ke;
    int pe;
    unsigned band;
};

// |r| of phi = n pi/2 + r in N words, with n and the sign of r in out: a first try takes the reduction in N words where
// it can, which sets only these two fields
template <typename T, int N>
constexpr auto ellint_reduce(const T phi, trig::trig_arg<trig::trig_traits<T>::words>& out) noexcept -> fx<N>
{
    BOOST_DECIMAL_IF_CONSTEXPR (N < fx_words(precision_v<T>))
    {
        trig::trig_arg<N> r {};
        if (trig::trig_prepare_try(phi, r))
        {
            out.n = r.n;
            out.neg = r.neg;
            return trig::fixed_r(r);
        }
    }
    out = trig::trig_prepare(phi);
    return trig::decimal_r_with_tail<N>(out, 38 - trig::sig_digits(trig::trig_traits<T>::words));
}

// The digits of E(1, phi) next to an odd multiple of pi/2, from d near them and with up units in 1: E is that multiple
// m plus or minus 1 - cos(r), which is r^2 (1/2 - z q) in its own precision, for the polynomial q of cos(r) and z = r^2
template <int N>
constexpr auto ellint_e_one_digits(const boost::int128::uint128_t& d, const boost::int128::uint128_t& up, const fx<N>& rf, const fx<N>& z, const int lead, const bool below) noexcept -> boost::int128::uint128_t
{
    const auto m {(d + up / 2U) / up * up};
    fx<N> half {};
    half.w[N - 1] = UINT64_C(1) << 61U;
    const int l {lead - 1};
    const auto w {fx_shift(rf, l)};
    const auto rest {fx_scale(up, fx_mul(fx_mul(w, w), fx_sub(half, fx_mul(z, trig::cos_poly(z)))))};
    const auto q {2 * l < 128 ? rest >> static_cast<unsigned>(2 * l) : boost::int128::uint128_t {}};
    return below ? m - q - 1U : m + q;
}

// F or E for 0 < k <= 1 and phi > 0 in N words: phi = m pi + psi gives F(phi) = 2K/pi phi + (F(psi) - 2K/pi psi), where
// the first term keeps the precision of phi. For k = 1, F(psi) = ln((1 + s) / c) and E(psi) = s.
// In fewer words, a first try for k <= 0.99 and phi >= 1/2: false where it cannot round, and out is not set
template <typename T, int N>
constexpr auto ellint_words(const ellint_args<T>& x, const ellint_kind kind, const bool neg, T& out, const bool taylor = false) noexcept -> bool
{
    const bool second {kind == ellint_kind::second};
    constexpr int p {precision_v<T>};
    constexpr int words {N};
    constexpr bool first_try {N < fx_words(p)};
    using sig_type = boost::int128::uint128_t;

    if (first_try && (x.band > 3U || !decimal_at_least(x.ps, x.pe, p, 50U, 0)))
    {
        return false;
    }

    const bool k_one {x.band == 5U};
    const auto m {ellint_modulus_of<words>(k_one ? sig_type {0U} : x.ks, k_one ? 0 : x.ke, x.band)};
    const auto k2 {k_one ? fx_one<words>() : fx_mul(m.k, m.k)};

    // phi = sig * 10^e, with 37 digits for F and 38 for E, and F / 32 but in the first try F < 4
    const int digits {second ? 38 : 37};
    const int scale {second || first_try ? 0 : 5};
    const auto sig {x.ps * pow10(static_cast<sig_type>(digits - p))};
    const int e {x.pe + p - digits};
    // 10^-e for phi > pi/4, where e > -39
    const auto up {(e < 0 ? pow10(static_cast<sig_type>(e > -39 ? -e : 0)) : sig_type {1U}) << scale};

    trig::trig_arg<trig::trig_traits<T>::words> r {};
    const auto rf {ellint_reduce<T, words>(x.phi, r)};
    const auto z {fx_mul(rf, rf)};
    const bool unreduced {r.n == 0U && x.pe + p <= 0};
    const bool near {unreduced || (r.n == 1U && r.neg && !decimal_at_least(x.ps, x.pe, p, 20U, 1))};
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
        const auto sinc {ellint_sinc(z, first_try, taylor)};
        const auto cr {ellint_cos(z, first_try, taylor)};
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
                    out = neg ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
                    return true;
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
            const auto v {ellint_phase<words, first_try>(c, s, ellint_kp<words, first_try>(m), m.k, true, kind)};
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
            if (k_one && second && odd && z.w[N - 1] == 0U && e < 0)
            {
                digits_of = ellint_e_one_digits(digits_of, up, rf, z, lead, r.neg);
            }
            else if (!taylor && fx_near_boundary(digits_of))
            {
                return ellint_words<T, N>(x, kind, neg, out, true);
            }
        }
    }

    // F > phi and E < phi
    if (second ? digits_of >= sig : digits_of < sig)
    {
        if (first_try)
        {
            return false;
        }
        digits_of = second ? sig - 1U : sig;
    }
    BOOST_DECIMAL_IF_CONSTEXPR (first_try && words == 2)
    {
        // digits_of has an error below R digits_of 2^-126 + 2 units: R is 175, 210 or 265 for F and 290, 310 or 350
        // for E, with k up to 0.7, 0.9 or 0.99
        const std::uint64_t err {x.band == 0U ? (second ? 290U : 175U) : x.band == 1U ? (second ? 310U : 210U) : (second ? 350U : 265U)};
        if (fx_near_boundary(digits_of, ((digits_of >> 64U) * err >> 62U).low + 2U))
        {
            return false;
        }
        out = fx_round<T>(digits_of, e, neg);
        return true;
    }
    BOOST_DECIMAL_IF_CONSTEXPR (first_try)
    {
        // 19 digits q with an error below R q 2^-62 + 2 units, with R for F and E by k up to 0.7, 0.9, 0.97 or 0.99, for
        // phi below pi/2 and above
        constexpr std::uint8_t bound[2][4][2] {{{40U, 62U}, {50U, 66U}, {53U, 66U}, {66U, 77U}}, {{20U, 29U}, {29U, 46U}, {38U, 62U}, {63U, 77U}}};
        const std::uint64_t err {bound[near ? 0 : 1][x.band][second ? 1 : 0]};
        const int drop {num_digits(digits_of) - 19};
        const auto q {impl::divmod_pow10_dispatch(digits_of, drop, pow10(static_cast<sig_type>(drop))).quotient.low};
        return fx_try_round(q, e + drop, neg, (static_cast<sig_type>(q) * err >> 62U).low + 2U, out);
    }
    out = fx_round<T>(digits_of, e, neg);
    return true;
}

// The digits of |k| and |phi| for ellint_fixed
template <typename T>
constexpr auto ellint_args_of(const T k, const T phi) noexcept -> ellint_args<T>
{
    ellint_args<T> x {fabs(phi), 0U, 0U, 0, 0, 0U};
    x.ks = static_cast<boost::int128::uint128_t>(frexp10(k, &x.ke));
    x.ps = static_cast<boost::int128::uint128_t>(frexp10(phi, &x.pe));
    x.band = ellint_band(x.ks, x.ke, precision_v<T>);
    return x;
}

// ellint_words, with a first try in one word fewer for decimal64_t and decimal128_t
template <typename T>
constexpr auto ellint_fixed(const ellint_args<T> x, const ellint_kind kind, const bool neg) noexcept -> T
{
    constexpr int words {fx_words(precision_v<T>)};
    T r {};
    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        if (ellint_words<T, (words > 1 ? words - 1 : 1)>(x, kind, neg, r))
        {
            return r;
        }
    }
    ellint_words<T, words>(x, kind, neg, r);
    return r;
}

template <bool b>
struct ellint_fit_table_imp
{
    // K / 2 and E / 2 in u = 8m - i on m = k^2 in [i / 8, (i + 1) / 8) for i < 6, and P_K / 2 and P_E / 2 in u = 4c on
    // c = 1 - m below 1/4: polynomials as two's complement words times 2^62, by interval from the constant term
    static constexpr std::uint64_t k[96] = {
        UINT64_C(0x3243F6A8885A308D),
        UINT64_C(0x01921FB54442D185),
        UINT64_C(0x001C463ABECCB2B3),
        UINT64_C(0x000274518B3AA939),
        UINT64_C(0x00003C21CDD31237),
        UINT64_C(0x000006169F66AC89),
        UINT64_C(0x000000A3B59FE57D),
        UINT64_C(0x00000011A511C3FC),
        UINT64_C(0x00000001F0334732),
        UINT64_C(0x00000000376D8A59),
        UINT64_C(0x000000000620D16B),
        UINT64_C(0x0000000000CD38A9),
        UINT64_C(0x000000000009233F),
        UINT64_C(0x000000000005AB90),
        UINT64_C(0x33F513DA0D2D70A4),
        UINT64_C(0x01D31C7B385AA857),
        UINT64_C(0x00255425F139F50F),
        UINT64_C(0x0003B170EB760C0D),
        UINT64_C(0x0000673D47882025),
        UINT64_C(0x00000BEEF230A3BD),
        UINT64_C(0x0000016E6EE3EED0),
        UINT64_C(0x0000002D1C7BAB46),
        UINT64_C(0x00000005A8EAF87E),
        UINT64_C(0x00000000B953365C),
        UINT64_C(0x0000000016F6DFD4),
        UINT64_C(0x0000000003C14302),
        UINT64_C(0x000000000013798F),
        UINT64_C(0x000000000024AE9D),
        UINT64_C(0x35F1AABA647EF1EA),
        UINT64_C(0x022ABBC0F3FB8ED9),
        UINT64_C(0x003365271DB762EE),
        UINT64_C(0x0005EA0E39E920F9),
        UINT64_C(0x0000C08B07B044D8),
        UINT64_C(0x000019EF629A5412),
        UINT64_C(0x000003A05849D9C1),
        UINT64_C(0x0000008542852410),
        UINT64_C(0x000000137C1D2DE3),
        UINT64_C(0x00000002ECE1EEE4),
        UINT64_C(0x00000000670BE4E9),
        UINT64_C(0x00000000178F8200),
        UINT64_C(0xFFFFFFFFFF4D4532),
        UINT64_C(0x000000000144E215),
        UINT64_C(0x385694679D61465C),
        UINT64_C(0x02A6E23512E49BF8),
        UINT64_C(0x004AEDE3B9A6E711),
        UINT64_C(0x000A4FCC5F9D4DD7),
        UINT64_C(0x0001920ED8F8D57F),
        UINT64_C(0x000040E7D4081A9C),
        UINT64_C(0x00000AE181286E11),
        UINT64_C(0x000001DF7B7DBBD5),
        UINT64_C(0x000000542833FEA1),
        UINT64_C(0x0000000F04686815),
        UINT64_C(0x00000002AB79EA2D),
        UINT64_C(0x0000000089B3536A),
        UINT64_C(0x000000000C532AA5),
        UINT64_C(0x000000000B33C962),
        UINT64_C(0xFFFFFFFFFDF9A7C0),
        UINT64_C(0x0000000000C12511),
        UINT64_C(0x3B54946AEB3C9661),
        UINT64_C(0x03638BD3AFBBD129),
        UINT64_C(0x0076A928D5D67919),
        UINT64_C(0x00145546F61E6A34),
        UINT64_C(0x0003DCD6FEF5B160),
        UINT64_C(0x0000C743B76FE387),
        UINT64_C(0x000029B77837ED51),
        UINT64_C(0x000008F84A2A8CAC),
        UINT64_C(0x000001F79090D290),
        UINT64_C(0x000000702C967F7B),
        UINT64_C(0x0000001923A7D7A4),
        UINT64_C(0x00000005FD410575),
        UINT64_C(0x00000000FA93E91B),
        UINT64_C(0x00000000A9E3DFCA),
        UINT64_C(0xFFFFFFFFCBA13A6F),
        UINT64_C(0x0000000029CE35AF),
        UINT64_C(0xFFFFFFFFF43D450E),
        UINT64_C(0x0000000002753971),
        UINT64_C(0x3F47F8014E5F4062),
        UINT64_C(0x04A283BE66731D2F),
        UINT64_C(0x00D619C870043311),
        UINT64_C(0x0030B1AC9F01ED8C),
        UINT64_C(0x000C4DAF5533A0B7),
        UINT64_C(0x00034CFB3EB2789D),
        UINT64_C(0x0000EB9F639D72CE),
        UINT64_C(0x00004380461C4288),
        UINT64_C(0x000013B9950E6A43),
        UINT64_C(0x000005DC1886ECDC),
        UINT64_C(0x000001BCEB5CB1D3),
        UINT64_C(0x00000095999B3747),
        UINT64_C(0x0000001025625B84),
        UINT64_C(0x00000033966DC978),
        UINT64_C(0xFFFFFFD70548D022),
        UINT64_C(0x00000028A5E7FE13),
        UINT64_C(0xFFFFFFE7092157AC),
        UINT64_C(0x0000000B7FCD7FE3),
        UINT64_C(0xFFFFFFFCCF6CC928),
        UINT64_C(0x0000000076F5C963)
    };
    static constexpr std::uint64_t e[92] = {
        UINT64_C(0x3243F6A8885A308D),
        UINT64_C(0xFE6DE04ABBBD2E78),
        UINT64_C(0xFFF69341C0666FD5),
        UINT64_C(0xFFFF8256175A9D57),
        UINT64_C(0xFFFFF768E29936DB),
        UINT64_C(0xFFFFFF52D1D5454B),
        UINT64_C(0xFFFFFFF11E12F989),
        UINT64_C(0xFFFFFFFEA470A428),
        UINT64_C(0xFFFFFFFFDF102E60),
        UINT64_C(0xFFFFFFFFFC98B8E9),
        UINT64_C(0xFFFFFFFFFFC3BB81),
        UINT64_C(0xFFFFFFFFFFEF4056),
        UINT64_C(0x30A7E3366E266927),
        UINT64_C(0xFE5967AE307C7C41),
        UINT64_C(0xFFF4DEF5A5CA36DA),
        UINT64_C(0xFFFF57D42F913BB1),
        UINT64_C(0xFFFFF2ED44D679F4),
        UINT64_C(0xFFFFFED3C8A905AF),
        UINT64_C(0xFFFFFFE29552B94D),
        UINT64_C(0xFFFFFFFCF076F028),
        UINT64_C(0xFFFFFFFFAADA7366),
        UINT64_C(0xFFFFFFFFF66FC178),
        UINT64_C(0xFFFFFFFFFEE602AF),
        UINT64_C(0xFFFFFFFFFFDE90B4),
        UINT64_C(0xFFFFFFFFFFFC3DC2),
        UINT64_C(0xFFFFFFFFFFFF601B),
        UINT64_C(0x0000000000000025),
        UINT64_C(0xFFFFFFFFFFFFFA8F),
        UINT64_C(0x2EF5734EA751E1DD),
        UINT64_C(0xFE40F22510B4BBFD),
        UINT64_C(0xFFF28A433F69F6B7),
        UINT64_C(0xFFFF15014300E2F1),
        UINT64_C(0xFFFFEACEB7713FE4),
        UINT64_C(0xFFFFFDCA4C35267D),
        UINT64_C(0xFFFFFFBF7F548BB2),
        UINT64_C(0xFFFFFFF82F560E5D),
        UINT64_C(0xFFFFFFFF02D6EC3F),
        UINT64_C(0xFFFFFFFFDEA9F901),
        UINT64_C(0xFFFFFFFFFBC434A9),
        UINT64_C(0xFFFFFFFFFF3743AE),
        UINT64_C(0x0000000000004108),
        UINT64_C(0xFFFFFFFFFFF6FE48),
        UINT64_C(0x2D27ED07C93614D9),
        UINT64_C(0xFE22E41AB1F8CD12),
        UINT64_C(0xFFEF2F795A42E2C6),
        UINT64_C(0xFFFEA3C2ADC0D400),
        UINT64_C(0xFFFFDA8EED27BA41),
        UINT64_C(0xFFFFFB559F43879D),
        UINT64_C(0xFFFFFF5D0B5AF941),
        UINT64_C(0xFFFFFFE85AE4E416),
        UINT64_C(0xFFFFFFFC6B7C3C39),
        UINT64_C(0xFFFFFFFF6D7CCFED),
        UINT64_C(0xFFFFFFFFEBBA4F01),
        UINT64_C(0xFFFFFFFFF9E7CD0F),
        UINT64_C(0x00000000007A95ED),
        UINT64_C(0xFFFFFFFFFF9A19CF),
        UINT64_C(0x2B387984348D8FD6),
        UINT64_C(0xFDFC7CA3292A1F2E),
        UINT64_C(0xFFE9FF789271A121),
        UINT64_C(0xFFFDCE59E4C8D717),
        UINT64_C(0xFFFFB517BC900353),
        UINT64_C(0xFFFFF4634569D564),
        UINT64_C(0xFFFFFE06A498002C),
        UINT64_C(0xFFFFFFA48F59E50D),
        UINT64_C(0xFFFFFFEEB33EF7DA),
        UINT64_C(0xFFFFFFFC988C51A3),
        UINT64_C(0xFFFFFFFF5895EA56),
        UINT64_C(0xFFFFFFFFD30BF5B3),
        UINT64_C(0x0000000000ECD6F6),
        UINT64_C(0xFFFFFFFFF915F7DD),
        UINT64_C(0x0000000001C92287),
        UINT64_C(0xFFFFFFFFFF7C6585),
        UINT64_C(0x291C6B0A7D936593),
        UINT64_C(0xFDC871E7518536EB),
        UINT64_C(0xFFE10D515066C899),
        UINT64_C(0xFFFBF546854AEA05),
        UINT64_C(0xFFFF49B2DE90364B),
        UINT64_C(0xFFFFDA874AA7CA3F),
        UINT64_C(0xFFFFF789A6F9E058),
        UINT64_C(0xFFFFFDF6B84B6E5C),
        UINT64_C(0xFFFFFF7CB54FC6DE),
        UINT64_C(0xFFFFFFDDC25A9DD9),
        UINT64_C(0xFFFFFFF6E01CD0E0),
        UINT64_C(0xFFFFFFFD62AD41BD),
        UINT64_C(0xFFFFFFFF82617A13),
        UINT64_C(0xFFFFFFFF7CBFC944),
        UINT64_C(0x000000004F5C8FC6),
        UINT64_C(0xFFFFFFFFA95C1894),
        UINT64_C(0x0000000033FD0093),
        UINT64_C(0xFFFFFFFFE7A704F4),
        UINT64_C(0x0000000006BEA3A1),
        UINT64_C(0xFFFFFFFFFEFECC5C)
    };
    static constexpr std::uint64_t pk[16] = {
        UINT64_C(0x2C5C85FDF473DE6B),
        UINT64_C(0x00C5C85FDF473DE7),
        UINT64_C(0x000FD02D7B660486),
        UINT64_C(0x0001E979AC88D8A9),
        UINT64_C(0x000047D04A05C5F0),
        UINT64_C(0x00000BC937F7919E),
        UINT64_C(0x00000214C330A5F2),
        UINT64_C(0x000000630B866336),
        UINT64_C(0x000000131F42637D),
        UINT64_C(0x00000003D2DF1902),
        UINT64_C(0x00000000BDF85580),
        UINT64_C(0x000000003412B5DA),
        UINT64_C(0xFFFFFFFFFE9F46CE),
        UINT64_C(0x000000000837BF62),
        UINT64_C(0xFFFFFFFFFDDC4A64),
        UINT64_C(0x00000000009CCF01)
    };
    static constexpr std::uint64_t pe[16] = {
        UINT64_C(0x2000000000000000),
        UINT64_C(0x038B90BFBE8E7BCE),
        UINT64_C(0x001D1591F9DD5B6A),
        UINT64_C(0x0002CB5ECF0A9CBC),
        UINT64_C(0x00005E929DBD7D7E),
        UINT64_C(0x00000EA07712D831),
        UINT64_C(0x0000027C520783FA),
        UINT64_C(0x000000732BF0807A),
        UINT64_C(0x00000015CD14AC52),
        UINT64_C(0x000000044ADB3497),
        UINT64_C(0x00000000D33099E4),
        UINT64_C(0x0000000038C6B5B1),
        UINT64_C(0xFFFFFFFFFECF1200),
        UINT64_C(0x0000000008BF56E5),
        UINT64_C(0xFFFFFFFFFDBD52F6),
        UINT64_C(0x0000000000A68663)
    };
    // 2 / pi
    static constexpr std::uint64_t two_over_pi {UINT64_C(0x28BE60DB9391054A)};
};

template <bool b>
constexpr std::uint64_t ellint_fit_table_imp<b>::k[96];

template <bool b>
constexpr std::uint64_t ellint_fit_table_imp<b>::e[92];

template <bool b>
constexpr std::uint64_t ellint_fit_table_imp<b>::pk[16];

template <bool b>
constexpr std::uint64_t ellint_fit_table_imp<b>::pe[16];

template <bool b>
constexpr std::uint64_t ellint_fit_table_imp<b>::two_over_pi;

using ellint_fit_table = ellint_fit_table_imp<true>;

// The digits of K or E from one word, 18 to 20 of them, and their error bound in units of 2^-62
struct ellint_fit_digits
{
    boost::int128::uint128_t digits;
    std::uint64_t err;
    int exp10;
};

// K(k) or E(k) for k = sig * 10^-q < 1 in one word: their polynomials for m = k^2 below 3/4, else K = P_K(c) + K(c)
// ln(1 / c) / pi and E = P_E(c) + (K(c) - E(c)) ln(1 / c) / pi with the exact c = 1 - m
constexpr auto ellint_fit(const boost::int128::uint128_t& sig, const int q, const ellint_kind kind) noexcept -> ellint_fit_digits
{
    using sig_type = boost::int128::uint128_t;
    using table = ellint_fit_table;
    constexpr int kf[7] {0, 14, 28, 42, 58, 76, 96};
    constexpr int ef[7] {0, 12, 28, 42, 56, 72, 92};

    const bool second {kind == ellint_kind::second};
    const auto s2 {sig * sig};
    const auto m {fx_from_decimal<1>(s2, 2 * q)};
    if (m.w[0] < (UINT64_C(3) << 60U))
    {
        // The value / 2 has an error below 4 units, and m one unit: K or E is below 8 + 2 units off, plus a margin
        const auto i {static_cast<int>(m.w[0] >> 59U)};
        const std::uint64_t u {m.w[0] << 5U};
        fx<1> v {};
        v.w[0] = second ? fx_poly_signed(table::e, ef[i], ef[i + 1], u) : fx_poly_signed(table::k, kf[i], kf[i + 1], u);
        return {fx_scale(pow10(static_cast<sig_type>(18)) * 2U, v), 16U, -18};
    }

    // Here k > 0.86 has p digits, so 10^(2q) fits
    const auto cs {pow10(static_cast<sig_type>(2 * q)) - s2};
    const auto c {fx_from_decimal<1>(cs, 2 * q)};
    const auto ic {static_cast<int>(c.w[0] >> 59U)};
    const std::uint64_t uc {c.w[0] << 5U};
    const std::uint64_t kc {fx_poly_signed(table::k, kf[ic], kf[ic + 1], uc)};
    fx<1> qv {};
    qv.w[0] = second ? kc - fx_poly_signed(table::e, ef[ic], ef[ic + 1], uc) : kc;
    fx<1> two_pi {};
    two_pi.w[0] = table::two_over_pi;
    fx<1> pv {};
    pv.w[0] = fx_poly_signed(second ? table::pe : table::pk, 0, 16, c.w[0] << 4U) << 1U;
    // ln(1 / c) = 2q ln(10) - ln(cs) below 74
    const auto l {fx_wide_add(fx_wide<1> {}, log_detail::log_wide_of<1>(cs, -2 * q, log_detail::log_base::e), true)};
    const auto w {fx_wide_add(fx_widen(pv), fx_wide_mul(l, fx_mul(qv, two_pi)))};
    // P is 8 units off and Q 6, which ln(1 / c) multiplies, and ln(1 / c) is 4 units off
    const std::uint64_t whole {(w.w[1] << 2U) | (w.w[0] >> 62U)};
    const std::uint64_t lw {(l.w[1] << 2U) | (l.w[0] >> 62U)};
    const int s {whole >= 10U ? 17 : 18};
    fx<1> frac {};
    frac.w[0] = w.w[0] & ((UINT64_C(1) << 62U) - 1U);
    const auto p10 {pow10(static_cast<sig_type>(s))};
    return {p10 * whole + fx_scale(p10, frac), 16U + 6U * (lw + 1U), -s};
}

// K(k) or E(k) for 0 <= k < 1 in N words as digits * 10^e, truncated, where the value / 2^scale is below 2
template <int N, bool Precise = false>
constexpr auto comp_ellint_digits(const boost::int128::uint128_t& sig, const int q, const unsigned band, const ellint_kind kind, const int scale, const int e) noexcept -> boost::int128::uint128_t
{
    const bool second {kind == ellint_kind::second};
    const auto m {ellint_modulus_of<N>(sig, q, band)};
    const auto v {ellint_phase<N, Precise>(fx_one<N>(), fx<N> {}, ellint_kp<N, Precise>(m), m.k, false, kind)};
    const auto f {fx_mul(fx_top<N>(fx_table::half_pi), second ? fx_mul(v.ratio, v.r) : v.r)};
    return fx_scale(pow10(static_cast<boost::int128::uint128_t>(-e)) << scale, fx_shift(f, v.s - scale));
}

// K(k) or E(k) for 0 <= k < 1, rounded once in the current mode
template <typename T>
constexpr auto comp_ellint_fixed(const boost::int128::uint128_t& sig, const int e, const unsigned band, const ellint_kind kind) noexcept -> T
{
    const bool second {kind == ellint_kind::second};
    constexpr int words {fx_words(precision_v<T>)};

    BOOST_DECIMAL_IF_CONSTEXPR (words <= 2)
    {
        // A first try in one word from the polynomials of K and E
        const auto d {ellint_fit(sig, -e, kind)};
        T r {};
        if (fx_try_round_word(d.digits, d.exp10, false, d.err, r))
        {
            return r;
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // For k <= 0.99, where K < 4 and k' needs no u, a first try in two words: its error is below 15, 22 or 50 units
        // of the 38 digits for K and 13, 17 or 23 for E, with k up to 0.7, 0.9 or 0.99
        if (band <= 3U)
        {
            const auto d2 {comp_ellint_digits<2, true>(sig, e, band, kind, second ? 0 : 2, -37)};
            if (!fx_near_boundary(d2, band == 0U ? (second ? 13U : 15U) : band == 1U ? (second ? 17U : 22U) : (second ? 23U : 50U)))
            {
                return fx_round<T>(d2, -37, false);
            }
        }
    }
    // K stays below 48 and E below 2: K / 32 keeps K in fx, and K takes one digit more before the point
    const int ed {second ? -37 : -36};
    return fx_round<T>(comp_ellint_digits<words>(sig, e, band, kind, second ? 0 : 5, ed), ed, false);
}

} //namespace ellint_detail

} //namespace detail
} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_IMPL_ELLINT_IMPL_HPP
