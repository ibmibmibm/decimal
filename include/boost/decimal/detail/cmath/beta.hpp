// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/impl/gamma_impl.hpp>

namespace boost {
namespace decimal {

namespace detail {

namespace beta_detail {

using gamma_detail::gamma_float;
using gamma_detail::gamma_parts;

// x + y < 1 for x = sx * 10^ex and y = sy * 10^ey with p digits, as v * 2^n2 * 10^q with q = max(ex, ey) + p: the
// halves of x / 10^q and y / 10^q add up below 1, and the smaller loses only its digits below 2^-64N of the sum
template <int N>
constexpr auto beta_small_sum(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const int p) noexcept -> gamma_float<N>
{
    const int q {(ex > ey ? ex : ey) + p};
    const auto h {fx_add(fx_shr(fx_from_decimal<N>(sx, q - ex), 1U), fx_shr(fx_from_decimal<N>(sy, q - ey), 1U))};
    unsigned s {};
    const auto v {fx_normalize(h, s)};
    return {v, 1 - static_cast<int>(s), q};
}

// acc * f, or false where it does not fit in 256 bits
constexpr auto beta_mul(u256& acc, const boost::int128::uint128_t& f) noexcept -> bool
{
    std::uint64_t r[6] {};
    const std::uint64_t fw[2] {f.low, f.high};
    for (int j {}; j < 2; ++j)
    {
        std::uint64_t carry {};
        for (int i {}; i < 4; ++i)
        {
            const auto t {boost::int128::uint128_t {acc[static_cast<std::size_t>(i)]} * fw[j] + r[i + j] + carry};
            r[i + j] = t.low;
            carry = t.high;
        }
        r[4 + j] = carry;
    }
    acc = u256 {r[3], r[2], r[1], r[0]};
    return (r[4] | r[5]) == 0U;
}

// acc * 10^d, or false where it does not fit in 256 bits
constexpr auto beta_mul_pow10(u256& acc, int d) noexcept -> bool
{
    bool ok {true};
    for (; d > 0; d -= 38)
    {
        ok = beta_mul(acc, pow10(static_cast<boost::int128::uint128_t>(d < 38 ? d : 38))) && ok;
    }
    return ok;
}

// The other argument z = a * 10^-k of beta(n, z) for an integer n in [2, 35], or n = 0 for none
struct beta_integer
{
    boost::int128::uint128_t a;
    int k;
    unsigned n;
};

// Whether beta(n, z) = (n - 1)! / (z (z + 1) ... (z + n - 1)) is r * 10^q exactly: r * prod (a + i 10^k) * 10^q =
// (n - 1)! * 10^(kn), in 256 bits. Where it does not fit, the result is not exact. From k = 36 on, the parts of the
// factors prime to 10 multiply past (n - 1)!, and below, the factors fit in 128 bits.
constexpr auto beta_is_exact(const beta_integer& b, const boost::int128::uint128_t& r, const int q) noexcept -> bool
{
    if (b.k >= 36)
    {
        return false;
    }
    u256 lhs {r};
    u256 rhs {gamma_detail::gamma_at(b.n)};
    const auto step {pow10(static_cast<boost::int128::uint128_t>(b.k))};
    for (unsigned i {}; i < b.n; ++i)
    {
        if (!beta_mul(lhs, b.a + step * i))
        {
            return false;
        }
    }
    const int d {q - b.k * static_cast<int>(b.n)};
    return (d >= 0 ? beta_mul_pow10(lhs, d) : beta_mul_pow10(rhs, -d)) && lhs == rhs;
}

// Whether digits with a relative error below 2^-bits can round another way than the exact value: the boundaries
// are the multiples of half a unit in the last of p digits
template <typename digits_type>
constexpr auto beta_near_rounding(const digits_type& digits, const int p, const int bits) noexcept -> bool
{
    const auto m {pow10(static_cast<digits_type>(num_digits(digits) - p - 1)) * 5U};
    const auto r {digits % m};
    const auto w {(digits >> bits) + 1U};
    return r < w || r > m - w;
}

// beta(m, M) = 1 / m + d for a tiny m = sm * 10^em, with d > 0 for M < 1, d < 0 for M > 1 and |d| < 10^bound. An exact
// 1 / m moves to the side of d, else 1 / m + d rounds as 1 / m where |d| is below its gap to the boundaries, or NaN.
template <typename T>
constexpr auto beta_tiny(const boost::int128::uint128_t& sm, const int em, const bool above, const int bound) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto r {impl::div_mod(umul256(pow10(static_cast<boost::int128::uint128_t>(p)), pow10(static_cast<boost::int128::uint128_t>(p))), sm)};
    const auto q {static_cast<boost::int128::uint128_t>(r.quotient)};
    const auto rem {static_cast<boost::int128::uint128_t>(r.remainder)};
    const int e {-em - 2 * p};
    if (rem == 0U)
    {
        return bound < e ? T {above ? q * 10U + 1U : q * 10U - 1U, e - 1} : std::numeric_limits<T>::quiet_NaN();
    }
    // |d| < 10^bound <= 10^e * gap / sm, with gap = min(rem, sm - rem) >= 10^(digits of gap - 1)
    const auto gap {rem < sm - rem ? rem : sm - rem};
    return bound + num_digits(sm) <= e + num_digits(gap) - 1 ? T {q * 10U + 5U, e - 1} : std::numeric_limits<T>::quiet_NaN();
}

// beta(m, M) = gamma(m) M^-m e^-(F + d) for M = sM * 10^eM >= 10^15 and 0 < m = sm * 10^em < 10^4, with r = m / M:
// F = (M + m - 1/2) ln(1 + r) - m = r m A - ln(1 + r) / 2 for A = sum (-r)^(j - 1) / (j (j + 1)), d = -r (1 - r) / (12 M)
template <typename T, int N>
constexpr auto beta_large(const boost::int128::uint128_t& sm, const int em, const boost::int128::uint128_t& sM, const int eM) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto gm {gamma_detail::gamma_positive<N>(sm, em, p)};
    const auto mw {gamma_detail::gamma_wide<N>(sm, em)};
    const auto m14 {gamma_detail::gamma_shr(mw, 14U)};
    const int nd {num_digits(sM)};
    const auto ln_m {log_detail::log_wide_of<N>(sM, nd, eM, log_detail::log_base::e)};

    // m ln(M) = m * (ln(M) / 2^14) * 2^14
    const auto m_ln {gamma_detail::gamma_shl(fx_wide_mul(mw, gamma_detail::gamma_shr(ln_m, 14U)), 14U)};

    // 1 / M, r = m / M and a = r m
    const auto inv {gamma_detail::gamma_inv<N>(sM, nd, eM)};
    const auto r {fx_shl(fx_mul(inv, m14), 14U)};
    const auto a {fx_shl(fx_mul(r, m14), 14U)};

    // A and ln(1 + r) / 2 as alternating sums, which stay between 0 and their first terms
    fx<N> sum_a {};
    fx<N> half_ln {};
    auto pw {fx_one<N>()};
    for (std::uint64_t j {1}; fx_lead(pw) < 64 * N; ++j)
    {
        const auto ta {fx_div_small(pw, j * (j + 1U))};
        pw = fx_mul(pw, r);
        const auto tb {fx_div_small(pw, 2U * j)};
        sum_a = (j & 1U) != 0U ? fx_add(sum_a, ta) : fx_sub(sum_a, ta);
        half_ln = (j & 1U) != 0U ? fx_add(half_ln, tb) : fx_sub(half_ln, tb);
    }

    // st = st(m) - m ln(M) - F - d = st(m) - m ln(M) - a A - q r + ln(1 + r) / 2 + q, with q = r / (12 M)
    const auto q {fx_div_small(fx_mul(r, inv), 12U)};
    auto st {fx_wide_add(gm.st, m_ln, true)};
    st = fx_wide_add(st, fx_widen(fx_add(fx_mul(a, sum_a), fx_mul(q, r))), true);
    st = fx_wide_add(st, fx_widen(fx_add(half_ln, q)));
    const auto res {gamma_detail::gamma_exp_digits<T>(st, gamma_detail::gamma_frecip(gm.c))};
    if (res.kind != exp_detail::exp_kind::digits)
    {
        return exp_detail::exp10_beyond<T, N>(res.kind, false);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N < 3)
    {
        if (beta_near_rounding(res.digits, p, 72))
        {
            return beta_large<T, (N < 3 ? N + 1 : 3)>(sm, em, sM, eM);
        }
    }
    return fx_round<T>(res.digits, res.exp10, false);
}

// beta(x, y) = e^(st(x) + st(y) - st(x + y)) * c(x + y) / (c(x) c(y)) with gamma(z) = e^st(z) / c(z) in N words, and
// gamma(x + y + 1) / (x + y) for x + y < 1. One or two words compute it again in one more near a rounding boundary.
template <typename T, int N>
constexpr auto beta_words(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const beta_integer& b) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};

    const auto gx {gamma_detail::gamma_positive<N>(sx, ex, p)};
    const auto gy {gamma_detail::gamma_positive<N>(sy, ey, p)};
    const auto w {fx_wide_add(gamma_detail::gamma_wide<N>(sx, ex), gamma_detail::gamma_wide<N>(sy, ey))};
    gamma_parts<N> gs {};
    if (gamma_detail::gamma_whole(w) >= 1U)
    {
        gs = gamma_detail::gamma_shift(w);
    }
    else
    {
        gs = gamma_detail::gamma_shift(fx_wide_add(w, fx_wide_int<N>(1U)));
        gs.c = gamma_detail::gamma_fmul(gs.c, beta_small_sum<N>(sx, ex, sy, ey, p));
    }
    const auto c {gamma_detail::gamma_fmul(gs.c, gamma_detail::gamma_frecip(gamma_detail::gamma_fmul(gx.c, gy.c)))};
    const auto st {fx_wide_add(fx_wide_add(gx.st, gy.st), gs.st, true)};
    const auto r {gamma_detail::gamma_exp_digits<T>(st, c)};
    if (r.kind != exp_detail::exp_kind::digits)
    {
        return exp_detail::exp10_beyond<T, N>(r.kind, false);
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N < 3)
    {
        if (beta_near_rounding(r.digits, p, N == 1 ? 48 : 72))
        {
            return beta_words<T, (N < 3 ? N + 1 : 3)>(sx, ex, sy, ey, b);
        }
    }
    if (b.n != 0U)
    {
        // The value next to the digits, without its trailing zeros
        const int drop {num_digits(r.digits) - p};
        const auto unit {pow10(static_cast<boost::int128::uint128_t>(drop))};
        auto sig {(r.digits + unit / 2U) / unit};
        const auto t {remove_trailing_zeros(sig)};
        if (beta_is_exact(b, t.trimmed_number, r.exp10 + drop + static_cast<int>(t.number_of_removed_zeros)))
        {
            return T {sig, r.exp10 + drop};
        }
    }
    return fx_round<T>(r.digits, r.exp10, false);
}

} // namespace beta_detail

template <typename T>
constexpr auto beta_impl(const T x, const T y) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    // beta is defined for x > 0 and y > 0, and other arguments are a domain error
    if (isnan(x) || isnan(y) || x <= T {0} || y <= T {0})
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (isinf(x) || isinf(y))
    {
        return T {0};
    }
    #else
    if (x <= T {0} || y <= T {0})
    {
        return T {0};
    }
    #endif

    // beta(x, 1) = 1 / x
    if (x == T {1} || y == T {1})
    {
        return T {1} / (x == T {1} ? y : x);
    }
    constexpr int p {std::numeric_limits<T>::digits10};
    int ex {};
    int ey {};
    const auto sx {static_cast<boost::int128::uint128_t>(frexp10(x, &ex))};
    const auto sy {static_cast<boost::int128::uint128_t>(frexp10(y, &ey))};

    // A tiny m = min(x, y) next to M = max(x, y), with m < 10^mag_m and M < 10^mag_M
    const bool x_min {x < y};
    const auto s_m {x_min ? sx : sy};
    const int e_m {x_min ? ex : ey};
    const auto s_M {x_min ? sy : sx};
    const int e_M {x_min ? ey : ex};
    const int mag_m {e_m + p};
    const int mag_M {e_M + p};
    if (mag_m <= -20 && mag_m - mag_M <= -20)
    {
        const int bound {mag_M <= 0 ? 2 - mag_M : num_digits(static_cast<unsigned>(3 + 3 * mag_M)) + 1};
        const T r {beta_detail::beta_tiny<T>(s_m, e_m, x_min ? y < T {1} : x < T {1}, bound)};
        if (!isnan(r))
        {
            return r;
        }
    }

    constexpr int words {p <= 16 ? 2 : 3};
    if (ex + p <= 15 && ey + p <= 15)
    {
        // beta(n, z) can be exact for an integer n
        const bool ix {gamma_detail::gamma_is_integer(sx, ex)};
        const bool iy {gamma_detail::gamma_is_integer(sy, ey)};
        const bool use_x {ix && (!iy || x < y)};
        beta_detail::beta_integer b {0U, 0, 0U};
        if ((ix || iy) && (use_x ? x : y) <= T {35})
        {
            const auto sn {use_x ? sx : sy};
            const int en {use_x ? ex : ey};
            const auto tz {remove_trailing_zeros(use_x ? sy : sx)};
            const int ez {(use_x ? ey : ex) + static_cast<int>(tz.number_of_removed_zeros)};
            b.n = static_cast<unsigned>(gamma_detail::gamma_integer_of(sn, en));
            b.a = ez >= 0 ? tz.trimmed_number * pow10(static_cast<boost::int128::uint128_t>(ez)) : tz.trimmed_number;
            b.k = ez >= 0 ? 0 : -ez;
        }
        // One word for decimal32_t below 1000, else two words, or three for decimal128_t
        if (p <= 7 && ex + p <= 3 && ey + p <= 3)
        {
            return beta_detail::beta_words<T, (p <= 7 ? 1 : words)>(sx, ex, sy, ey, b);
        }
        return beta_detail::beta_words<T, words>(sx, ex, sy, ey, b);
    }

    // M >= 10^15, where m >= 10^4 gives beta(m, M) <= beta(10^4, 10^15) < e^-250000
    if (mag_m <= 4)
    {
        return beta_detail::beta_large<T, words>(s_m, e_m, s_M, e_M);
    }
    return T {1U, -2 * std::numeric_limits<T>::max_exponent10};
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto beta(const T y, const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::beta_impl(static_cast<evaluation_type>(y), static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_BETA_HPP
