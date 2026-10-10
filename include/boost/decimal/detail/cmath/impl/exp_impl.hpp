// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_EXP_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_EXP_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace exp_detail {

// The coefficients 1/(First + Step n)! for n < M, in N words
template <std::size_t First = 0U, std::size_t Step = 1U, int N, std::size_t M>
constexpr auto exp_head(fx<N> (&c)[M]) noexcept -> void
{
    static_assert(First + Step * (M - 1U) < sizeof(fx_table::inv_factorial) / sizeof(fx_table::inv_factorial[0]), "1/n! is in the table");

    for (std::size_t n {}; n < M; ++n)
    {
        c[n] = fx_top<N>(fx_table::inv_factorial[First + Step * n]);
    }
}

// exp(-z) for 0 <= z <= ln(2) / 256 by the Taylor series, whose terms alternate in sign. The terms from
// z^4 on need fewer words. The error is below 2^-60, 2^-100 and 2^-160 for one, two and three words.
constexpr auto exp_neg_series(const fx<1>& z) noexcept -> fx<1>
{
    fx<1> head[7] {};
    exp_head(head);
    return fx_alternating(z, head);
}

constexpr auto exp_neg_series(const fx<2>& z) noexcept -> fx<2>
{
    fx<2> head[4] {};
    exp_head(head);
    return fx_alternating_split<fx_table::exp_d64_shift>(z, head, fx_table::exp_d64_tail);
}

constexpr auto exp_neg_series(const fx<3>& z) noexcept -> fx<3>
{
    fx<3> head[4] {};
    exp_head(head);
    return fx_alternating_split<fx_table::exp_d128_mid_shift, fx_table::exp_d128_tail_shift>(z, head, fx_table::exp_d128_mid, fx_table::exp_d128_tail);
}

// exp(-z) as above in two words with the error below 2^-124: terms to z^6 in two words, and z^7 to z^11 in one
constexpr auto exp_neg_series_wide(const fx<2>& z) noexcept -> fx<2>
{
    fx<2> head[7] {};
    exp_head(head);
    return fx_alternating_split<fx_table::exp_d128_fast_shift>(z, head, fx_table::exp_d128_fast_tail);
}

// exp(-z) as above to z^3 only, with the error below z^4 / 24 < 2^-38
template <int N>
constexpr auto exp_neg_series_short(const fx<N>& z) noexcept -> fx<N>
{
    fx<N> head[4] {};
    exp_head(head);
    return fx_alternating(z, head);
}

template <int N>
constexpr auto exp_neg_series_of(const fx<N>& z, std::false_type) noexcept -> fx<N>
{
    return exp_neg_series(z);
}

constexpr auto exp_neg_series_of(const fx<2>& z, std::true_type) noexcept -> fx<2>
{
    return exp_neg_series_wide(z);
}

// exp(-x) is below half of the smallest subnormal value and exp(x), sinh(x) and cosh(x) overflow for |x|
// at or above this bound, and the same for exp2(x) with Base2. As sig has p digits, |x| < bound gives e < 0.
template <typename T, bool Base2 = false>
constexpr auto exp_bound() noexcept -> std::uint64_t
{
    // log2(10) and ln(10) times 10^6, rounded up
    const auto c {Base2 ? UINT64_C(3321929) : UINT64_C(2302586)};
    return static_cast<std::uint64_t>(std::numeric_limits<T>::max_exponent10 + std::numeric_limits<T>::digits10 + 1) * c / UINT64_C(1000000) + 1U;
}

// The integer part of |x| = sig * 10^-k, or a value above each bound for k <= 0, in 64 bits where sig fits
constexpr auto exp_whole(const boost::int128::uint128_t& sig, const int k) noexcept -> boost::int128::uint128_t
{
    if (k > 0 && k < 20 && sig.high == 0U)
    {
        return sig.low / pow10(static_cast<std::uint64_t>(k));
    }
    return k <= 0 ? boost::int128::uint128_t {UINT64_MAX} : k < 39 ? sig / pow10(static_cast<boost::int128::uint128_t>(k)) : boost::int128::uint128_t {};
}

// The fraction of |x| = sig * 10^-k, for the integer part whole
template <int words>
constexpr auto exp_frac(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole) noexcept -> fx<words>
{
    return fx_from_decimal<words>(k < 39 ? sig - static_cast<boost::int128::uint128_t>(whole) * pow10(static_cast<boost::int128::uint128_t>(k)) : sig, k);
}

// Whether |x| = whole + frac is below 1/16, where the Taylor series of sinh, cosh, tanh and expm1 converge fast
template <int words>
constexpr auto exp_is_small(const std::uint64_t whole, const fx<words>& frac) noexcept -> bool
{
    return whole == 0U && frac.w[words - 1] < (UINT64_C(1) << 58U);
}

// The number of terms in x^2 of these series for |x| < 1/16 in words words
constexpr auto exp_small_terms(const int words) noexcept -> std::size_t
{
    return words == 1 ? 5U : words == 2 ? 8U : 12U;
}

// expm1(x) / x = 1 + x R for x = frac, or 1 - x R for x = -frac with neg, where R = even(x^2) + x odd(x^2) and both
// series start from 1/2! and 1/3!. Returns R, and the ratio in s.
template <int words, std::size_t terms>
constexpr auto expm1_rest(const fx<words>& frac, const bool neg, fx<words>& s) noexcept -> fx<words>
{
    fx<words> even[terms - 1U] {};
    fx<words> odd[terms - 1U] {};
    exp_head<2U, 2U>(even);
    exp_head<3U, 2U>(odd);
    const auto u {fx_mul(frac, frac)};
    const auto h {fx_horner(u, even)};
    const auto o {fx_mul(frac, fx_horner(u, odd))};
    const auto r {neg ? fx_sub(h, o) : fx_add(h, o)};
    s = neg ? fx_sub(fx_one<words>(), fx_mul(frac, r)) : fx_add(fx_one<words>(), fx_mul(frac, r));
    return r;
}

// y = exp10 + f with 0 <= f < 1
template <int words>
struct exp_parts
{
    fx<words> f;
    int exp10;
};

// y = |x| * c for |x| = whole + frac, with c = log10(e) or log10(2). In one word, the error of f is below 3.5 units
// of 2^-62.
template <int words>
constexpr auto exp_log10(const std::uint64_t whole, const fx<words>& frac, const fx<4>& c) noexcept -> exp_parts<words>
{
    // whole * c + frac * c in words + 2 words, scaled by 2^(64 words + 62)
    std::uint64_t y[static_cast<std::size_t>(words + 2)] {};
    y[words + 1] = fx_mul_add(fx_make_span<words + 1>(y), fx_make_span<words + 1>(c.w, 3 - words), whole);
    const auto low {fx_mul(frac, fx_top<words>(c))};
    y[words + 1] += fx_add_to(fx_make_span<words>(y, 1), fx_make_span(low.w));

    exp_parts<words> r {{}, static_cast<int>((y[words + 1] << 2U) | (y[words] >> 62U))};
    for (int i {}; i < words; ++i)
    {
        r.f.w[i] = y[i + 1];
    }
    r.f.w[words - 1] &= (UINT64_C(1) << 62U) - 1U;
    return r;
}

// 10^f = 2^a * m with 1 <= m < 2
template <int words>
struct exp_mantissa
{
    fx<words> m;
    unsigned a;
};

// 10^f = 2^(a + i/16 + j/256) * exp(-v ln(2)) for 0 <= f <= 1 and 0 < v <= 1/256, with Short the series to z^3. In
// one word otherwise, m has a relative error below 6.4 units of 2^-62 plus ln(10) times the error of f.
template <int words, bool Wide = false, bool Short = false>
constexpr auto exp10_frac(const fx<words>& f) noexcept -> exp_mantissa<words>
{
    // t = f * log2(10) = a + i/16 + (j - 1)/256 + r, and v = 1/256 - r
    auto t {fx_mul(f, fx_top<words>(fx_table::log2_10))};
    const auto a {static_cast<unsigned>(t.w[words - 1] >> 62U)};
    const auto i {static_cast<std::size_t>((t.w[words - 1] >> 58U) & 15U)};
    const auto j {static_cast<std::size_t>((t.w[words - 1] >> 54U) & 15U) + 1U};
    t.w[words - 1] &= (UINT64_C(1) << 54U) - 1U;
    fx<words> step {};
    step.w[words - 1] = UINT64_C(1) << 54U;
    const auto v {fx_sub(step, t)};

    const auto coarse {fx_mul(fx_top<words>(fx_table::pow2[i]), fx_top<words>(fx_table::pow2_fine[j]))};
    const auto z {fx_mul(v, fx_top<words>(fx_table::ln2))};
    return {fx_mul(coarse, Short ? exp_neg_series_short(z) : exp_neg_series_of(z, std::integral_constant<bool, Wide> {})), a};
}

// For decimal64_t, m * s * 10^e10 for an exact m and s with an error below err units of 2^-62, rounded once in out
// where it cannot be next to a rounding boundary: 19 digits from m * s with m * s in [10^17, 10^20)
template <typename T>
constexpr auto exp_try_round(boost::int128::uint128_t m, const fx<1>& s, int e10, const bool neg, const std::uint64_t err, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;

    auto d {fx_scale(m, s)};
    if (d < pow10(static_cast<sig_type>(18)))
    {
        m *= 10U;
        --e10;
        d = fx_scale(m, s);
    }
    auto w {((m * err) >> 62U).low + 2U};
    if (d >= pow10(static_cast<sig_type>(19)))
    {
        d /= 10U;
        ++e10;
        w = w / 10U + 2U;
    }
    return fx_try_round(d.low, e10, neg, w, out);
}

// -y = (-exp10 - 1) + (1 - f) for y = exp10 + f with 0 <= f < 1
template <int words>
BOOST_DECIMAL_FORCE_INLINE constexpr auto exp_negate(int& exp10, fx<words>& f) noexcept -> void
{
    exp10 = -exp10 - 1;
    f = fx_sub(fx_one<words>(), f);
}

// y = |x| * c for |x| = sig * 10^-k with the integer part whole, or -y = (-exp10 - 1) + (1 - f) for neg
template <int words>
constexpr auto exp_log10_signed(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole, const fx<4>& c, const bool neg) noexcept -> exp_parts<words>
{
    auto y {exp_log10(whole, exp_frac<words>(sig, k, whole), c)};
    if (neg)
    {
        exp_negate(y.exp10, y.f);
    }
    return y;
}

// 10^f = 2^a * m as scale10 + 1 digits: m / 2 < 1 for fx_scale, and 10^scale10 * 2^(a + 1) * m / 2 is 10^f
template <int words>
constexpr auto exp10_digits(const exp_mantissa<words>& t, const int scale10) noexcept -> boost::int128::uint128_t
{
    return fx_scale(pow10(static_cast<boost::int128::uint128_t>(scale10)) << (t.a + 1U), fx_shr(t.m, 1U));
}

// A value in (1 - 10^-(scale10 + 1), 1), rounded once in the current mode with the sign neg
template <typename T, typename digits_type>
constexpr auto below_one(const int scale10, const bool neg) noexcept -> T
{
    return fx_round<T>(static_cast<digits_type>(pow10(static_cast<boost::int128::uint128_t>(scale10 + 1)) - 1U), -scale10 - 1, neg);
}

enum class exp_kind { digits, overflow, underflow, above_one, below_one };

// A result as digits * 10^exp10 with scale10 + 1 digits, or the kind of a result beyond them
template <int words>
struct exp_digits
{
    fx_digits_t<words> digits;
    int exp10;
    exp_kind kind;
};

// exp(x), or exp2(x) with Base2, for a finite x != 0 where the result is not exact, rounded once in
// the current mode. With y = x * log10(e) or x * log10(2) = k + f and 0 <= f < 1, the result is 10^k * 10^f.
template <typename T, bool Base2 = false>
constexpr auto exp_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    // One word gives 19 digits of the result, and two or three words give 38
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    const int k {-e};
    const auto wide_whole {exp_whole(sig, k)};
    const auto whole {static_cast<std::uint64_t>(wide_whole)};
    if (wide_whole >= exp_bound<T, Base2>())
    {
        return neg ? T {1, -2 * std::numeric_limits<T>::max_exponent10} : T {1, 2 * std::numeric_limits<T>::max_exponent10};
    }

    // For |x| < 10^-(p + 1), exp(x) rounds as a value just above or just below 1
    if (k >= 2 * p + 1)
    {
        return neg ? below_one<T, digits_type>(scale10, false)
                   : fx_round<T>(static_cast<digits_type>(pow10(static_cast<boost::int128::uint128_t>(scale10))), -scale10, false);
    }
    const auto& c {Base2 ? fx_table::log10_2 : fx_table::log10_e};

    BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
    {
        // A first try with Short: its error is below 2^-38 of the result, so below d1 * 2^-37 + 2 units of d1
        const auto y1 {exp_log10_signed<1>(sig, k, whole, c, neg)};
        const auto d1 {static_cast<std::uint64_t>(exp10_digits(exp10_frac<1, false, true>(y1.f), fx_scale10(1)))};
        T r {};
        if (fx_try_round(d1, y1.exp10 - fx_scale10(1), false, (d1 >> 37U) + 2U, r))
        {
            return r;
        }
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        // A first try in one word: m has a relative error below 14.5 units of 2^-62 and m / 2 below 2 more, so the
        // error is below d1 * 2^-58 + 3 units of d1
        const auto y1 {exp_log10_signed<1>(sig, k, whole, c, neg)};
        const auto d1 {static_cast<std::uint64_t>(exp10_digits(exp10_frac(y1.f), fx_scale10(1)))};
        T r {};
        if (fx_try_round(d1, y1.exp10 - fx_scale10(1), false, (d1 >> 58U) + 3U, r))
        {
            return r;
        }
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // A first try in two words: its error is below 46 units of 2^-126, so below d2 * 2^-119 + 2 units of d2
        const auto y2 {exp_log10_signed<2>(sig, k, whole, c, neg)};
        const auto d2 {exp10_digits(exp10_frac<2, true>(y2.f), fx_scale10(2))};
        if (!fx_near_boundary(d2, static_cast<std::uint64_t>(d2 >> 119U) + 2U))
        {
            return fx_round<T>(static_cast<digits_type>(d2), y2.exp10 - fx_scale10(2), false);
        }
    }

    const auto y {exp_log10_signed<words>(sig, k, whole, c, neg)};
    const auto d {exp10_digits(exp10_frac(y.f), scale10)};
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // Next to a rounding boundary for a small x, exp(x) = 1 + x S with S = expm1(x) / x, and 2^x = 1 + x ln(2) S at
        // x ln(2): the rest past 1 takes its own precision, with 1 as 10^37, or 10^38 below 1
        if (whole == 0U && fx_near_boundary(d, 2U))
        {
            auto frac {exp_frac<words>(sig, k, whole)};
            if (Base2)
            {
                frac = fx_mul(frac, fx_top<words>(fx_table::ln2));
            }
            if (exp_is_small(whole, frac))
            {
                fx<words> s {};
                static_cast<void>(expm1_rest<words, exp_small_terms(words)>(frac, neg, s));
                int a {};
                auto g {fx_mul(fx_lead_t<words>(sig, k, a), s)};
                if (Base2)
                {
                    g = fx_mul(g, fx_top<words>(fx_table::ln2));
                }
                const auto one {pow10(static_cast<boost::int128::uint128_t>(neg ? 38 : 37))};
                return fx_round<T>(static_cast<digits_type>(fx_lead_digits(one, g, a, neg)), neg ? -38 : -37, false);
            }
        }
    }
    return fx_round<T>(static_cast<digits_type>(d), y.exp10 - scale10, false);
}

} // namespace exp_detail
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_EXP_IMPL_HPP
