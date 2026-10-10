// Copyright John Maddock 2006.
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP

#include <boost/decimal/detail/cmath/impl/exp_impl.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace erf_detail {

// The polynomial of g on interval i for the given words: its first coefficient in the table, the number of
// coefficients, its center c = cw + 1/2 for half, log2 of its half width h, and the first coefficient of the tail
struct erf_interval
{
    int first;
    int count;
    unsigned cw;
    int log2_h;
    int split;
    bool half;
};

constexpr auto erf_interval_of(const int words, const unsigned i) noexcept -> erf_interval
{
    // [1, 2], [2, 4] and [4, 8] for all types, then [8, 10] for two words, and [8, 12] and [12, 16] for three words
    return words == 1 ? (i == 0U   ? erf_interval {0, 20, 1U, -1, 20, true}
                         : i == 1U ? erf_interval {20, 23, 3U, 0, 23, false}
                                   : erf_interval {43, 26, 6U, 1, 26, false})
         : words == 2 ? (i == 0U   ? erf_interval {0, 31, 1U, -1, 23, true}
                         : i == 1U ? erf_interval {31, 36, 3U, 0, 28, false}
                         : i == 2U ? erf_interval {67, 40, 6U, 1, 33, false}
                                   : erf_interval {107, 27, 9U, 0, 19, false})
                      : (i == 0U   ? erf_interval {0, 43, 1U, -1, 23, true}
                         : i == 1U ? erf_interval {43, 50, 3U, 0, 28, false}
                         : i == 2U ? erf_interval {93, 57, 6U, 1, 33, false}
                         : i == 3U ? erf_interval {150, 49, 10U, 1, 25, false}
                                   : erf_interval {199, 44, 14U, 1, 21, false});
}

// Whether each first is the sum of the counts before it, and the last interval ends at the size of the table
constexpr auto erf_intervals_fit(const int words, const unsigned last, const int size) noexcept -> bool
{
    for (unsigned i {1}; i <= last; ++i)
    {
        if (erf_interval_of(words, i).first != erf_interval_of(words, i - 1U).first + erf_interval_of(words, i - 1U).count)
        {
            return false;
        }
    }
    return erf_interval_of(words, last).first + erf_interval_of(words, last).count == size;
}

static_assert(erf_intervals_fit(1, 2U, static_cast<int>(sizeof(fx_table::erfc_d32) / sizeof(fx_table::erfc_d32[0]))), "the intervals cover the table");
static_assert(erf_intervals_fit(2, 3U, static_cast<int>(sizeof(fx_table::erfc_d64) / sizeof(fx_table::erfc_d64[0]))), "the intervals cover the table");
static_assert(erf_intervals_fit(3, 4U, static_cast<int>(sizeof(fx_table::erfc_d128) / sizeof(fx_table::erfc_d128[0]))), "the intervals cover the table");

// Whether the tail of each interval, the coefficients from its split on, sums to below 2^-62: for |u| <= 1 its
// partial sums then fit in N - 1 words scaled by 2^64
template <int N, std::size_t M>
constexpr auto erf_splits_fit(const fx<N> (&c)[M], const unsigned last) noexcept -> bool
{
    for (unsigned i {}; i <= last; ++i)
    {
        const auto iv {erf_interval_of(N, i)};
        boost::int128::uint128_t sum {};
        for (int k {iv.first + iv.split}; k < iv.first + iv.count; ++k)
        {
            if (c[k].w[N - 1] != 0U)
            {
                return false;
            }
            sum += c[k].w[N - 2] + UINT64_C(1);
        }
        if (sum.high != 0U)
        {
            return false;
        }
    }
    return true;
}

static_assert(erf_splits_fit(fx_table::erfc_d64, 3U) && erf_splits_fit(fx_table::erfc_d128, 4U), "the tails fit in fewer words");

// The coefficients of the polynomial of the interval iv in the table of N words
constexpr auto erf_poly(const fx<1>&, const erf_interval& iv) noexcept -> fx_span<const fx<1>> { return fx_make_span(fx_table::erfc_d32, iv.first, iv.count); }
constexpr auto erf_poly(const fx<2>&, const erf_interval& iv) noexcept -> fx_span<const fx<2>> { return fx_make_span(fx_table::erfc_d64, iv.first, iv.count); }
constexpr auto erf_poly(const fx<3>&, const erf_interval& iv) noexcept -> fx_span<const fx<3>> { return fx_make_span(fx_table::erfc_d128, iv.first, iv.count); }

// erf(x) / x for z = x^2 <= 1 by the Taylor series, whose terms alternate in sign. The terms from z^20, z^33 and z^45
// on are below 2^-64, 2^-128 and 2^-192. As z can be near 1, all terms take N words.
constexpr auto erf_series(const fx<1>& z) noexcept -> fx<1>
{
    return fx_even_odd(z, fx_make_span<20>(fx_table::erf_taylor), false);
}

template <int N>
constexpr auto erf_series(const fx<N>& z) noexcept -> fx<N>
{
    constexpr std::size_t terms {N == 2 ? 33U : 45U};
    return fx_series<true, terms, terms, terms>(z, fx_table::erf_taylor);
}

// x = whole + frac as x / 4 in fixed point, for x < 16
template <int N>
constexpr auto erf_quarter(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    auto r {fx_shr(frac, 2U)};
    r.w[N - 1] += whole << 60U;
    return r;
}

// The polynomial c in u >= 0, alternating where above, else with positive terms
constexpr auto erf_g_sum(const fx<1>& u, const fx_span<const fx<1>> c, const bool above, std::size_t) noexcept -> fx<1>
{
    return fx_even_odd(u, c, !above);
}

// The tail from c[split] on runs in N - 1 words scaled by 2^64, the low words of c, with the same last bit
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto erf_g_sum(const fx<N>& u, const fx_span<const fx<N>> c, const bool above, const std::size_t split) noexcept -> fx<N>
{
    const auto ut {fx_top<N - 1>(u)};
    auto t {fx_low<N - 1>(c[c.size() - 1U].w)};
    for (std::size_t k {c.size() - 1U}; k-- > split;)
    {
        const auto ck {fx_low<N - 1>(c[k].w)};
        t = above ? fx_sub(ck, fx_mul(ut, t)) : fx_add(ck, fx_mul(ut, t));
    }
    fx<N> s {};
    for (int i {}; i < N - 1; ++i)
    {
        s.w[i] = t.w[i];
    }
    for (std::size_t k {split}; k-- > 0U;)
    {
        s = above ? fx_sub(c[k], fx_mul(u, s)) : fx_add(c[k], fx_mul(u, s));
    }
    return s;
}

// g(x) = erfc(x) exp(x^2) for 1 <= x = whole + frac below the asymptotic range: the polynomial of the interval in
// u = (x - c) / h, which is alternating for u >= 0 and has positive terms in |u| for u < 0
template <int N>
constexpr auto erf_g_arg(const std::uint64_t whole, const fx<N>& frac, const erf_interval& iv, bool& above) noexcept -> fx<N>
{
    const auto x4 {erf_quarter(whole, frac)};
    fx<N> c4 {};
    c4.w[N - 1] = (static_cast<std::uint64_t>(iv.cw) << 60U) | (iv.half ? UINT64_C(1) << 59U : 0U);
    above = !fx_less(x4, c4);
    return fx_shl(above ? fx_sub(x4, c4) : fx_sub(c4, x4), static_cast<unsigned>(2 - iv.log2_h));
}

// g(x) from the polynomial of its interval in N words
template <int N>
constexpr auto erf_g_table(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    const unsigned i {whole < 2U ? 0U : whole < 4U ? 1U : whole < 8U ? 2U : N == 3 && whole >= 12U ? 4U : 3U};
    const auto iv {erf_interval_of(N, i)};
    bool above {};
    const auto u {erf_g_arg(whole, frac, iv, above)};
    return erf_g_sum(u, erf_poly(u, iv), above, static_cast<std::size_t>(iv.split));
}

// g(x) = m * 2^-(64N - 66 + l) in N - 1 words, for the table of N words of decimal64_t or decimal128_t
template <int N>
using erf_mant_t = std::conditional_t<N == 2, std::uint64_t, boost::int128::uint128_t>;

template <int N>
struct erf_float
{
    erf_mant_t<N> m;
    int l;
};

// The top N - 1 words of c * 2^l
template <int N>
constexpr auto erf_float_m(const fx<N>& c, const int l) noexcept -> erf_mant_t<N>
{
    std::uint64_t w[2] {};
    for (int i {N - 2}; i < N; ++i)
    {
        const int j {i - l / 64};
        const auto r {static_cast<unsigned>(l % 64)};
        const auto hi {j >= 0 ? c.w[j] : UINT64_C(0)};
        const auto lo {j >= 1 ? c.w[j - 1] : UINT64_C(0)};
        w[i - (N - 2)] = r == 0U ? hi : (hi << r) | (lo >> (64U - r));
    }
    return static_cast<erf_mant_t<N>>(boost::int128::uint128_t {w[1], w[0]} >> (N == 2 ? 64 : 0));
}

// The words of g as a mantissa
constexpr auto erf_mant(const fx<1>& g) noexcept -> std::uint64_t { return g.w[0]; }
constexpr auto erf_mant(const fx<2>& g) noexcept -> boost::int128::uint128_t { return {g.w[1], g.w[0]}; }

// l for c in [2^(64N - 66 - l), 2^(64N - 65 - l))
template <int N>
constexpr auto erf_float_l(const fx<N>& c) noexcept -> int
{
    int z {};
    for (int i {N - 1}; i >= 0 && c.w[i] == 0U; --i)
    {
        z += 64;
    }
    return z + (z < 64 * N ? countl_zero(c.w[N - 1 - z / 64]) : 0) - 1;
}

// For the first M coefficients of the table of N words, m in [2^(64N - 67), 2^(64N - 65)) with each l below the next
// one, the shift l_(k+1) - l_k - 1 that takes 2 u s_(k+1) times 2^-(64N - 64) to the scale of c_k, l of each first
// coefficient, and for |u| < (j + 1) / 8 the last term and for three words the first term below 2^-64 of the first
template <int N, std::size_t M>
struct erf_float_table
{
    erf_mant_t<N> m[M];
    int l[4];
    unsigned char last[4][8];
    unsigned char split[4][8];
    unsigned char shift[M];
};

template <int N>
using erf_float_table_t = erf_float_table<N, N == 2 ? 134U : 199U>;

template <int N, std::size_t M, std::size_t S>
constexpr auto erf_float_table_of(const fx<N> (&c)[S]) noexcept -> erf_float_table<N, M>
{
    erf_float_table<N, M> t {};
    for (unsigned i {}; i <= 3U; ++i)
    {
        const auto iv {erf_interval_of(N, i)};
        int next {};
        for (int k {iv.first + iv.count - 1}; k >= iv.first; --k)
        {
            const int own {erf_float_l(c[k])};
            const int l {k + 1 == iv.first + iv.count || own < next ? own : next - 1};
            t.m[k] = erf_float_m(c[k], l);
            t.shift[k] = static_cast<unsigned char>(k + 1 == iv.first + iv.count ? 0 : next - l - 1);
            next = l;
        }
        t.l[i] = next;
        for (unsigned j {}; j < 8U; ++j)
        {
            // The terms after the last sum to below 2^-130 of the scale 2^-l of the first
            fx<N> um {};
            um.w[N - 1] = static_cast<std::uint64_t>(j + 1U) << 59U;
            fx<N> pw[64] {fx_one<N>()};
            for (int k {1}; k < iv.count; ++k)
            {
                pw[k] = fx_mul(pw[k - 1], um);
            }
            fx<N> tail {};
            fx<N> bound {};
            bound.w[0] = UINT64_C(1) << (60 - next);
            int last {iv.first + iv.count - 1};
            while (last > iv.first + 1 && !fx_less(bound, tail = fx_add(tail, fx_mul(c[last], pw[last - iv.first]))))
            {
                --last;
            }
            int split {iv.first};
            while (N == 3 && split < last && erf_float_l(fx_mul(c[split], pw[split - iv.first])) < next + 65)
            {
                ++split;
            }
            t.last[i][j] = static_cast<unsigned char>(last);
            t.split[i][j] = static_cast<unsigned char>(split);
        }
    }
    return t;
}

// The tables of decimal64_t on [1, 10] and of decimal128_t on [1, 12]
constexpr auto erf_float_table_for(std::integral_constant<int, 2>) noexcept -> erf_float_table_t<2> { return erf_float_table_of<2, 134>(fx_table::erfc_d64); }
constexpr auto erf_float_table_for(std::integral_constant<int, 3>) noexcept -> erf_float_table_t<3> { return erf_float_table_of<3, 199>(fx_table::erfc_d128); }

template <int N>
struct erf_float_imp
{
    static constexpr erf_float_table_t<N> t {erf_float_table_for(std::integral_constant<int, N> {})};
};

#if !(defined(__cpp_inline_variables) && __cpp_inline_variables >= 201606L) && (!defined(_MSC_VER) || _MSC_VER != 1900)

template <int N>
constexpr erf_float_table_t<N> erf_float_imp<N>::t;

#endif

// The high half of a b, in two words without the low product and so below the floor by at most 1
constexpr auto erf_mul_high(const std::uint64_t a, const std::uint64_t b) noexcept -> std::uint64_t
{
    std::uint64_t hi {};
    boost::int128::detail::umul(a, b, hi);
    return hi;
}

constexpr auto erf_mul_high(const boost::int128::uint128_t& a, const boost::int128::uint128_t& b) noexcept -> boost::int128::uint128_t
{
    const boost::int128::uint128_t lh {boost::int128::uint128_t {a.low} * b.high};
    const boost::int128::uint128_t hl {boost::int128::uint128_t {a.high} * b.low};
    const boost::int128::uint128_t mid {boost::int128::uint128_t {lh.low} + boost::int128::uint128_t {hl.low}};
    return boost::int128::uint128_t {a.high} * b.high + boost::int128::uint128_t {lh.high} + boost::int128::uint128_t {hl.high} + boost::int128::uint128_t {mid.high};
}

// The top word of m, and the mantissa in N - 1 words with the top word m
constexpr auto erf_top(const std::uint64_t m) noexcept -> std::uint64_t { return m; }
constexpr auto erf_top(const boost::int128::uint128_t& m) noexcept -> std::uint64_t { return m.high; }

template <int N>
constexpr auto erf_widen(const std::uint64_t m) noexcept -> erf_mant_t<N>
{
    return static_cast<erf_mant_t<N>>(boost::int128::uint128_t {m, 0U} >> (N == 2 ? 64 : 0));
}

// The mantissa m in the words of s
constexpr auto erf_as(const std::uint64_t m, std::uint64_t) noexcept -> std::uint64_t { return m; }
constexpr auto erf_as(const boost::int128::uint128_t& m, std::uint64_t) noexcept -> std::uint64_t { return m.high; }
constexpr auto erf_as(const boost::int128::uint128_t& m, const boost::int128::uint128_t&) noexcept -> boost::int128::uint128_t { return m; }

// From s = s_hi to s_lo, with u and w = u^2 times 2^(64W - 1) in the W words of s: s_k = c_k -+ u c_(k+1) + w s_(k+2),
// and first s_k = c_k -+ u s_(k+1) where hi - lo is odd, so that the two products of a step do not wait on each other
template <typename Table, typename M>
constexpr auto erf_run(const Table& t, int k, const int lo, M s, const M& u, const M& w, const bool above) noexcept -> M
{
    if (((k - lo) & 1) != 0)
    {
        --k;
        const auto p {static_cast<M>(erf_mul_high(u, s) >> t.shift[k])};
        s = above ? erf_as(t.m[k], s) - p : erf_as(t.m[k], s) + p;
    }
    while (k > lo)
    {
        k -= 2;
        const auto p1 {static_cast<M>(erf_mul_high(u, erf_as(t.m[k + 1], s)) >> t.shift[k])};
        const auto p2 {static_cast<M>(erf_mul_high(w, s) >> (t.shift[k] + t.shift[k + 1] + 1U))};
        s = above ? erf_as(t.m[k], s) - p1 + p2 : erf_as(t.m[k], s) + p1 + p2;
    }
    return s;
}

// g(x) for 1 <= x < 10 or 12 in N - 1 words by Horner's rule on the table of N words, each sum at the scale of its
// coefficient: g is below 4.25 units of 2^-62 or 8.34 units of 2^-126 off, from the steps, the tail and u
template <int N>
constexpr auto erf_g_float(const std::uint64_t whole, const fx<2>& frac) noexcept -> erf_float<N>
{
    const unsigned i {whole < 2U ? 0U : whole < 4U ? 1U : whole < 8U ? 2U : 3U};
    const auto iv {erf_interval_of(N, i)};
    bool above {};
    const auto u2 {static_cast<erf_mant_t<N>>(erf_mant(fx_top<N - 1>(erf_g_arg(whole, frac, iv, above))) << 1U)};
    const auto& t {erf_float_imp<N>::t};
    const auto u1 {erf_top(u2)};
    const auto w2 {static_cast<erf_mant_t<N>>(erf_mul_high(u2, u2) << 1U)};
    const auto w1 {erf_top(w2)};
    // |u| <= 1, and |u| = 1 at the left end of an interval takes the row of |u| < 1
    const auto j {static_cast<std::size_t>((u1 - (u1 >> 63U)) >> 60U)};
    const int last {t.last[i][j]};
    const auto s1 {erf_run(t, last, t.split[i][j], erf_top(t.m[last]), u1, w1, above)};
    const auto s {erf_run(t, t.split[i][j], iv.first, erf_widen<N>(s1), u2, w2, above)};
    return {s, t.l[i]};
}

// g(x) for x >= 10 or 12 in N - 1 words: 1 / x = r 2^-s with r = 1 / d for x = d 2^s, the series S in v = 1 / (2x^2) by
// Horner's rule from term 9 + 80 / whole or 9 + 384 / whole, one word above 6 + 120 / whole, and g = S r (2 / sqrt(pi))
// 2^-(s + 1), below 10 units of 2^-62 or 43.2 units of 2^-126 off
template <int N>
constexpr auto erf_g_asymptotic_float(const std::uint64_t whole, const fx<2>& frac) noexcept -> erf_float<N>
{
    const auto s {static_cast<unsigned>(63 - countl_zero(whole))};
    auto d {fx_shr(frac, s)};
    d.w[1] += whole << (62U - s);
    const auto dt {fx_top<N - 1>(d)};
    auto r {fx_recip(dt)};
    BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
    {
        r = fx_mul(r, fx_sub(fx_int<N - 1>(2U), fx_mul(dt, r)));
    }
    const auto v {fx_shr(fx_mul(r, r), 2U * s + 1U)};
    const auto v1 {fx_top<1>(v)};
    auto t1 {fx_one<1>()};
    std::uint64_t k {9U + (N == 2 ? 80U : 384U) / whole};
    const std::uint64_t k2 {N == 2 ? 0U : 6U + 120U / whole};
    while (k-- > k2)
    {
        t1 = fx_sub(fx_one<1>(), fx_mul(fx_mul_small(v1, 2U * k + 1U), t1));
    }
    auto t {fx_extend<N - 1>(t1)};
    for (++k; k-- > 0U;)
    {
        t = fx_sub(fx_one<N - 1>(), fx_mul(fx_mul_small(v, 2U * k + 1U), t));
    }
    const auto g {fx_mul(fx_mul(t, r), fx_shl(fx_top<N - 1>(fx_table::inv_sqrt_pi), 1U))};
    return {erf_mant(g), static_cast<int>(s) + 1};
}

// g(x) for x >= 8, 10 or 16 in one, two or three words by the asymptotic series 1 / (x sqrt(pi)) * sum (-1)^n
// (2n - 1)!! / (2x^2)^n, where the first term left out is below 2^-(64N - 2)
template <int N>
constexpr auto erf_g_asymptotic(const std::uint64_t whole, const fx<N>& frac) noexcept -> fx<N>
{
    // 1 / x = 1 / d / 2^s with x = d * 2^s and d in [1, 2)
    const auto s {static_cast<unsigned>(63 - countl_zero(whole))};
    auto d {fx_shr(frac, s)};
    d.w[N - 1] += whole << (62U - s);
    const auto inv_x {fx_shr(fx_recip(d), s)};
    const auto v {fx_shr(fx_mul(inv_x, inv_x), 1U)};

    // The sum stops where a term is below 4 units of the last bit; 200 terms bound the loop
    auto term {fx_one<N>()};
    auto pos {term};
    fx<N> neg {};
    for (std::uint64_t n {1U}; n < 200U; ++n)
    {
        const auto t {fx_mul(term, v)};
        term = fx_mul_small(t, 2U * n - 1U);
        bool small {term.w[0] < 4U};
        for (int i {1}; i < N; ++i)
        {
            small = small && term.w[i] == 0U;
        }
        if (small)
        {
            break;
        }
        if ((n & 1U) != 0U)
        {
            neg = fx_add(neg, term);
        }
        else
        {
            pos = fx_add(pos, term);
        }
    }
    return fx_mul(fx_mul(fx_sub(pos, neg), inv_x), fx_top<N>(fx_table::inv_sqrt_pi));
}

// erf(x), or erfc(x) with complement, as digits * 10^exp10 for 0 < |x| < 1: erf(x) = x * s(x^2) by the Taylor series
// keeps all digits, and erfc(x) = 1 - erf(x)
template <typename T, int N>
constexpr auto erf_small_x(const boost::int128::uint128_t& sig, const int e, const bool complement, const bool neg) noexcept -> exp_detail::exp_digits<N>
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int scale10 {fx_scale10(N)};
    using digits_type = fx_digits_t<N>;
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};

    const auto ax {fx_from_decimal<N>(sig, -e)};
    const auto s {erf_series(fx_mul(ax, ax))};
    if (!complement)
    {
        // m has scale10 + 1 digits where s < 1, else scale10, and one more where m * s has a digit fewer
        int len {s.w[N - 1] >= (UINT64_C(1) << 62U) ? scale10 : scale10 + 1};
        const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(len - p))};
        auto d {fx_scale(m, s)};
        if (d < unit)
        {
            ++len;
            d = fx_scale(m * 10U, s);
        }
        return {static_cast<digits_type>(d), e - (len - p), exp_detail::exp_kind::digits};
    }
    const auto er {fx_mul(ax, s)};
    if (neg)
    {
        return {static_cast<digits_type>(fx_scale(unit * 2U, fx_shr(fx_add(fx_one<N>(), er), 1U))), -scale10, exp_detail::exp_kind::digits};
    }
    return {static_cast<digits_type>(fx_scale(unit * 10U, fx_sub(fx_one<N>(), er))), -scale10 - 1, exp_detail::exp_kind::digits};
}

// x^2 = whole^2 + 2 whole frac + frac^2, and exp(-x^2) = 10^(-n - 1) * 10^(1 - f) for x^2 log10(e) = n + f
template <int N>
constexpr auto erf_square_log10(const std::uint64_t whole, const fx<N>& frac) noexcept -> exp_detail::exp_parts<N>
{
    std::uint64_t sq[static_cast<std::size_t>(N + 1)] {};
    sq[N] = fx_mul_add(fx_make_span<N>(sq), fx_make_span(frac.w), 2U * whole);
    const auto f2 {fx_mul(frac, frac)};
    sq[N] += fx_add_to(fx_make_span<N>(sq), fx_make_span(f2.w));
    return exp_detail::exp_log10(whole * whole + ((sq[N] << 2U) | (sq[N - 1] >> 62U)), fx_frac<N>(sq), fx_table::log10_e);
}

// For decimal64_t and x >= 1, erfc(x) in one word with n + f from two words, rounded once in out away from a rounding
// boundary: 10^(1 - f) = 2^a m is below 8.8 units of 2^-62 off, so with g and q below 17 for x < 10, else 23
template <typename T>
constexpr auto erfc_above_one_try(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole, T& out) noexcept -> bool
{
    using sig_type = boost::int128::uint128_t;

    const auto frac {exp_detail::exp_frac<2>(sig, k, whole)};
    const auto y {erf_square_log10(whole, frac)};
    const auto r10 {exp_detail::exp10_frac<1>(fx_top<1>(fx_sub(fx_one<2>(), y.f)))};
    const auto g {whole < 10U ? erf_g_float<2>(whole, frac) : erf_g_asymptotic_float<2>(whole, frac)};

    // erfc(x) = q 2^-(61 + l - a) 10^(-n - 1) with q = m g / 2^63, and d its digits times 10^19, or 10^20 below 10^18
    const auto q {((static_cast<sig_type>(r10.m.w[0]) * g.m) >> 63U).low};
    const auto t {static_cast<sig_type>(q) * pow10(static_cast<sig_type>(19))};
    const auto sh {static_cast<unsigned>(61 + g.l) - r10.a};
    int e10 {-y.exp10 - 20};
    auto d {t >> sh};
    if (d < pow10(static_cast<sig_type>(18)))
    {
        d = ((t >> 4U) * 10U) >> (sh - 4U);
        --e10;
    }
    auto w {((d * (whole < 10U ? 17U : 23U)) >> 62U).low + 2U};
    if (d >= pow10(static_cast<sig_type>(19)))
    {
        d /= 10U;
        ++e10;
        w = w / 10U + 2U;
    }
    return fx_try_round(d.low, e10, false, w, out);
}

// For decimal128_t and x >= 1, the digits D of erfc(x) in two words with n + f from three words: m g / 2^252 =
// q 2^(2 - lz) with q in [2, 4), and D = q C with C = 10^P 2^e in (2^126 / 10, 2^126]
constexpr auto erfc_above_one_digits2(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole) noexcept -> exp_detail::exp_digits<3>
{
    using sig_type = boost::int128::uint128_t;

    const auto frac3 {exp_detail::exp_frac<3>(sig, k, whole)};
    const auto y {erf_square_log10(whole, frac3)};
    const auto r10 {exp_detail::exp10_frac<2, true>(fx_top<2>(fx_sub(fx_one<3>(), y.f)))};
    const auto frac {fx_top<2>(frac3)};
    const auto g {whole < 12U ? erf_g_float<3>(whole, frac) : erf_g_asymptotic_float<3>(whole, frac)};
    const auto pr {umul256(erf_mant(r10.m), g.m)};
    const auto lz {static_cast<unsigned>(countl_zero(pr[3]))};
    const fx<2> q {{(pr[2] << lz) | (pr[1] >> (64U - lz)), (pr[3] << lz) | (pr[2] >> (64U - lz))}};
    const int e {2 - static_cast<int>(lz) + static_cast<int>(r10.a) - g.l};
    const int p {((126 - e) * 78913) >> 18};
    const auto c {e < 0 ? (pow10(static_cast<sig_type>(p - 3)) >> -e) * 1000U : pow10(static_cast<sig_type>(p)) << e};
    return {fx_scale(c, q), -y.exp10 - 1 - p, exp_detail::exp_kind::digits};
}

// g in N words, and for the first try of decimal128_t from its floating g, below 8.34 units of 2^-126 off relative
// and 1 more unit off
template <int N>
constexpr auto erf_g_fixed(const std::uint64_t whole, const fx<N>& frac, std::false_type) noexcept -> fx<N>
{
    // The end of the last interval of erf_interval_of, where g becomes asymptotic
    constexpr std::uint64_t asymptotic {N == 1 ? 8U : N == 2 ? 10U : 16U};
    return whole < asymptotic ? erf_g_table(whole, frac) : erf_g_asymptotic(whole, frac);
}

constexpr auto erf_g_fixed(const std::uint64_t whole, const fx<2>& frac, std::true_type) noexcept -> fx<2>
{
    const auto g {erf_g_float<3>(whole, frac)};
    const auto m {g.m >> g.l};
    return fx<2> {{m.low, m.high}};
}

// erf(x), or erfc(x) with complement, in N words for |x| = whole + frac >= 1: erfc(|x|) = exp(-x^2) * g(|x|), erf(x) =
// 1 - erfc(x) and erfc(-|x|) = 2 - erfc(|x|). Kind below_one where erfc(|x|) cannot change the rounding of 1 or 2.
template <typename T, int N>
constexpr auto erf_above_one(const boost::int128::uint128_t& sig, const int k, const std::uint64_t whole, const bool complement, const bool neg) noexcept -> exp_detail::exp_digits<N>
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int scale10 {fx_scale10(N)};
    using digits_type = fx_digits_t<N>;
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};
    const auto frac {exp_detail::exp_frac<N>(sig, k, whole)};
    const auto y {erf_square_log10(whole, frac)};
    int exp10 {-y.exp10 - 1};
    // A first try in two words for decimal128_t takes the precise series of exp(-z), with an error below 2^-124
    const auto r10 {exp_detail::exp10_frac<N, (N == 2 && fx_words(p) == 3)>(fx_sub(fx_one<N>(), y.f))};
    const auto g {erf_g_fixed(whole, frac, std::integral_constant<bool, N == 2 && fx_words(p) == 3> {})};

    // erfc(|x|) = 16 h * 10^exp10 with h in [1/16, 10/16)
    const auto mg {fx_mul(r10.m, g)};
    auto h {fx_shr(mg, 4U - r10.a)};
    while (h.w[N - 1] < (UINT64_C(1) << 58U))
    {
        h = fx_mul_small(h, 10U);
        --exp10;
    }
    if (complement && !neg)
    {
        return {static_cast<digits_type>(fx_scale(unit * 16U, h)), exp10 - scale10, exp_detail::exp_kind::digits};
    }

    // erfc(|x|) in fixed point. Below 10^-(p + 1) it cannot change the rounding of 1 or 2.
    if (exp10 < -(p + 1))
    {
        return {0U, 0, exp_detail::exp_kind::below_one};
    }
    // erfc(|x|) = 2^a m g * 10^(-y.exp10 - 1) in units of 10^-scale10, or of 10^-(scale10 + 1) for erf(x). It is not
    // an integer, so the floor of 2 - erfc or 1 - erfc is 1 below the difference with the floor of erfc.
    const int s {complement ? scale10 - 1 - y.exp10 : scale10 - y.exp10};
    const auto c {static_cast<digits_type>(fx_scale(pow10(static_cast<boost::int128::uint128_t>(s)) << r10.a, mg))};
    if (complement)
    {
        return {static_cast<digits_type>(unit * 2U) - c - 1U, -scale10, exp_detail::exp_kind::digits};
    }
    return {static_cast<digits_type>(unit * 10U) - c - 1U, -scale10 - 1, exp_detail::exp_kind::digits};
}

// A value in (2 - 10^-scale10, 2), rounded once in the current mode
template <typename T, typename digits_type>
constexpr auto just_below_two(const int scale10) noexcept -> T
{
    return fx_round<T>(static_cast<digits_type>(pow10(static_cast<boost::int128::uint128_t>(scale10)) * 2U - 1U), -scale10, false);
}

// erf(x), or erfc(x) with complement, for a finite x != 0, rounded once in the current mode. For |x| >= 1, erfc(|x|) =
// exp(-x^2) * g(|x|), and erf(x) = 1 - erfc(x) and erfc(-|x|) = 2 - erfc(|x|) in fixed point.
template <typename T>
constexpr auto erf_fixed(const T x, const bool complement) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};
    const int k {-e};
    const auto unit {pow10(static_cast<boost::int128::uint128_t>(scale10))};

    if (k >= p)
    {
        // For |x| < 10^-(p + 1), erfc(x) rounds as a value next to 1
        if (complement && k - p >= p + 1)
        {
            return neg ? fx_round<T>(static_cast<digits_type>(unit), -scale10, false) : exp_detail::just_below_one<T, digits_type>(scale10, false);
        }
        BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
        {
            // One word: the sum is below 22, 3.5 or 5.6 units of 2^-62 off for x < 1, 0.7 or 0.86, so d1 * 2^-57, 2^-60
            // or 2^-59 + 2 units of d1 for erf(x), and 54, 11 or 15 units for erfc(x) with x > 0, or 8 with x < 0
            const auto r1 {erf_small_x<T, 1>(sig, e, complement, neg)};
            const auto tenth {pow10(static_cast<boost::int128::uint128_t>(p - 1))};
            const bool below_07 {k > p || sig < tenth * 7U};
            const bool below_086 {below_07 || sig * 10U < tenth * 86U};
            const auto w {complement ? (neg ? 8U : below_07 ? 11U : below_086 ? 15U : 54U) : (r1.digits >> (below_07 ? 60U : below_086 ? 59U : 57U)) + 2U};
            T r {};
            if (fx_try_round(r1.digits, r1.exp10, !complement && neg, w, r))
            {
                return r;
            }
        }
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // A first try in two words: the error of the sum is below 210 units of 2^-126, so below d2 * 2^-118 + 2
            // units of d2 for erf(x), and 250 units for erfc(x)
            const auto r2 {erf_small_x<T, 2>(sig, e, complement, neg)};
            if (!fx_near_boundary(r2.digits, complement ? 250U : static_cast<std::uint64_t>(r2.digits >> 118U) + 2U))
            {
                return fx_round<T>(r2.digits, r2.exp10, !complement && neg);
            }
        }
        const auto r {erf_small_x<T, words>(sig, e, complement, neg)};
        return fx_round<T>(r.digits, r.exp10, !complement && neg);
    }

    // Above these bounds, erfc(x) is below half of the smallest subnormal value, and erf(x) rounds as a value next to 1
    const auto wide_whole {exp_detail::exp_whole(sig, k)};
    constexpr std::uint64_t tiny_bound {words == 1 ? 16U : words == 2 ? 31U : 120U};
    constexpr std::uint64_t one_bound {words == 1 ? 7U : 10U};
    if (wide_whole >= (complement && !neg ? tiny_bound : one_bound))
    {
        if (!complement)
        {
            return exp_detail::just_below_one<T, digits_type>(scale10, neg);
        }
        return neg ? just_below_two<T, digits_type>(scale10) : T {1U, -2 * std::numeric_limits<T>::max_exponent10};
    }
    const auto whole {static_cast<std::uint64_t>(wide_whole)};
    BOOST_DECIMAL_IF_CONSTEXPR (words == 2)
    {
        // For erf(x) and erfc(-|x|), a first try in one word: the digits c of erfc(|x|) are below c * 104 / 2^62 + 4
        // units off. Where it is too small in one word, it is below 10^-(p + 1) * (1 + 2^-50).
        if (!complement || neg)
        {
            const auto r1 {erf_above_one<T, 1>(sig, k, whole, complement, neg)};
            if (r1.kind != exp_detail::exp_kind::digits)
            {
                return complement ? just_below_two<T, digits_type>(scale10) : exp_detail::just_below_one<T, digits_type>(scale10, neg);
            }
            T r {};
            const auto c {(complement ? UINT64_C(2000000000000000000) : UINT64_C(10000000000000000000)) - r1.digits};
            if (fx_try_round(r1.digits, r1.exp10, !complement && neg, ((static_cast<boost::int128::uint128_t>(c) * 104U) >> 62U).low + 4U, r))
            {
                return r;
            }
        }
        else
        {
            T r {};
            if (erfc_above_one_try(sig, k, whole, r))
            {
                return r;
            }
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        // The same in two words: the error of erfc(|x|) is below 6.7 units of 2^-126, so below 10 or 3 units of d2
        if (!complement || neg)
        {
            const auto r2 {erf_above_one<T, 2>(sig, k, whole, complement, neg)};
            if (r2.kind != exp_detail::exp_kind::digits)
            {
                return complement ? just_below_two<T, digits_type>(scale10) : exp_detail::just_below_one<T, digits_type>(scale10, neg);
            }
            if (!fx_near_boundary(r2.digits, complement ? 3U : 10U))
            {
                return fx_round<T>(r2.digits, r2.exp10, !complement && neg);
            }
        }
        else
        {
            // For erfc(x), D is below 32 or 67 units of 2^-126 off for x < 12 or above, and the result is not subnormal
            const auto r2 {erfc_above_one_digits2(sig, k, whole)};
            const std::uint64_t err {whole < 12U ? 32U : 67U};
            if (r2.exp10 + 36 >= std::numeric_limits<T>::min_exponent10 - 1 && !fx_near_boundary(r2.digits, static_cast<std::uint64_t>(((r2.digits >> 64U) * err) >> 62U) + 2U))
            {
                return fx_round<T>(r2.digits, r2.exp10, false);
            }
        }
    }
    const auto r {erf_above_one<T, words>(sig, k, whole, complement, neg)};
    if (r.kind != exp_detail::exp_kind::digits)
    {
        return complement ? just_below_two<T, digits_type>(scale10) : exp_detail::just_below_one<T, digits_type>(scale10, neg);
    }
    return fx_round<T>(r.digits, r.exp10, !complement && neg);
}

} // namespace erf_detail

template <typename T>
constexpr auto erf_impl(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // Edge cases
    const auto fp {fpclassify(z)};

    if (fp == FP_ZERO || fp == FP_NAN)
    {
        return issignaling(z) ? nan_conversion(z) : z;
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fp == FP_INFINITE)
    {
        return z < T{0} ? T{-1} : T{1};
    }
    #endif

    return erf_detail::erf_fixed(z, false);
}

template <typename T>
constexpr auto erfc_impl(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    // Edge cases
    const auto fp {fpclassify(z)};

    if (fp == FP_NAN)
    {
        return issignaling(z) ? nan_conversion(z) : z;
    }
    else if (fp == FP_ZERO)
    {
        return T{1};
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fp == FP_INFINITE)
    {
        return z < T{0} ? T{2} : T{0};
    }
    #endif

    return erf_detail::erf_fixed(z, true);
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto erf(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::erf_impl(static_cast<evaluation_type>(z)));
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto erfc(const T z) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::erfc_impl(static_cast<evaluation_type>(z)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_ERF_HPP
