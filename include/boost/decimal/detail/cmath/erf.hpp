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
#include <boost/decimal/detail/check_non_finite.hpp>
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

// The words of g as a mantissa
constexpr auto erf_mant(const fx<1>& g) noexcept -> std::uint64_t { return g.w[0]; }
constexpr auto erf_mant(const fx<2>& g) noexcept -> boost::int128::uint128_t { return {g.w[1], g.w[0]}; }

// For the first M coefficients of the table of N words, m in [2^(64N - 67), 2^(64N - 65)) with each l below the next
// one, the shift l_(k+1) - l_k - 1 that takes 2 u s_(k+1) times 2^-(64N - 64) to the scale of c_k, l of each first
// coefficient, and for |u| < (j + 1) / 8 the last term and for three words the first term below 2^-64 of the first
template <int N, std::size_t M>
struct erf_float_table
{
    int l[4];
    unsigned char last[4][8];
    unsigned char split[4][8];
    unsigned char shift[M];
};

template <int N>
using erf_float_table_t = erf_float_table<N, N == 2 ? 134U : 199U>;

// m of each coefficient of the table, in a struct of its own: the smaller members after 128-bit words would pad it
template <int N, std::size_t M>
struct erf_float_mant
{
    erf_mant_t<N> m[M];
};

template <int N>
using erf_float_mant_t = erf_float_mant<N, N == 2 ? 134U : 199U>;

// The tables of decimal64_t on [1, 10] and of decimal128_t on [1, 12], from erfc_d64 and erfc_d128. A test in
// github_issue_erf.cpp makes them again from the coefficients.
constexpr auto erf_float_mant_for(std::integral_constant<int, 2>) noexcept -> erf_float_mant_t<2>
{
    return {
        {
            UINT64_C(0x52536BFE8AB5A139), UINT64_C(0x53C66211132BDADB), UINT64_C(0x4DFA89C6F152F454),
            UINT64_C(0x4371474A5079B449), UINT64_C(0x6D96513CD3DEB477), UINT64_C(0x544FB5D8034701F5),
            UINT64_C(0x7B9C179CD91932A7), UINT64_C(0x56C2AFF9968CB47D), UINT64_C(0x75142743505F5692),
            UINT64_C(0x4C303B1D0EF378EF), UINT64_C(0x5FE65E48D50DF963), UINT64_C(0x75055C622BCC4FBF),
            UINT64_C(0x455B24DAB0D724A3), UINT64_C(0x50013BFDBA0D7CF4), UINT64_C(0x59F3D7ACAB970FCF),
            UINT64_C(0x62B706A436D5A9B0), UINT64_C(0x69DE6A5E2E0B94D5), UINT64_C(0x6F15DF8938A5755E),
            UINT64_C(0x72277D6866741649), UINT64_C(0x72FE6CA4AFB19AEA), UINT64_C(0x71A68776BD036D54),
            UINT64_C(0x6E491B5D78F554F4), UINT64_C(0x69277087E135E7D7), UINT64_C(0x6293DFEF2AFE2D42),
            UINT64_C(0x5AEA96FE72B114D8), UINT64_C(0x528C2C8658DB3500), UINT64_C(0x49CEF956C348CC00),
            UINT64_C(0x40C8BCE5E5459000), UINT64_C(0x708FA43BAFD80000), UINT64_C(0x6966D3BB38D80000),
            UINT64_C(0x58E4321C8BC00000), UINT64_C(0x5BA609F48A74B337), UINT64_C(0x6F5AB93452C5BB9C),
            UINT64_C(0x410FF86A63033412), UINT64_C(0x4964186EA15A4804), UINT64_C(0x502730BB4FFBF078),
            UINT64_C(0x54F7B27422554AA5), UINT64_C(0x579CE4C1213FD7EE), UINT64_C(0x5809440F0661DE9C),
            UINT64_C(0x5657C6D771D9C3E5), UINT64_C(0x52C518A191C1F0AD), UINT64_C(0x4DA6412DA81A97A4),
            UINT64_C(0x475E45155085CE63), UINT64_C(0x405423A47490A23E), UINT64_C(0x71D481BAF1B7C84B),
            UINT64_C(0x62EFC48ABFA349BC), UINT64_C(0x5489FC6C8CC71D66), UINT64_C(0x47108E72AC1A98C7),
            UINT64_C(0x759C7E72FED960D3), UINT64_C(0x5FDBFC8C9EE76C78), UINT64_C(0x4CFF2A2B0081487B),
            UINT64_C(0x79F5295A9A914E6B), UINT64_C(0x5F46CB276125267A), UINT64_C(0x497475D2B632177A),
            UINT64_C(0x6FD1F190CC99D7BE), UINT64_C(0x540F479D19A11D67), UINT64_C(0x7CDD62C5E1ADFE7E),
            UINT64_C(0x5BA82A5AD151267B), UINT64_C(0x428414D66CD4D3D6), UINT64_C(0x5F75AFEAAD22917E),
            UINT64_C(0x43C57C3583061D6C), UINT64_C(0x5F7C2947415A7AC0), UINT64_C(0x426080F2F9217380),
            UINT64_C(0x57EE065EF8662000), UINT64_C(0x780824FC11894000), UINT64_C(0x6BC0EF4DF3EC8000),
            UINT64_C(0x47BB1B2FD1D80000), UINT64_C(0x5F00D212275C85CA), UINT64_C(0x7B5FD9B67427F761),
            UINT64_C(0x4F1DD92A07D1882C), UINT64_C(0x6441A9F35E78ED2F), UINT64_C(0x7D9336701EDAC927),
            UINT64_C(0x4DC341905DD7841B), UINT64_C(0x5F4380FB6BA5753A), UINT64_C(0x737933B338D0DED3),
            UINT64_C(0x4544D1A8084670E4), UINT64_C(0x5246EC3302DDEE21), UINT64_C(0x60C8067D1B331D5F),
            UINT64_C(0x70C4B26AD6E4012A), UINT64_C(0x4118039A8AD5ED45), UINT64_C(0x4A797B2296D440BA),
            UINT64_C(0x5475F92A0FAB5831), UINT64_C(0x5EF779B391B9CA77), UINT64_C(0x69E2EF1B12FCEE4B),
            UINT64_C(0x7518A872ECFBBEDB), UINT64_C(0x403A748321505D3F), UINT64_C(0x45E85680C94761E8),
            UINT64_C(0x4B814A7691E529A7), UINT64_C(0x50EFF7CDE5131C64), UINT64_C(0x561EEC4A62C699FE),
            UINT64_C(0x5AF92FDB8ED431FE), UINT64_C(0x5F6AD6FBEED550EC), UINT64_C(0x63619054E1F5FC55),
            UINT64_C(0x66CD2CFA1D4E72FE), UINT64_C(0x699FE1B2AD9A4302), UINT64_C(0x6BCE050302585AA8),
            UINT64_C(0x6D51B5DF85D0351C), UINT64_C(0x6E32DA10258470C2), UINT64_C(0x6E550A616113257C),
            UINT64_C(0x6D32E0A6DD4C6BE7), UINT64_C(0x6C0606DC4FFBE0DC), UINT64_C(0x6F051CC21995C6D8),
            UINT64_C(0x6C560561F276F7E0), UINT64_C(0x4E0AF103C7EBC800), UINT64_C(0x4BFB6E294753D000),
            UINT64_C(0x5395C8EBD4D3FC00), UINT64_C(0x4F57710B95D28000), UINT64_C(0x7F9B31283E54E1C9),
            UINT64_C(0x70119A675DF54575), UINT64_C(0x61D8DCFD283CE95F), UINT64_C(0x54F04C6D72CA7C64),
            UINT64_C(0x4950E040E323ADCC), UINT64_C(0x7DDCD4170D44F94A), UINT64_C(0x6B717DFC9CD04303),
            UINT64_C(0x5B3998C6A1E7D9A8), UINT64_C(0x4D0B21D26AB6E254), UINT64_C(0x40B9EDE2533A641B),
            UINT64_C(0x6C3238614EEB5F16), UINT64_C(0x59F88F995945AD07), UINT64_C(0x4A70F3886FDFA29A),
            UINT64_C(0x7A93827356F678BF), UINT64_C(0x646D4BF7ED4AF008), UINT64_C(0x51E290F282F607C9),
            UINT64_C(0x4273473898B68FE3), UINT64_C(0x6B5884E4B5FD38EA), UINT64_C(0x564E06D7A21A01DE),
            UINT64_C(0x4512722995DE87A0), UINT64_C(0x6E10374F76D32100), UINT64_C(0x574DD5C83FFB2E00),
            UINT64_C(0x44F1E4AE2BDA9000), UINT64_C(0x6C3648DDCAB10000), UINT64_C(0x54B9695435580000),
            UINT64_C(0x46835319A5800000), UINT64_C(0x6D6D8E76E8000000)
        }
    };
}

constexpr auto erf_float_table_for(std::integral_constant<int, 2>) noexcept -> erf_float_table_t<2>
{
    return {
        {2, 3, 4, 5},
        {{12, 14, 16, 18, 20, 21, 23, 24}, {44, 47, 50, 52, 54, 56, 59, 61},
         {81, 84, 87, 90, 93, 96, 99, 102}, {117, 119, 121, 122, 124, 125, 126, 127}},
        {{0, 0, 0, 0, 0, 0, 0, 0}, {31, 31, 31, 31, 31, 31, 31, 31},
         {67, 67, 67, 67, 67, 67, 67, 67}, {107, 107, 107, 107, 107, 107, 107, 107}},
        {
            1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 0, 1, 0, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1, 0, 1, 0, 1, 1, 0,
            1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 2,
            2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 3, 0
        }
    };
}

constexpr auto erf_float_mant_for(std::integral_constant<int, 3>) noexcept -> erf_float_mant_t<3>
{
    return {
        {
            boost::int128::uint128_t {UINT64_C(0x52536BFE8AB5A139), UINT64_C(0xBDF0F0D8E20F3F6B)},
            boost::int128::uint128_t {UINT64_C(0x53C66211132BDADB), UINT64_C(0x53AC84E8639F908D)},
            boost::int128::uint128_t {UINT64_C(0x4DFA89C6F152F454), UINT64_C(0xFCBE34AA5D5E4C03)},
            boost::int128::uint128_t {UINT64_C(0x4371474A5079B449), UINT64_C(0x904FA3C1A4423EC7)},
            boost::int128::uint128_t {UINT64_C(0x6D96513CD3DEB477), UINT64_C(0x4209E76488B273B9)},
            boost::int128::uint128_t {UINT64_C(0x544FB5D8034701F5), UINT64_C(0xE4F2F65969977138)},
            boost::int128::uint128_t {UINT64_C(0x7B9C179CD91932A7), UINT64_C(0x90DF2703EEACFCD2)},
            boost::int128::uint128_t {UINT64_C(0x56C2AFF9968CB47D), UINT64_C(0xD8477F5B49582A15)},
            boost::int128::uint128_t {UINT64_C(0x75142743505F5692), UINT64_C(0x5D530EFEEF55BA85)},
            boost::int128::uint128_t {UINT64_C(0x4C303B1D0EF378EF), UINT64_C(0xECF266A2FF4724DB)},
            boost::int128::uint128_t {UINT64_C(0x5FE65E48D50DF963), UINT64_C(0xDF62033AB30099D5)},
            boost::int128::uint128_t {UINT64_C(0x75055C622BCC4FC1), UINT64_C(0x6084C9F73A9B847A)},
            boost::int128::uint128_t {UINT64_C(0x455B24DAB0D724A4), UINT64_C(0x79959F52A6AF63E1)},
            boost::int128::uint128_t {UINT64_C(0x50013BFDBA0D7B69), UINT64_C(0xDF2D030B65614DDA)},
            boost::int128::uint128_t {UINT64_C(0x59F3D7ACAB970E8A), UINT64_C(0xE8E644FC53F20299)},
            boost::int128::uint128_t {UINT64_C(0x62B706A436D71DC8), UINT64_C(0x10B8E053C0B4EE13)},
            boost::int128::uint128_t {UINT64_C(0x69DE6A5E2E0CC6BF), UINT64_C(0xC541E1A88B4090AA)},
            boost::int128::uint128_t {UINT64_C(0x6F15DF8937AA061E), UINT64_C(0x670F72E3AD964769)},
            boost::int128::uint128_t {UINT64_C(0x72277D6865A55D40), UINT64_C(0xA832079AFA6D40FB)},
            boost::int128::uint128_t {UINT64_C(0x72FE6CA52B06C1EF), UINT64_C(0x6BB07EFB9A19DF4B)},
            boost::int128::uint128_t {UINT64_C(0x71A68777226AEDA4), UINT64_C(0x993219539B48E3CF)},
            boost::int128::uint128_t {UINT64_C(0x6E491B31833AE8DE), UINT64_C(0x1528D5BEC0BB0031)},
            boost::int128::uint128_t {UINT64_C(0x69277063BBEFC35F), UINT64_C(0x184C8B6C5FA7D35D)},
            boost::int128::uint128_t {UINT64_C(0x6293EB323E442BEB), UINT64_C(0xE4A9D5586DB04FAD)},
            boost::int128::uint128_t {UINT64_C(0x5AEAA04130C843CF), UINT64_C(0x1C382AC7C10D6890)},
            boost::int128::uint128_t {UINT64_C(0x528A27CF1990C24D), UINT64_C(0x78B41203EA33DEC0)},
            boost::int128::uint128_t {UINT64_C(0x49CD5066CDC76444), UINT64_C(0xF4854A4C9F0DFC00)},
            boost::int128::uint128_t {UINT64_C(0x41061E68F23D061D), UINT64_C(0x896EC20F3311E000)},
            boost::int128::uint128_t {UINT64_C(0x70F49DB598FA5AE7), UINT64_C(0xB6484E1A1D950000)},
            boost::int128::uint128_t {UINT64_C(0x60C2C40F3A7083EB), UINT64_C(0x7FEA63D8D9180000)},
            boost::int128::uint128_t {UINT64_C(0x51C815886BF4C1BC), UINT64_C(0x3BE15495C6C00000)},
            boost::int128::uint128_t {UINT64_C(0x4439883D74556416), UINT64_C(0x6DDDBE5C6C000000)},
            boost::int128::uint128_t {UINT64_C(0x706504E18A2A50B5), UINT64_C(0x883E7164E0000000)},
            boost::int128::uint128_t {UINT64_C(0x5B71A1D51D0CA333), UINT64_C(0xCC2C971200000000)},
            boost::int128::uint128_t {UINT64_C(0x4982443270493816), UINT64_C(0x11C5630800000000)},
            boost::int128::uint128_t {UINT64_C(0x74CE27F5953D033C), UINT64_C(0xB1BE5D8000000000)},
            boost::int128::uint128_t {UINT64_C(0x5BBF3B6ADAEFE385), UINT64_C(0x459AD00000000000)},
            boost::int128::uint128_t {UINT64_C(0x4746127EB9BE2A4A), UINT64_C(0xCF61800000000000)},
            boost::int128::uint128_t {UINT64_C(0x6D863AA364280B8D), UINT64_C(0x9870000000000000)},
            boost::int128::uint128_t {UINT64_C(0x52E6B78A121D29C9), UINT64_C(0xFA30000000000000)},
            boost::int128::uint128_t {UINT64_C(0x7CC2DE8757F77AB4), UINT64_C(0xE600000000000000)},
            boost::int128::uint128_t {UINT64_C(0x662212C9A4AC7816), UINT64_C(0x1000000000000000)},
            boost::int128::uint128_t {UINT64_C(0x4B2D854F453DA4AE), UINT64_C(0x8000000000000000)},
            boost::int128::uint128_t {UINT64_C(0x5BA609F48A74B337), UINT64_C(0xC5610C765368D086)},
            boost::int128::uint128_t {UINT64_C(0x6F5AB93452C5BB9C), UINT64_C(0x9C2F7CDEEE18A79C)},
            boost::int128::uint128_t {UINT64_C(0x410FF86A63033412), UINT64_C(0x81EB767906B2968C)},
            boost::int128::uint128_t {UINT64_C(0x4964186EA15A4804), UINT64_C(0x86F6E63215997435)},
            boost::int128::uint128_t {UINT64_C(0x502730BB4FFBF078), UINT64_C(0xE5924E9BB3FBFB23)},
            boost::int128::uint128_t {UINT64_C(0x54F7B27422554AA5), UINT64_C(0x783AAE552A4FCBDD)},
            boost::int128::uint128_t {UINT64_C(0x579CE4C1213FD7EE), UINT64_C(0xE776E9E9C1560BF1)},
            boost::int128::uint128_t {UINT64_C(0x5809440F0661DE9C), UINT64_C(0xC2E2441ABCD80D49)},
            boost::int128::uint128_t {UINT64_C(0x5657C6D771D9C3E5), UINT64_C(0x5534DB56CED007F8)},
            boost::int128::uint128_t {UINT64_C(0x52C518A191C1F0AD), UINT64_C(0x5FECE2223F0E5310)},
            boost::int128::uint128_t {UINT64_C(0x4DA6412DA81A97A4), UINT64_C(0x2A709F29FE7739A0)},
            boost::int128::uint128_t {UINT64_C(0x475E45155085CE63), UINT64_C(0xA32FC2368C3C8D4C)},
            boost::int128::uint128_t {UINT64_C(0x405423A47490A243), UINT64_C(0xD577795838E0B70F)},
            boost::int128::uint128_t {UINT64_C(0x71D481BAF1B7C852), UINT64_C(0x85595E8BB98750C5)},
            boost::int128::uint128_t {UINT64_C(0x62EFC48ABFA3453A), UINT64_C(0xEB3FD1C82B1A06A8)},
            boost::int128::uint128_t {UINT64_C(0x5489FC6C8CC71A73), UINT64_C(0xE8588B1CA032457A)},
            boost::int128::uint128_t {UINT64_C(0x47108E72AC1BE2C7), UINT64_C(0xF9FAD300804F08B6)},
            boost::int128::uint128_t {UINT64_C(0x759C7E72FEDB1070), UINT64_C(0xA8E37C218DB5DDC4)},
            boost::int128::uint128_t {UINT64_C(0x5FDBFC8C9E5B8230), UINT64_C(0xD0854C4EFEE29F45)},
            boost::int128::uint128_t {UINT64_C(0x4CFF2A2B0025C814), UINT64_C(0x4AF1864014908582)},
            boost::int128::uint128_t {UINT64_C(0x79F5295AC7308D38), UINT64_C(0x4DDF1892A4E8C350)},
            boost::int128::uint128_t {UINT64_C(0x5F46CB277E540983), UINT64_C(0x97B69E2CB003C103)},
            boost::int128::uint128_t {UINT64_C(0x497475CD52197D36), UINT64_C(0xAF6E1C5987C63A58)},
            boost::int128::uint128_t {UINT64_C(0x6FD1F189BF507955), UINT64_C(0x063B8150D2A8D812)},
            boost::int128::uint128_t {UINT64_C(0x540F4899C648D494), UINT64_C(0x726FE3F1CDF5A99F)},
            boost::int128::uint128_t {UINT64_C(0x7CDD64106D0EFAA5), UINT64_C(0x09C6D9ED5ED58258)},
            boost::int128::uint128_t {UINT64_C(0x5BA80792BE7E0F5C), UINT64_C(0x9396404FED1D5E8E)},
            boost::int128::uint128_t {UINT64_C(0x4283FE15C4F7072F), UINT64_C(0xCA6690D3D82C4A35)},
            boost::int128::uint128_t {UINT64_C(0x5F79392A189D7D6B), UINT64_C(0xCBE1B68498EA180C)},
            boost::int128::uint128_t {UINT64_C(0x43C7CC7722B9ACBF), UINT64_C(0xFC10A0088F01F6C0)},
            boost::int128::uint128_t {UINT64_C(0x5F3A6601967D01A6), UINT64_C(0xFAC3B91C1D3D7980)},
            boost::int128::uint128_t {UINT64_C(0x4235782E88C66CA8), UINT64_C(0xAEEF38440F5FA300)},
            boost::int128::uint128_t {UINT64_C(0x5B2497BD63198EAF), UINT64_C(0x39446A62390A0800)},
            boost::int128::uint128_t {UINT64_C(0x7C3D147A452628E0), UINT64_C(0xF8837E305CF0C000)},
            boost::int128::uint128_t {UINT64_C(0x53DD87A43B99128C), UINT64_C(0xDD9790F965628000)},
            boost::int128::uint128_t {UINT64_C(0x7029EF12EBAC1477), UINT64_C(0x0286DD1788400000)},
            boost::int128::uint128_t {UINT64_C(0x4A516E46413778BF), UINT64_C(0x85F9EBCFE0100000)},
            boost::int128::uint128_t {UINT64_C(0x6199B4B4AC372A4F), UINT64_C(0xFB30ED18DB000000)},
            boost::int128::uint128_t {UINT64_C(0x7F0CC12233463B32), UINT64_C(0xAE7C5E0B2C000000)},
            boost::int128::uint128_t {UINT64_C(0x51FAEE9B5C3B32ED), UINT64_C(0x757338D9B0000000)},
            boost::int128::uint128_t {UINT64_C(0x68E79F83D8768C6F), UINT64_C(0x7030219300000000)},
            boost::int128::uint128_t {UINT64_C(0x42901D4168E5C871), UINT64_C(0xC0E242EE00000000)},
            boost::int128::uint128_t {UINT64_C(0x53C2D79933096A68), UINT64_C(0xC9FDD68000000000)},
            boost::int128::uint128_t {UINT64_C(0x68971F0BEDA5A80C), UINT64_C(0x1AE84F8000000000)},
            boost::int128::uint128_t {UINT64_C(0x410EA6E2A5EF27F7), UINT64_C(0x585F360000000000)},
            boost::int128::uint128_t {UINT64_C(0x4FF672B56972A989), UINT64_C(0x987D000000000000)},
            boost::int128::uint128_t {UINT64_C(0x5CD7A586A12A7C08), UINT64_C(0x16A7800000000000)},
            boost::int128::uint128_t {UINT64_C(0x709739B6DBBE59F0), UINT64_C(0x3780000000000000)},
            boost::int128::uint128_t {UINT64_C(0x5D5B72C8EEE0EE7B), UINT64_C(0xB120000000000000)},
            boost::int128::uint128_t {UINT64_C(0x6EFCA640A20D59D0), UINT64_C(0x3180000000000000)},
            boost::int128::uint128_t {UINT64_C(0x5F00D212275C85CA), UINT64_C(0x9CEBB9A1C4424B85)},
            boost::int128::uint128_t {UINT64_C(0x7B5FD9B67427F761), UINT64_C(0x9CBD033166D47ECA)},
            boost::int128::uint128_t {UINT64_C(0x4F1DD92A07D1882C), UINT64_C(0xEBBEE796E45D8DC8)},
            boost::int128::uint128_t {UINT64_C(0x6441A9F35E78ED2F), UINT64_C(0x97E351415B58D7D9)},
            boost::int128::uint128_t {UINT64_C(0x7D9336701EDAC927), UINT64_C(0x3A8D54BBFB5D7CA1)},
            boost::int128::uint128_t {UINT64_C(0x4DC341905DD7841B), UINT64_C(0xCC774904BDBC84CA)},
            boost::int128::uint128_t {UINT64_C(0x5F4380FB6BA5753A), UINT64_C(0xFE6EF789481CE4BC)},
            boost::int128::uint128_t {UINT64_C(0x737933B338D0DED3), UINT64_C(0x429362208BEAF222)},
            boost::int128::uint128_t {UINT64_C(0x4544D1A8084670E4), UINT64_C(0x64036F86F965791A)},
            boost::int128::uint128_t {UINT64_C(0x5246EC3302DDEE21), UINT64_C(0x3CD6F7B97EBE1DD3)},
            boost::int128::uint128_t {UINT64_C(0x60C8067D1B331D5F), UINT64_C(0x84E916315124DAC2)},
            boost::int128::uint128_t {UINT64_C(0x70C4B26AD6E40129), UINT64_C(0xDEEA80EEA194D01F)},
            boost::int128::uint128_t {UINT64_C(0x4118039A8AD5ED55), UINT64_C(0xF3DBC7F72A815804)},
            boost::int128::uint128_t {UINT64_C(0x4A797B2296D440C9), UINT64_C(0x8F668EB0A43AFB50)},
            boost::int128::uint128_t {UINT64_C(0x5475F92A0FAB4FF4), UINT64_C(0xDE3EBF8631870C2A)},
            boost::int128::uint128_t {UINT64_C(0x5EF779B391B9C2BA), UINT64_C(0xF50CC4ACBBA7F615)},
            boost::int128::uint128_t {UINT64_C(0x69E2EF1B12FFEF45), UINT64_C(0x33A960336C8A6081)},
            boost::int128::uint128_t {UINT64_C(0x7518A872ECFE9126), UINT64_C(0x11585CD725132547)},
            boost::int128::uint128_t {UINT64_C(0x403A748320E72502), UINT64_C(0x508B04A92F504B97)},
            boost::int128::uint128_t {UINT64_C(0x45E85680C8E48C3E), UINT64_C(0xA3791A7A67FDA70A)},
            boost::int128::uint128_t {UINT64_C(0x4B814A76A7E4B215), UINT64_C(0x5936A290F9487CC3)},
            boost::int128::uint128_t {UINT64_C(0x50EFF7CDF9BD1C89), UINT64_C(0xEAE17A9B3DAFD852)},
            boost::int128::uint128_t {UINT64_C(0x561EEC46D4E21B39), UINT64_C(0xA827E485ACFA55CC)},
            boost::int128::uint128_t {UINT64_C(0x5AF92FD83813A79F), UINT64_C(0xA40B187CA2EFC479)},
            boost::int128::uint128_t {UINT64_C(0x5F6AD76E39BC9D5F), UINT64_C(0xE3F40173231B7EF2)},
            boost::int128::uint128_t {UINT64_C(0x636190C04167644C), UINT64_C(0x38790625798EFB2A)},
            boost::int128::uint128_t {UINT64_C(0x66CD21CEEFACA80D), UINT64_C(0xD32AAA663158D585)},
            boost::int128::uint128_t {UINT64_C(0x699FD73466A0B65E), UINT64_C(0x59AF55E8A3A6185D)},
            boost::int128::uint128_t {UINT64_C(0x6BCEDD477A12DC5F), UINT64_C(0x2193D86DF310428A)},
            boost::int128::uint128_t {UINT64_C(0x6D528115280B2284), UINT64_C(0xD6FE2A07ADC79BEF)},
            boost::int128::uint128_t {UINT64_C(0x6E265753AAD6BA67), UINT64_C(0x69379D2D754E27EC)},
            boost::int128::uint128_t {UINT64_C(0x6E49489E9452BE83), UINT64_C(0xC36C151BF505376E)},
            boost::int128::uint128_t {UINT64_C(0x6DBD8372EE81690A), UINT64_C(0x4734164CAD892FEC)},
            boost::int128::uint128_t {UINT64_C(0x6C885575D5F6DA9E), UINT64_C(0xC894CD47BD1E2ABC)},
            boost::int128::uint128_t {UINT64_C(0x6AB1EE8212BCB2D1), UINT64_C(0x8849999D2747FF20)},
            boost::int128::uint128_t {UINT64_C(0x684510BF52996936), UINT64_C(0xA9C003543E321F40)},
            boost::int128::uint128_t {UINT64_C(0x654EB197D4DD2704), UINT64_C(0xF0C1BA7153B8CA00)},
            boost::int128::uint128_t {UINT64_C(0x61DD8FBCCF195DEA), UINT64_C(0xD87FA655B1834E00)},
            boost::int128::uint128_t {UINT64_C(0x5E01C29F18CA98E4), UINT64_C(0x601444DCB2405000)},
            boost::int128::uint128_t {UINT64_C(0x59CC47B1B1CF1454), UINT64_C(0x6470A801E9C66000)},
            boost::int128::uint128_t {UINT64_C(0x554E904AE8A5FD93), UINT64_C(0x8757F6D515A98000)},
            boost::int128::uint128_t {UINT64_C(0x509A134D17F282BC), UINT64_C(0xE0E5F58E8D880000)},
            boost::int128::uint128_t {UINT64_C(0x4BBFF9025A91224F), UINT64_C(0x1FD2E2A49D280000)},
            boost::int128::uint128_t {UINT64_C(0x46D0E0649CCE262A), UINT64_C(0x8B2453361A000000)},
            boost::int128::uint128_t {UINT64_C(0x41DBC1B26D98ADF1), UINT64_C(0xD3B0642243000000)},
            boost::int128::uint128_t {UINT64_C(0x79D93E6B78320E83), UINT64_C(0xB3B154D880000000)},
            boost::int128::uint128_t {UINT64_C(0x702AD8ADBBA37153), UINT64_C(0xB9AF9EC060000000)},
            boost::int128::uint128_t {UINT64_C(0x66E9276845AE2009), UINT64_C(0x70AFDA5C80000000)},
            boost::int128::uint128_t {UINT64_C(0x5DC9704F38BC8E09), UINT64_C(0xB3E42DA000000000)},
            boost::int128::uint128_t {UINT64_C(0x53D6F9CC72FBAB55), UINT64_C(0x47FA3CB000000000)},
            boost::int128::uint128_t {UINT64_C(0x4BB6C18FD61164BB), UINT64_C(0x9A35377000000000)},
            boost::int128::uint128_t {UINT64_C(0x4A8421014B0E288D), UINT64_C(0xBE6F590000000000)},
            boost::int128::uint128_t {UINT64_C(0x4276B7EEF8F6497F), UINT64_C(0xD510EC0000000000)},
            boost::int128::uint128_t {UINT64_C(0x455BA8EB63138500), UINT64_C(0xFE1FD80000000000)},
            boost::int128::uint128_t {UINT64_C(0x3E3941B2ED80E511), UINT64_C(0x9328600000000000)},
            boost::int128::uint128_t {UINT64_C(0x557212A0ABFB4A97), UINT64_C(0xF4F8400000000000)},
            boost::int128::uint128_t {UINT64_C(0x4A8C81732EDCD742), UINT64_C(0xF4CA000000000000)},
            boost::int128::uint128_t {UINT64_C(0x72FA0C7E65C48F09), UINT64_C(0x3F50ED58EB5A1E3E)},
            boost::int128::uint128_t {UINT64_C(0x5B15720AD2DAE0F1), UINT64_C(0x57B0E822E70EBAE8)},
            boost::int128::uint128_t {UINT64_C(0x47CF7C378CBD76E4), UINT64_C(0x6CF2CB52A1ED2703)},
            boost::int128::uint128_t {UINT64_C(0x70B241CBA4AEC69A), UINT64_C(0x8078A934376E0BEB)},
            boost::int128::uint128_t {UINT64_C(0x5804C611740EB0F7), UINT64_C(0x29D8647FD217E433)},
            boost::int128::uint128_t {UINT64_C(0x446DD22EC2D756F4), UINT64_C(0x6A9DF53E7B37DDE0)},
            boost::int128::uint128_t {UINT64_C(0x69EA8E8015960B9B), UINT64_C(0x6DD2FD9EADAD248E)},
            boost::int128::uint128_t {UINT64_C(0x519A7219E839199C), UINT64_C(0x484C63E6AACEB7E2)},
            boost::int128::uint128_t {UINT64_C(0x7D2FFBF669D57302), UINT64_C(0x6E7017CB0557D6B7)},
            boost::int128::uint128_t {UINT64_C(0x5F9B2E3602303132), UINT64_C(0x244256329B890B20)},
            boost::int128::uint128_t {UINT64_C(0x48B34F585ADC4696), UINT64_C(0x74A89A30256CB161)},
            boost::int128::uint128_t {UINT64_C(0x6E184681829F5C30), UINT64_C(0x38B07F2128EBDC8E)},
            boost::int128::uint128_t {UINT64_C(0x5302F09E1A0E6A06), UINT64_C(0xC4DB8F9E5351A46D)},
            boost::int128::uint128_t {UINT64_C(0x7CA95D38C5F30CE7), UINT64_C(0x9600E48EAB67CB06)},
            boost::int128::uint128_t {UINT64_C(0x5D38742FE19B9192), UINT64_C(0xCB125737D0803D42)},
            boost::int128::uint128_t {UINT64_C(0x456D3318EE6B6E67), UINT64_C(0x349951781A9F442F)},
            boost::int128::uint128_t {UINT64_C(0x66FF450B7954791C), UINT64_C(0x952B1858E092A9D4)},
            boost::int128::uint128_t {UINT64_C(0x4C18F01A8154DD7D), UINT64_C(0x9C23CA956DF963D4)},
            boost::int128::uint128_t {UINT64_C(0x7001625A34220388), UINT64_C(0x1C58631418460B5B)},
            boost::int128::uint128_t {UINT64_C(0x521BED95EF231EC3), UINT64_C(0xEAD8F6CB10A8C9CE)},
            boost::int128::uint128_t {UINT64_C(0x77EC77F075E7DA8F), UINT64_C(0x8B771A90066A8989)},
            boost::int128::uint128_t {UINT64_C(0x573E99E34E94DD70), UINT64_C(0x202AE3FA4F95F34F)},
            boost::int128::uint128_t {UINT64_C(0x7E76E3D54214A026), UINT64_C(0xA37BB520A720DEDF)},
            boost::int128::uint128_t {UINT64_C(0x5B5132003AC80EFA), UINT64_C(0xC49901F6019468D9)},
            boost::int128::uint128_t {UINT64_C(0x41B21C6FD88D9D1D), UINT64_C(0xB9FE98B518A760AB)},
            boost::int128::uint128_t {UINT64_C(0x5E2EA87F341F224D), UINT64_C(0x0516A2635FAD59A8)},
            boost::int128::uint128_t {UINT64_C(0x43446D8D3ADA13F9), UINT64_C(0x4592E0CA98CCB7C8)},
            boost::int128::uint128_t {UINT64_C(0x5FBECC6521B9FEBB), UINT64_C(0x9A772556B4111000)},
            boost::int128::uint128_t {UINT64_C(0x43E635127DF0F9BB), UINT64_C(0x8D0FCA69F82B9100)},
            boost::int128::uint128_t {UINT64_C(0x5FF784E5AA2CFDB1), UINT64_C(0x60DF4771B802E800)},
            boost::int128::uint128_t {UINT64_C(0x4395667F94DD90C9), UINT64_C(0x2D16B625E119E000)},
            boost::int128::uint128_t {UINT64_C(0x5EDD0EA6C0F80CBC), UINT64_C(0xDD5CBE74F2650000)},
            boost::int128::uint128_t {UINT64_C(0x4259EABB1D0B0F0E), UINT64_C(0x3748908ECACC0000)},
            boost::int128::uint128_t {UINT64_C(0x5C813FE5DCFBA9C8), UINT64_C(0x43D86B1C96400000)},
            boost::int128::uint128_t {UINT64_C(0x4044C9B616067C5D), UINT64_C(0xC4F54CCB40800000)},
            boost::int128::uint128_t {UINT64_C(0x5901EFDAC3A80B57), UINT64_C(0x0CD16F02A0000000)},
            boost::int128::uint128_t {UINT64_C(0x7ADD8DE861C08DFB), UINT64_C(0x3FFB8F91E0000000)},
            boost::int128::uint128_t {UINT64_C(0x5486A71B80C9E953), UINT64_C(0xA965A5C680000000)},
            boost::int128::uint128_t {UINT64_C(0x73ED21D64CD473E7), UINT64_C(0xF742E8E000000000)},
            boost::int128::uint128_t {UINT64_C(0x4F3E9B80A4325118), UINT64_C(0x7B12855000000000)},
            boost::int128::uint128_t {UINT64_C(0x6BFE65EDC2618290), UINT64_C(0x7306CE0000000000)},
            boost::int128::uint128_t {UINT64_C(0x4950C34F0AE28DFA), UINT64_C(0xE9EE4E0000000000)},
            boost::int128::uint128_t {UINT64_C(0x634AD9EEC02F4379), UINT64_C(0xC41CD00000000000)},
            boost::int128::uint128_t {UINT64_C(0x438B38735408E0E1), UINT64_C(0xD013C00000000000)},
            boost::int128::uint128_t {UINT64_C(0x5AE6CA44EB128FD7), UINT64_C(0xD70A000000000000)},
            boost::int128::uint128_t {UINT64_C(0x70C9E44CE0E77A5E), UINT64_C(0xF130000000000000)},
            boost::int128::uint128_t {UINT64_C(0x4B85BEA805F1B7F3), UINT64_C(0x73C0000000000000)},
            boost::int128::uint128_t {UINT64_C(0x4BF492B992EFF413), UINT64_C(0x4500000000000000)},
            boost::int128::uint128_t {UINT64_C(0x64D0BA7EA76484A9), UINT64_C(0x0000000000000000)}
        }
    };
}

constexpr auto erf_float_table_for(std::integral_constant<int, 3>) noexcept -> erf_float_table_t<3>
{
    return {
        {2, 3, 4, 5},
        {{22, 27, 30, 33, 36, 39, 41, 42}, {68, 74, 78, 82, 86, 90, 92, 92},
         {120, 126, 132, 138, 143, 148, 149, 149}, {174, 179, 183, 187, 191, 195, 198, 198}},
        {{12, 15, 17, 18, 20, 22, 23, 24}, {56, 60, 62, 64, 67, 69, 71, 73},
         {107, 111, 114, 117, 120, 122, 125, 128}, {163, 165, 168, 170, 172, 173, 175, 177}},
        {
            1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 3, 2,
            2, 3, 2, 3, 2, 2, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1,
            2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 2, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1,
            1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1,
            1, 2, 1, 0, 1, 0, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1,
            2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 0
        }
    };
}

template <int N>
struct erf_float_imp
{
    static constexpr erf_float_mant_t<N> mant {erf_float_mant_for(std::integral_constant<int, N> {})};
    static constexpr erf_float_table_t<N> t {erf_float_table_for(std::integral_constant<int, N> {})};
};

#if !(defined(__cpp_inline_variables) && __cpp_inline_variables >= 201606L) && (!defined(_MSC_VER) || _MSC_VER != 1900)

template <int N>
constexpr erf_float_mant_t<N> erf_float_imp<N>::mant;

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
template <int N, typename M>
constexpr auto erf_run(int k, const int lo, M s, const M& u, const M& w, const bool above) noexcept -> M
{
    const auto& t {erf_float_imp<N>::t};
    const auto& m {erf_float_imp<N>::mant.m};
    if (((k - lo) & 1) != 0)
    {
        --k;
        const auto p {static_cast<M>(erf_mul_high(u, s) >> t.shift[k])};
        s = above ? erf_as(m[k], s) - p : erf_as(m[k], s) + p;
    }
    while (k > lo)
    {
        k -= 2;
        const auto p1 {static_cast<M>(erf_mul_high(u, erf_as(m[k + 1], s)) >> t.shift[k])};
        const auto p2 {static_cast<M>(erf_mul_high(w, s) >> (t.shift[k] + t.shift[k + 1] + 1U))};
        s = above ? erf_as(m[k], s) - p1 + p2 : erf_as(m[k], s) + p1 + p2;
    }
    return s;
}

// g(x) for 1 <= x < 10 or 12 in N - 1 words by Horner's rule on the table of N words, each sum at the scale of its
// coefficient: g is below 4.25 units of 2^-62 or 8.34 units of 2^-126 off, from the steps, the tail and u
template <int N>
constexpr auto erf_g_float(const std::uint64_t whole, const fx<2>& frac, int& l) noexcept -> erf_mant_t<N>
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
    const auto s1 {erf_run<N>(last, t.split[i][j], erf_top(erf_float_imp<N>::mant.m[last]), u1, w1, above)};
    l = t.l[i];
    return erf_run<N>(t.split[i][j], iv.first, erf_widen<N>(s1), u2, w2, above);
}

// g(x) for x >= 10 or 12 in N - 1 words: 1 / x = r 2^-s with r = 1 / d for x = d 2^s, the series S in v = 1 / (2x^2) by
// Horner's rule from term 9 + 80 / whole or 9 + 384 / whole, one word above 6 + 120 / whole, and g = S r (2 / sqrt(pi))
// 2^-(s + 1), below 10 units of 2^-62 or 43.2 units of 2^-126 off
template <int N>
constexpr auto erf_g_asymptotic_float(const std::uint64_t whole, const fx<2>& frac, int& l) noexcept -> erf_mant_t<N>
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
    l = static_cast<int>(s) + 1;
    return erf_mant(g);
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
    const auto carry {fx_mul_add(fx_make_span<N>(sq), fx_make_span(frac.w), 2U * whole)};
    sq[N] = carry;
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
    int gl {};
    const auto gm {whole < 10U ? erf_g_float<2>(whole, frac, gl) : erf_g_asymptotic_float<2>(whole, frac, gl)};

    // erfc(x) = q 2^-(61 + l - a) 10^(-n - 1) with q = m g / 2^63, and d its digits times 10^19, or 10^20 below 10^18
    const auto q {((static_cast<sig_type>(r10.m.w[0]) * gm) >> 63U).low};
    const auto t {static_cast<sig_type>(q) * pow10(static_cast<sig_type>(19))};
    const auto sh {static_cast<unsigned>(61 + gl) - r10.a};
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
    int gl {};
    const auto gm {whole < 12U ? erf_g_float<3>(whole, frac, gl) : erf_g_asymptotic_float<3>(whole, frac, gl)};
    const auto pr {umul256(erf_mant(r10.m), gm)};
    const auto lz {static_cast<unsigned>(countl_zero(pr[3]))};
    const fx<2> q {{(pr[2] << lz) | (pr[1] >> (64U - lz)), (pr[3] << lz) | (pr[2] >> (64U - lz))}};
    const int e {2 - static_cast<int>(lz) + static_cast<int>(r10.a) - gl};
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
    int gl {};
    const auto m {erf_g_float<3>(whole, frac, gl) >> gl};
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
        return detail::check_non_finite(z);
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
