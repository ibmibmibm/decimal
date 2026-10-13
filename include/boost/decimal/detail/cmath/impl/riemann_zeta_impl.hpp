// Copyright 2024 Matt Borland
// Copyright 2024 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_RIEMANN_ZETA_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_RIEMANN_ZETA_IMPL_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/u256.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/impl/exp_impl.hpp>
#include <boost/decimal/detail/cmath/impl/gamma_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {

namespace riemann_zeta_detail {

template <bool b>
struct zeta_table_imp
{
    // |c_(k+1) / c_k| for c_k = B_2k / (2k)!, from k = 1
    static constexpr fx<3> ratio[46] = {
        {{UINT64_C(0x1111111111111111), UINT64_C(0x1111111111111111), UINT64_C(0x0111111111111111)}},
        {{UINT64_C(0x8618618618618618), UINT64_C(0x1861861861861861), UINT64_C(0x0186186186186186)}},
        {{UINT64_C(0x9999999999999999), UINT64_C(0x9999999999999999), UINT64_C(0x0199999999999999)}},
        {{UINT64_C(0x9DBCC48676F31219), UINT64_C(0x19DBCC48676F3121), UINT64_C(0x019DBCC48676F312)}},
        {{UINT64_C(0x9EB384D1E6B80519), UINT64_C(0x19EB384D1E6B8051), UINT64_C(0x019EB384D1E6B805)}},
        {{UINT64_C(0x22157CD9D73F5A06), UINT64_C(0x2FCAA6B7041343CC), UINT64_C(0x019EEF5419314605)}},
        {{UINT64_C(0x9EFE1405647A6BCA), UINT64_C(0xCAE0D231473897AD), UINT64_C(0x019EFE1405647A6B)}},
        {{UINT64_C(0x7E6A797D99A16F53), UINT64_C(0xA56E852C59E12280), UINT64_C(0x019F01BE4FD46F5D)}},
        {{UINT64_C(0x981E4737E1EB1C08), UINT64_C(0x0BF31E770D711D18), UINT64_C(0x019F02A8419F5AE2)}},
        {{UINT64_C(0x49472CDB05AB4502), UINT64_C(0x8B53A8E1DF3801E7), UINT64_C(0x019F02E2AC422E8B)}},
        {{UINT64_C(0xF7DF86BCCB2AD84B), UINT64_C(0x272F2B7377D56BF2), UINT64_C(0x019F02F144F10D42)}},
        {{UINT64_C(0x96E589882A1A84C2), UINT64_C(0xC193959465B4447A), UINT64_C(0x019F02F4EAE49CED)}},
        {{UINT64_C(0x09762273C84130DA), UINT64_C(0xD1B17E42DBB8F2BB), UINT64_C(0x019F02F5D45B443D)}},
        {{UINT64_C(0x213D4C0C1B26E643), UINT64_C(0x99AD24D8A117E895), UINT64_C(0x019F02F60EB83CB6)}},
        {{UINT64_C(0x99ED581FDC4ED017), UINT64_C(0xAEF0EA53D6B19F6B), UINT64_C(0x019F02F61D4F6720)}},
        {{UINT64_C(0x8C38B0D47BB9F0CC), UINT64_C(0xCA3334E8F685E49B), UINT64_C(0x019F02F620F52F8A)}},
        {{UINT64_C(0x2332626AB869E4D6), UINT64_C(0x0D0484A35E1DECF2), UINT64_C(0x019F02F621DEA167)}},
        {{UINT64_C(0x155D14AC73F8B247), UINT64_C(0x32A73C2BF41EE95A), UINT64_C(0x019F02F62218FDD7)}},
        {{UINT64_C(0x887E4BD6089AC8DB), UINT64_C(0x7747840D7691E232), UINT64_C(0x019F02F6222794F2)}},
        {{UINT64_C(0x8784D6468B801A78), UINT64_C(0x3292400A7551A1BA), UINT64_C(0x019F02F6222B3AB9)}},
        {{UINT64_C(0x052083525F7F9DF4), UINT64_C(0xDEF70298781DA512), UINT64_C(0x019F02F6222C242A)}},
        {{UINT64_C(0xF16A26AFF4121478), UINT64_C(0x49CB190128A3CD9F), UINT64_C(0x019F02F6222C5E87)}},
        {{UINT64_C(0xDA4779C287F772C7), UINT64_C(0x647871073D63307E), UINT64_C(0x019F02F6222C6D1E)}},
        {{UINT64_C(0xC30C1EFA8FDDC772), UINT64_C(0x2B22ECA30382F758), UINT64_C(0x019F02F6222C70C4)}},
        {{UINT64_C(0x9AB6D5E04AC72F77), UINT64_C(0x9CCD7345C408C1E2), UINT64_C(0x019F02F6222C71AD)}},
        {{UINT64_C(0xF70FA198485B5861), UINT64_C(0xF938123C35DF4645), UINT64_C(0x019F02F6222C71E7)}},
        {{UINT64_C(0x77618D11AC9AAEAE), UINT64_C(0x9052B9AD20BF055F), UINT64_C(0x019F02F6222C71F6)}},
        {{UINT64_C(0xFC0F135096F87C38), UINT64_C(0x36196380D5F49157), UINT64_C(0x019F02F6222C71FA)}},
        {{UINT64_C(0x0798539B24BA21EA), UINT64_C(0x1F8B0DF4D0DE22CD), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xEBB9BC97DF17BA7E), UINT64_C(0x59E77891B4A9E1B3), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xD2AD61BDB06164C3), UINT64_C(0x687E9338EA9EBF0B), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x317D56E73A431E6A), UINT64_C(0x6C2459E2B7C6D7E5), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x3FC138F7DE56956D), UINT64_C(0x6D0DCB8D2B0768F1), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x908943C6084EE4B3), UINT64_C(0x6D4827F7C7D6802F), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x10DE80F322743368), UINT64_C(0x6D56BF126F0A281B), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xB5DB532854202F82), UINT64_C(0x6D5A64D918D70EC3), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x66C27A79D6894A3E), UINT64_C(0x6D5B4E4AC34A480F), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xD3D609FCE0A7C6A8), UINT64_C(0x6D5B88A72DE71657), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xC488759C603D73E3), UINT64_C(0x6D5B973E488E49E8), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xDF88587D0A729D7E), UINT64_C(0x6D5B9AE40F3816CC), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xE298AE89BDEBDFEF), UINT64_C(0x6D5B9BCD80E28A05), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x22F3EADD6BD29A49), UINT64_C(0x6D5B9C07DD4D26D4), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xB2FF1397BAF88D10), UINT64_C(0x6D5B9C167467CE07), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x97009266E80E1BF9), UINT64_C(0x6D5B9C1A1A2E77D4), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0xD000CD48FD3F8264), UINT64_C(0x6D5B9C1B03A02247), UINT64_C(0x019F02F6222C71FB)}},
        {{UINT64_C(0x9E40D7EA356D31AB), UINT64_C(0x6D5B9C1B3DFC8CE4), UINT64_C(0x019F02F6222C71FB)}}
    };

    // log2(p) for p = 3, 5, 7, 11 and 13
    static constexpr fx<3> log2_prime[5] = {
        {{UINT64_C(0x493CF9A8E8966C10), UINT64_C(0xE802C48281A2EB74), UINT64_C(0x6570068E7EF5A1E7)}},
        {{UINT64_C(0xD96C55FE37B3AD4E), UINT64_C(0x492BF6FF4DAFDB4C), UINT64_C(0x949A784BCD1B8AFE)}},
        {{UINT64_C(0xCE31A952005C59F2), UINT64_C(0xCAB240E904F96A11), UINT64_C(0xB3ABB3FAA02166CC)}},
        {{UINT64_C(0xCF9CC2EDD043A256), UINT64_C(0x3EBE1990555535AF), UINT64_C(0xDD6753E032EA0EFE)}},
        {{UINT64_C(0x873D20F4A4019DB1), UINT64_C(0x5DB68721CA60D447), UINT64_C(0xECD4011C8F11979A)}}
    };
};

template <bool b>
constexpr fx<3> zeta_table_imp<b>::ratio[46];

template <bool b>
constexpr fx<3> zeta_table_imp<b>::log2_prime[5];

using zeta_table = zeta_table_imp<true>;

using gamma_detail::wide;

// 2^-t for t >= 0: 2^-t = 2^(-whole - 1) * 2^g with g = 1 - f, and 2^g = 10^(g log10(2)) = 2^a * m, where one word
// takes the plain series of exp
template <int N, bool Wide = false>
constexpr auto zeta_pow2(const wide<N>& t) noexcept -> fx<N>
{
    const std::uint64_t whole {gamma_detail::gamma_whole(t)};
    if (whole >= static_cast<std::uint64_t>(64 * N))
    {
        return fx<N> {};
    }
    const auto m {exp_detail::exp10_frac<N, (Wide && N > 1)>(fx_mul(fx_sub(fx_one<N>(), fx_frac<N>(t.w)), fx_top<N>(fx_table::log10_2)))};
    return fx_shift(m.m, static_cast<int>(m.a) - static_cast<int>(whole) - 1);
}

// 1 / d for d >= 1/2, with one more Newton step for Wide
template <int N, bool Wide = false>
constexpr auto zeta_inv(const wide<N>& d) noexcept -> fx<N>
{
    const auto f {gamma_detail::gamma_whole(d) == 0U ? gamma_detail::gamma_float<N> {fx_shl(fx_low<N>(d.w), 1U), -1, 0} : gamma_detail::gamma_normalize(d)};
    auto r {fx_recip(f.v)};
    BOOST_DECIMAL_IF_CONSTEXPR (Wide)
    {
        r = fx_mul(r, fx_sub(fx_int<N>(2U), fx_mul(f.v, r)));
    }
    return fx_shift(r, -f.n2);
}

// x - 1 = d with d < 0 in neg: |d| < 1/2 in small, else 1 / |d| in inv
template <int N>
struct zeta_delta
{
    fx<N> d;
    fx<N> inv;
    bool small;
    bool neg;
};

// S(x) / 16 for S(x) = zeta(x) - 1 / (x - 1) and 0 < x < 256, by Euler-Maclaurin with m = 8 or 16: S = sum n^-x for n < m +
// (m^-d - 1) / d + m^-x / 2 + sum T_k with T_k = c_k (x)_(2k-1) m^(-x-2k+1), until T_k is below the fixed point
template <int N, bool Wide = false>
constexpr auto zeta_em(const wide<N>& x, const zeta_delta<N>& dl) noexcept -> fx<N>
{
    constexpr int m {N == 1 ? 8 : 16};
    constexpr unsigned lm {N == 1 ? 3U : 4U};
    constexpr int primes[5] {3, 5, 7, 11, 13};
    const auto one {fx_one<N>()};

    // n^-x from 2^-(x log2(p)) for the primes p, then their products
    fx<N> pw[16] {};
    pw[1] = one;
    pw[2] = zeta_pow2<N, Wide>(x);
    for (int i {}; i < 5 && primes[i] < m; ++i)
    {
        pw[primes[i]] = zeta_pow2<N, Wide>(fx_wide_mul(x, fx_top<N>(zeta_table::log2_prime[i])));
    }
    wide<N> total {};
    fx<N> sum0 {};
    for (int n {1}; n < m; ++n)
    {
        // n = 1 and the primes 2, 3, 5, 7, 11 and 13 (the bits of 0x28AC) are set above
        if (n >= 4 && ((0x28ACU >> n) & 1U) == 0U)
        {
            pw[n] = n % 2 == 0 ? fx_mul(pw[2], pw[n / 2]) : fx_mul(pw[3], pw[n / 3]);
        }
        BOOST_DECIMAL_IF_CONSTEXPR (Wide)
        {
            total = fx_wide_add(total, fx_widen(pw[n]));
        }
        else
        {
            sum0 = fx_add(sum0, fx_shr(pw[n], 4U));
        }
    }
    const auto mx {N == 1 ? fx_mul(pw[2], pw[4]) : Wide ? zeta_pow2<N, Wide>(gamma_detail::gamma_shl(x, 2U)) : fx_mul(pw[4], pw[4])};
    const auto sum {Wide ? gamma_detail::gamma_shr(fx_wide_add(total, fx_widen(fx_shr(mx, 1U))), 4U) : fx_add(sum0, fx_shr(mx, 5U))};

    // (m^-d - 1) / d < 0: with h = d log2(m) / 2, -ln(m) 2^-h sinh(t) / t for t = h ln(2), else |m mx - 1| / |d|
    fx<N> drop {};
    if (dl.small)
    {
        const auto h {N == 1 ? fx_mul_shr(dl.d, 3U, 1U) : fx_shl(dl.d, 1U)};
        const auto p2 {dl.neg ? fx_shl(zeta_pow2<N, Wide>(fx_widen(fx_sub(one, h))), 1U) : zeta_pow2<N, Wide>(fx_widen(h))};
        const auto t {fx_mul(h, fx_top<N>(fx_table::ln2))};
        const auto tt {fx_mul(t, t)};
        auto shc {one};
        // The terms of sinh(t) / t to below the last bit
        for (std::uint64_t k {N == 1 ? 8U : N == 2 ? (Wide ? 17U : 14U) : 22U}; k > 0U; --k)
        {
            shc = fx_add(one, fx_div_small(fx_mul(tt, shc), 2U * k * (2U * k + 1U)));
        }
        drop = fx_mul(fx_mul(fx_mul_shr(fx_top<N>(fx_table::ln2), lm, 4U), p2), shc);
    }
    else
    {
        const auto inv_m {fx_shr(one, lm)};
        const auto diff {fx_less(mx, inv_m) ? fx_sub(inv_m, mx) : fx_sub(mx, inv_m)};
        drop = fx_mul(N == 1 ? fx_shr(diff, 1U) : diff, dl.inv);
    }

    // T_1 = x m^-x / (12 m), and T_(k+1) = -T_k |r_k| (x + 2k - 1) (x + 2k) / m^2, with x / 256 and all T_k times 2^g
    // for T_1 2^g >= 1/4, to keep the truncation of each step below the fixed point
    const auto xq {gamma_detail::gamma_shr(x, 8U)};
    auto term {fx_div_small(fx_mul(xq, mx), 3U * static_cast<unsigned>(m) / 4U)};
    const int g {Wide ? fx_lead(term) - 3 : 0};
    term = fx_shift(term, g);
    auto pos {term};
    fx<N> neg {};
    constexpr int ratios {static_cast<int>(sizeof(zeta_table::ratio) / sizeof(zeta_table::ratio[0]))};
    for (int k {1}; k < ratios && fx_lead(term) < 64 * N - (Wide ? 0 : 4) - g; ++k)
    {
        fx<N> a {};
        a.w[N - 1] = static_cast<std::uint64_t>(2 * k - 1) << 54U;
        fx<N> b {};
        b.w[N - 1] = static_cast<std::uint64_t>(2 * k) << 54U;
        term = fx_mul(fx_mul(fx_mul(term, fx_add(xq, a)), fx_add(xq, b)), fx_top<N>(zeta_table::ratio[k - 1]));
        term = fx_shl(term, 16U - 2U * lm);
        if (k % 2 == 1)
        {
            neg = fx_add(neg, term);
        }
        else
        {
            pos = fx_add(pos, term);
        }
    }
    return fx_sub(fx_add(sum, fx_shift(fx_sub(pos, neg), -g)), drop);
}

// The digits of a result, times 10^exp10 and with the sign neg
struct zeta_digits
{
    boost::int128::uint128_t digits;
    int exp10;
    bool neg;
};

// zeta(x) for x = sig * 10^e with 10^-30 <= x <= 128 and x != 1 in N words, x within 1/2 of 1 in small and below 1 in neg:
// zeta = 1 / (x - 1) + S(x), where near the pole 1 / (x - 1) is an exact division and keeps its fraction
template <int N, bool Wide = false>
constexpr auto zeta_positive_digits(const boost::int128::uint128_t& sig, const int e, const bool small, const bool neg) noexcept -> zeta_digits
{
    using sig_type = boost::int128::uint128_t;

    const auto xw {gamma_detail::gamma_wide<N>(sig, e)};
    zeta_delta<N> dl {{}, {}, small, neg};
    sig_type dsig {};
    if (dl.small)
    {
        const auto one {pow10(static_cast<sig_type>(-e))};
        dsig = dl.neg ? one - sig : sig - one;
        dl.d = fx_from_decimal<N>(dsig, -e);
    }
    else
    {
        const auto one {fx_wide_int<N>(1U)};
        dl.inv = zeta_inv<N, Wide>(dl.neg ? fx_wide_add(one, xw, true) : fx_wide_add(xw, one, true));
    }
    const auto s {fx_shl(zeta_em<N, Wide>(xw, dl), 4U)};

    if (!dl.small)
    {
        // zeta > 1 for x > 1, and zeta < -1/2 for 0 < x < 1
        const auto digits {fx_scale(pow10(static_cast<sig_type>(37)), dl.neg ? fx_sub(dl.inv, s) : fx_add(dl.inv, s))};
        const auto bound {dl.neg ? pow10(static_cast<sig_type>(36)) * 5U : pow10(static_cast<sig_type>(37))};
        return {digits < bound ? bound : digits, -37, dl.neg};
    }

    // 10^(37 + nd) / dsig = q + rem / dsig and S * 10^(37 + nd + e), with the fractions of both
    const int nd {num_digits(dsig)};
    const auto qr {impl::div_mod(umul256(pow10(static_cast<sig_type>(37)), pow10(static_cast<sig_type>(nd))), dsig)};
    const auto q {static_cast<sig_type>(qr.quotient)};
    const auto rem {static_cast<sig_type>(qr.remainder)};
    const auto fq {rem == 0U ? UINT64_C(0) : static_cast<sig_type>(impl::div_mod(umul256(rem, sig_type {1U, 0U}), dsig).quotient).low};
    std::uint64_t fs {};
    const auto a {fx_scale_frac(pow10(static_cast<sig_type>(37 + nd + e)), s, fs)};
    return {dl.neg ? q - a - (fq < fs ? 1U : 0U) : q + a + (fq + fs < fq ? 1U : 0U), -37 - nd - e, dl.neg};
}

// zeta(x) for 10^-30 <= x <= 128, x != 1, rounded once in the current mode
template <typename T>
constexpr auto zeta_positive(const T x) noexcept -> T
{
    constexpr int N {fx_words(precision_v<T>)};

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool small {x > T {5, -1} && x < T {15, -1}};
    const bool neg {x < T {1}};
    BOOST_DECIMAL_IF_CONSTEXPR (N == 2)
    {
        // For 1/2 < x < 54, a first try in one word: its error is below 300 units of 2^-62 within 1/2 of the pole,
        // 250 for x < 4, else 100. From 54, zeta - 1 < 10^-16 is next to the rounding boundary at 1.
        if (small || (!neg && x < T {54}))
        {
            const auto d1 {zeta_positive_digits<1, true>(sig, e, small, neg)};
            T r {};
            if (fx_try_round_digits(d1.digits, d1.exp10, d1.neg, small ? 300U : x < T {4} ? 250U : 100U, r))
            {
                return r;
            }
        }
    }
    BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
    {
        // For x > 1/2, a first try in two words: its error is below 116 units of the digits within 1/2 of the
        // pole, else 44
        if (small || !neg)
        {
            const auto d2 {zeta_positive_digits<2, true>(sig, e, small, neg)};
            if (!fx_near_boundary(d2.digits, small ? 116U : 44U))
            {
                return fx_round<T>(d2.digits, d2.exp10, d2.neg);
            }
        }
    }
    const auto d {zeta_positive_digits<N>(sig, e, small, neg)};
    return fx_round<T>(d.digits, d.exp10, d.neg);
}

// zeta(x) for x = -a with a = sig * 10^e >= 10^-30, rounded once in the current mode: zeta(x) = 2^x pi^(x - 1)
// sin(pi x / 2) gamma(1 - x) zeta(1 - x) in logarithms, with sin from the distance of a to the next even integer
template <typename T>
constexpr auto zeta_negative(const boost::int128::uint128_t& sig, const int e) noexcept -> T
{
    constexpr int p {precision_v<T>};
    constexpr int N {fx_words(p)};
    using sig_type = boost::int128::uint128_t;

    // a = whole + f with f = fs * 10^e, and whole mod 4
    sig_type fs {sig};
    unsigned quad {};
    if (e >= 0)
    {
        fs = 0U;
        // 10 = 2 and 100 = 0 (mod 4)
        quad = e >= 2 ? 0U : static_cast<unsigned>(((e == 1 ? sig * 2U : sig) % 4U).low);
    }
    else if (e >= -38)
    {
        const auto unit {pow10(static_cast<sig_type>(-e))};
        quad = static_cast<unsigned>(((sig / unit) % 4U).low);
        fs = sig % unit;
    }
    const bool odd {(quad & 1U) != 0U};
    const bool neg {quad < 2U};
    if (fs == 0U && !odd)
    {
        return T {0};
    }
    // |zeta(-a)| overflows for a >= 10^4
    if (e + num_digits(sig) > 4)
    {
        return exp_detail::exp10_beyond<T, N>(exp_detail::exp_kind::overflow, neg);
    }

    // |sin(pi x / 2)| = sin(pi u / 2) for u = f, or 1 - f for an odd whole, and ln|sin| - ln(pi / 2) = ln(u) +
    // ln(sin(pi v) / (pi v)) with v = u / 2 = 5 us * 10^(ue - 1)
    const auto us {odd ? (e < 0 ? pow10(static_cast<sig_type>(-e)) - fs : sig_type {1U}) : fs};
    const int ue {odd && e >= 0 ? 0 : e};
    const auto sinc {gamma_detail::gamma_sinc(us * 5U, 1 - ue, std::integral_constant<int, N> {})};
    const auto lsin {fx_wide_add(log_detail::log_wide_of<N>(us, num_digits(us), ue, log_detail::log_base::e),
                                 log_detail::log_wide_of(fx_shl(fx_clamp_below_one(sinc), 1U), -1, 0, log_detail::log_base::e))};

    // ln gamma(1 + a) = st - ln(c)
    const auto aw {gamma_detail::gamma_wide<N>(sig, e)};
    const auto a1 {fx_wide_add(aw, fx_wide_int<N>(1U))};
    const auto g {gamma_detail::gamma_shift(a1)};

    // ln zeta(1 + a) = ln(1 + a S) - ln(a) for a < 1/2, else ln(1 / a + S), and 0 where 2^-a is below the fixed point
    wide<N> lz {};
    const auto whole {gamma_detail::gamma_whole(aw)};
    if (whole < static_cast<std::uint64_t>(64 * N))
    {
        fx<N> half {};
        half.w[N - 1] = UINT64_C(1) << 61U;
        const bool below_one {whole == 0U};
        const auto af {below_one ? fx_from_decimal<N>(sig, -e) : fx<N> {}};
        const bool small {below_one && fx_less(af, half)};
        zeta_delta<N> dl {af, {}, small, false};
        if (!small)
        {
            dl.inv = zeta_inv(aw);
        }
        const auto s {fx_shl(zeta_em(a1, dl), 4U)};
        if (small)
        {
            lz = fx_wide_add(log_detail::log_wide_of(fx_add(fx_one<N>(), fx_mul(af, s)), 0, 0, log_detail::log_base::e),
                             log_detail::log_wide_of<N>(sig, num_digits(sig), e, log_detail::log_base::e), true);
        }
        else
        {
            const auto z {fx_add(dl.inv, s)};
            const bool two {z.w[N - 1] >= (UINT64_C(1) << 63U)};
            lz = log_detail::log_wide_of(two ? fx_shr(z, 1U) : z, two ? 1 : 0, 0, log_detail::log_base::e);
        }
    }

    // ln|zeta| = -a ln(2 pi) - ln(pi) + ln|sin| + ln gamma(1 + a) + ln zeta(1 + a), with ln(pi / 2) - ln(pi) = -ln(2)
    // and ln(2 pi) = 2 ln(sqrt(2 pi))
    const auto l2pi {fx_shl(fx_top<N>(fx_table::ln_sqrt_2pi), 1U)};
    auto st {fx_wide_add(fx_wide_add(lsin, g.st), lz)};
    st = fx_wide_add(st, fx_wide_add(fx_wide_mul(aw, l2pi), fx_widen(fx_top<N>(fx_table::ln2))), true);
    return gamma_detail::gamma_exp_round<T>(gamma_detail::gamma_exp_digits<T>(st, gamma_detail::gamma_frecip(g.c)), neg);
}

} // namespace riemann_zeta_detail

} //namespace detail
} //namespace decimal
} //namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_RIEMANN_ZETA_IMPL_HPP
