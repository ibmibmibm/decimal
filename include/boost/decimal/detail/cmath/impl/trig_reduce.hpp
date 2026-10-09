// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_TRIG_REDUCE_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_TRIG_REDUCE_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>
#include <boost/decimal/detail/u256.hpp>
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
namespace trig {

// Words of 2/pi in the product window, words of the fraction that go into r, fixed-point words, and bin_frac
// 64-bit fraction words for |x| < 10^19, where the fraction has at most 27, 58 and 117 leading zero bits.
template <typename T>
struct trig_traits;

// The window leaves window - 1 - frac = 2, 2 and 7 words for the zero words that lead the fraction, whose worst cases
// have 9, 19 and 37 zero digits. The frac words after them keep at least 9 frac - 8 = 19, 46 and 55 digits.
template <>
struct trig_traits<decimal32_t> { static constexpr int window = 6; static constexpr int frac = 3; static constexpr int words = 1; static constexpr int bin_frac = 2; };

template <>
struct trig_traits<decimal_fast32_t> : trig_traits<decimal32_t> {};

template <>
struct trig_traits<decimal64_t> { static constexpr int window = 9; static constexpr int frac = 6; static constexpr int words = 2; static constexpr int bin_frac = 3; };

template <>
struct trig_traits<decimal_fast64_t> : trig_traits<decimal64_t> {};

template <>
struct trig_traits<decimal128_t> { static constexpr int window = 15; static constexpr int frac = 7; static constexpr int words = 3; static constexpr int bin_frac = 5; };

template <>
struct trig_traits<decimal_fast128_t> : trig_traits<decimal128_t> {};

constexpr std::uint64_t word_base {UINT64_C(1000000000)};

// The most words of 9 digits in m * 10^shift, which has at most 38 + 8 digits
constexpr std::size_t m_words_max {6};

// The words of 2/pi in the table
constexpr int two_over_pi_words {static_cast<int>(sizeof(fx_table::two_over_pi) / sizeof(fx_table::two_over_pi[0]))};

// The digits of sig in trig_arg: one word keeps 19, which fit in 64 bits.
constexpr auto sig_digits(const int words) noexcept -> int
{
    return fx_scale10(words) + 1;
}

// x = n*pi/2 + r (mod 2pi), with |r| <= pi/4.
template <int N>
struct trig_arg
{
    boost::int128::uint128_t sig; // |r| = (sig + tail / 2^64) * 10^-k, with sig of sig_digits(N) digits
    std::uint64_t tail;
    fx<N> rf;                     // |r| in fixed point, if fixed is set
    int k;
    unsigned n;                   // the quadrant 0 to 3, which includes the sign of x
    bool fixed;
    bool neg;                     // the sign of r, which includes the sign of x
    bool zero;                    // x = 0: sig and k are not set
};

// |r| in fixed point: the binary reduction gives it, else it comes from sig and k.
template <int N>
constexpr auto fixed_r(const trig_arg<N>& r) noexcept -> fx<N>
{
    return r.fixed ? r.rf : fx_from_decimal<N>(r.sig, r.k);
}

// (sig + tail / 2^64) * 10^-k in N words, with sig and k first scaled by 10^up, where tail * 10^-k * 2^(64N - 2) / 2^64
// is below one word
template <int N, int W>
constexpr auto decimal_r_with_tail(const trig_arg<W>& r, const int up) noexcept -> fx<N>
{
    const auto scale {pow10(static_cast<boost::int128::uint128_t>(up))};
    fx<N> t {};
    boost::int128::detail::umul(r.tail, fx_from<N>(scale, r.k + up).w[0], t.w[0]);
    return fx_add(fx_from<N>(r.sig * scale, r.k + up), t);
}

// The same with the tail for two or three words, or the binary reduction's |r|
template <int N>
constexpr auto fixed_r_with_tail(const trig_arg<N>& r) noexcept -> fx<N>
{
    return r.fixed ? r.rf : decimal_r_with_tail<N>(r, 0);
}

// fx_round for a sig of sig_digits digits
template <typename T>
constexpr auto trig_round(const boost::int128::uint128_t& sig, const int e, const bool neg) noexcept -> T
{
    using digits_type = fx_digits_t<trig_traits<T>::words>;
    return fx_round<T>(static_cast<digits_type>(sig), e, neg);
}

// fx_try_round for a v of 18 to 20 digits with an error below w units
template <typename T>
constexpr auto trig_try_round(const std::uint64_t v, const int e, const bool neg, const std::uint64_t w, T& out) noexcept -> bool
{
    constexpr auto d19 {pow10(static_cast<std::uint64_t>(19))};
    if (v < d19 / 10U)
    {
        return fx_try_round(v * 10U, e - 1, neg, 10U * w, out);
    }
    if (v >= d19)
    {
        // Dropping a digit adds below one unit, and w / 10 rounds down.
        return fx_try_round(v / 10U, e + 1, neg, w / 10U + 2U, out);
    }
    return fx_try_round(v, e, neg, w, out);
}

// The significand in words of 9 digits, low word first; returns the count of words.
constexpr auto to_words(std::uint64_t m, const fx_span<std::uint64_t, m_words_max> out) noexcept -> int
{
    int n {};
    while (m != 0U)
    {
        out[n++] = m % word_base;
        m /= word_base;
    }
    return n;
}

constexpr auto to_words(const boost::int128::uint128_t& m, const fx_span<std::uint64_t, m_words_max> out) noexcept -> int
{
    constexpr std::uint64_t base2 {UINT64_C(1000000000000000000)};
    const auto hi {static_cast<std::uint64_t>(m / base2)};
    const auto lo {static_cast<std::uint64_t>(m % base2)};
    out[0] = lo % word_base;
    out[1] = lo / word_base;
    out[2] = hi % word_base;
    out[3] = hi / word_base;
    // Zero high words do not change the product.
    return 4;
}

// m * 10^shift in words of 9 digits, low word first, for 0 <= shift < 9; returns the words.
template <typename Significand>
constexpr auto scaled_words(const Significand m, const int shift, const fx_span<std::uint64_t, m_words_max> out) noexcept -> fx_span<std::uint64_t>
{
    int count {to_words(m, out)};
    const auto scale {pow10(static_cast<std::uint64_t>(shift))};
    std::uint64_t carry {};
    for (int i {}; i < count; ++i)
    {
        const std::uint64_t v {out[i] * scale + carry};
        out[i] = v % word_base;
        carry = v / word_base;
    }
    if (carry != 0U)
    {
        out[count++] = carry;
    }
    return out.first(count);
}

// The low Window words of m * 10^(9 * word_exp) * 2/pi, low word first. Word Window - 1 is the unit
// word, and 10^9 = 0 (mod 4), thus its low two bits are the quadrant and the higher words are not needed.
template <std::size_t Window>
constexpr auto times_two_over_pi(const fx_span<const std::uint64_t> m, const int word_exp, const fx_span<std::uint64_t, Window> out) noexcept -> void
{
    constexpr int window {static_cast<int>(Window)};
    const int count {static_cast<int>(m.size())};
    // two_over_pi[i] has the weight 10^(-9(i + 1)), thus word u of the window is table word Window - 2 - u + word_exp.
    std::uint64_t w[Window] {};
    for (int u {}; u < window; ++u)
    {
        const int i {window - 2 - u + word_exp};
        w[u] = (i >= 0 && i < two_over_pi_words) ? static_cast<std::uint64_t>(fx_table::two_over_pi[i]) : 0U;
    }

    std::uint64_t carry {};
    for (int col {}; col < window; ++col)
    {
        std::uint64_t acc {carry};
        for (int i {}; i < count && i <= col; ++i)
        {
            acc += m[i] * w[col - i];
        }
        out[col] = acc % word_base;
        carry = acc / word_base;
    }
}

// The fraction f in [0, 1) of Words words, high word first, to 1 - f when f >= 1/2; returns whether it did.
template <std::size_t Words>
constexpr auto fold_half(const fx_span<std::uint64_t, Words> f) noexcept -> bool
{
    constexpr int words {static_cast<int>(Words)};
    if (f[0] < word_base / 2U)
    {
        return false;
    }
    for (int i {}; i < words; ++i)
    {
        f[i] = word_base - 1U - f[i];
    }
    for (int i {words - 1}; i >= 0; --i)
    {
        if (++f[i] < word_base)
        {
            break;
        }
        f[i] = 0U;
    }
    return true;
}

// The Frac words of the fraction from word lead on, high word first, times pi/2: 2 * Frac + 1 words,
// low word first.
template <int Frac, std::size_t Words>
constexpr auto times_half_pi(const fx_span<std::uint64_t, Words> f, const int lead, const fx_span<std::uint64_t, static_cast<std::size_t>(2 * Frac + 1)> out) noexcept -> void
{
    std::uint64_t a[static_cast<std::size_t>(Frac)] {};
    std::uint64_t b[static_cast<std::size_t>(Frac + 1)] {};
    for (int u {}; u < Frac; ++u)
    {
        const int i {lead + Frac - 1 - u};
        a[u] = i < static_cast<int>(Words) ? f[i] : 0U;
    }
    for (int v {}; v <= Frac; ++v)
    {
        b[v] = static_cast<std::uint64_t>(fx_table::pio2[Frac - v]);
    }

    std::uint64_t carry {};
    for (int col {}; col <= 2 * Frac; ++col)
    {
        std::uint64_t acc {carry};
        for (int u {col - Frac < 0 ? 0 : col - Frac}; u < Frac && u <= col; ++u)
        {
            acc += a[u] * b[col - u];
        }
        out[col] = acc % word_base;
        carry = acc / word_base;
    }
}

// The first Digits digits of the nonzero words r, low word first, as sig * 10^exp10 with word 0 at 10^0, and
// with Tail the digits after them as a fraction of 64 bits in tail.
template <int Digits, bool Tail>
constexpr auto first_digits(const fx_span<const std::uint64_t> r, int& exp10, std::uint64_t& tail) noexcept -> boost::int128::uint128_t
{
    int top {static_cast<int>(r.size()) - 1};
    while (r[top] == 0U)
    {
        --top;
    }
    boost::int128::uint128_t sig {r[top]};
    int digits {num_digits(r[top])};
    exp10 = 9 * top;
    for (int i {top - 1}; i >= 0; --i)
    {
        if (digits + 9 <= Digits)
        {
            sig = sig * word_base + r[i];
            digits += 9;
            exp10 -= 9;
        }
        else
        {
            const int take {Digits - digits};
            const auto low {pow10(static_cast<std::uint64_t>(9 - take))};
            sig = sig * (word_base / low) + r[i] / low;
            exp10 -= take;

            // The next 18 digits d, and tail = d * 2^64 / 10^18 = 18d + d * 0.446744073709551616
            BOOST_DECIMAL_IF_CONSTEXPR (Tail)
            {
                const std::uint64_t below {(i >= 1 ? r[i - 1] : 0U) * word_base + (i >= 2 ? r[i - 2] : 0U)};
                const std::uint64_t d {(r[i] % low) * pow10(static_cast<std::uint64_t>(9 + take)) + below / low};
                std::uint64_t hi {};
                boost::int128::detail::umul(d, UINT64_C(0x725DD1D243ABA0E7), hi);
                tail = 18U * d + hi;
            }
            break;
        }
    }
    return sig;
}

// -x = (-n)*pi/2 + (-r)
BOOST_DECIMAL_FORCE_INLINE constexpr auto negate_arg(const bool xneg, unsigned& n, bool& rneg) noexcept -> void
{
    if (xneg)
    {
        n = (4U - n) & 3U;
        rneg = !rneg;
    }
}

// A sig below 10^(digits - 1) takes one more digit from the tail, and k goes up by one.
template <int N>
BOOST_DECIMAL_FORCE_INLINE constexpr auto fill_digit(trig_arg<N>& out, const int digits, int& k) noexcept -> void
{
    if (out.sig < pow10(static_cast<boost::int128::uint128_t>(digits - 1)))
    {
        const auto t {static_cast<boost::int128::uint128_t>(out.tail) * 10U};
        out.sig = out.sig * 10U + t.high;
        out.tail = t.low;
        ++k;
    }
}

// Payne-Hanek reduction in words of 9 decimal digits for x = m * 10^e: multiply m by the window of
// 2/pi which starts at the exponent, then the integer part mod 4 is the quadrant.
template <typename T, typename Significand>
constexpr auto trig_reduce(const Significand m, const int e, const bool xneg) noexcept -> trig_arg<trig_traits<T>::words>
{
    constexpr int words {trig_traits<T>::words};
    constexpr int window {trig_traits<T>::window};
    constexpr int frac {trig_traits<T>::frac};
    static_assert(window - 2 + std::numeric_limits<T>::max_exponent10 / 9 < two_over_pi_words, "2/pi must cover the largest exponent");

    // m * 10^e = (m * 10^shift) * 10^(9 * word_exp)
    const int shift {((e % 9) + 9) % 9};
    const int word_exp {(e - shift) / 9};
    std::uint64_t m_words[m_words_max] {};
    std::uint64_t prod[static_cast<std::size_t>(window)] {};
    times_two_over_pi(scaled_words(m, shift, fx_make_span(m_words)), word_exp, fx_make_span(prod));

    // The fraction, high word first. Above one half, go to the next quadrant with r < 0.
    std::uint64_t f[static_cast<std::size_t>(window - 1)] {};
    for (int i {}; i < window - 1; ++i)
    {
        f[i] = prod[window - 2 - i];
    }
    unsigned n {static_cast<unsigned>(prod[window - 1] & 3U)};
    bool rneg {fold_half(fx_make_span(f))};
    if (rneg)
    {
        n = (n + 1U) & 3U;
    }

    negate_arg(xneg, n, rneg);

    trig_arg<words> out {};
    out.n = n;
    out.neg = rneg;

    // The fraction is not zero: its worst case has 9, 19 and 37 leading zero digits (see trig_traits).
    int lead {};
    while (lead < window - 1 && f[lead] == 0U)
    {
        ++lead;
    }

    std::uint64_t r[static_cast<std::size_t>(2 * frac + 1)] {};
    times_half_pi<frac>(fx_make_span(f), lead, fx_make_span(r));

    // Word 0 of r has the weight 10^-9(lead + 2 * frac).
    int exp10 {};
    out.sig = first_digits<sig_digits(words), words == 3>(fx_make_span(r), exp10, out.tail);
    out.k = 9 * (lead + 2 * frac) - exp10;
    return out;
}

// The first count words of a times v, low word first, in place; returns the new count.
template <std::size_t M>
constexpr auto times_word(const fx_span<std::uint64_t, M> a, int count, const std::uint64_t v) noexcept -> int
{
    std::uint64_t carry {};
    for (int i {}; i < count; ++i)
    {
        a[i] = mul_add(a[i], v, 0U, carry);
    }
    if (carry != 0U)
    {
        a[count++] = carry;
    }
    return count;
}

// x = m * 10^e < 10^19 as Frac + 1 words, low word first, with word Frac as the integer part.
// For e < 0, x = (m * 10^b) * 10^(-9a) with 9a = b - e, and the table error is below 2^-300.
template <int Frac, typename Significand>
constexpr auto binary_words(const Significand m, const int e, const fx_span<std::uint64_t, static_cast<std::size_t>(Frac + 1)> x) noexcept -> void
{
    constexpr int cw {Frac + 2};
    const auto m128 {static_cast<boost::int128::uint128_t>(m)};
    if (e >= 0)
    {
        x[Frac] = m128.low * pow10(static_cast<std::uint64_t>(e));
        return;
    }

    const int a {(8 - e) / 9};
    std::uint64_t g[4] {m128.low, m128.high};
    const int count {times_word(fx_make_span(g), m128.high != 0U ? 2 : 1, pow10(static_cast<std::uint64_t>(9 * a + e)))};

    const auto& c {fx_table::pow10_neg9_bin[a - 1]};
    std::uint64_t p[static_cast<std::size_t>(cw + 3)] {};
    for (int i {}; i < count; ++i)
    {
        std::uint64_t carry {};
        for (int j {}; j < cw; ++j)
        {
            p[i + j] = mul_add(g[i], c[cw - 1 - j], p[i + j], carry);
        }
        p[i + cw] = carry;
    }
    for (int i {}; i <= Frac; ++i)
    {
        x[i] = p[cw - Frac + i];
    }
}

// One-word Payne-Hanek reduction, where m < 10^16 and pi/4 < x < 10^19 give e >= -16: x = m * 10^e, or m with
// 2/pi * 10^e, is an integer below 2^64, and three words of the constant leave an error below 2^-128.
template <int N, typename Significand>
constexpr auto trig_reduce_word(const Significand m, const int e, const bool xneg) noexcept -> trig_arg<N>
{
    const bool scaled {e < 0};
    const std::uint64_t x {static_cast<std::uint64_t>(m) * pow10(static_cast<std::uint64_t>(scaled ? 0 : e))};
    const auto c {fx_make_span<3>(scaled ? fx_table::two_over_pi_pow10_bin[-e - 1] : fx_table::two_over_pi_bin)};

    // q = x * c: q[3] is the integer part, and q[2] to q[0] are the fraction f
    std::uint64_t q[4] {};
    std::uint64_t carry {};
    for (int i {}; i < 3; ++i)
    {
        q[i] = mul_add(x, c[2 - i], 0U, carry);
    }
    q[3] = carry;

    // Above one half, go to the next quadrant with r = -(1 - f).
    unsigned n {static_cast<unsigned>(q[3] & 3U)};
    boost::int128::uint128_t f {q[1], q[0]};
    bool rneg {(q[2] >> 63U) != 0U};
    if (rneg)
    {
        f = ~f + 1U;
        q[2] = ~q[2] + (f == 0U ? 1U : 0U); // with the carry from f
        n = (n + 1U) & 3U;
    }
    negate_arg(xneg, n, rneg);
    trig_arg<N> out {};
    out.n = n;
    out.neg = rneg;

    // f = g * 2^-(128 + sh) with the top bit of g set, where sh <= 58 (see trig_traits), and
    // r = f * pi/2 = big * 2^-(255 + sh) with big = g * pi/4.
    const auto sh {static_cast<unsigned>(countl_zero(q[2]))};
    const boost::int128::uint128_t head {q[2], f.high};
    const auto g {sh == 0U ? head : (head << sh) | (f.low >> (64U - sh))};
    constexpr boost::int128::uint128_t pi4 {(UINT64_C(1) << 63U) | (fx_table::pio2_bin[1] >> 1U), (fx_table::pio2_bin[1] << 63U) | (fx_table::pio2_bin[2] >> 1U)};
    std::uint64_t big[4] {};
    mul_words(g, pi4, big);
    out.rf.w[0] = big[3] >> (1U + sh);
    out.fixed = true;

    // r = h * 2^-(128 + z) with the top bit of h set, and sig + tail / 2^64 = r * 10^k
    const bool top {(big[3] >> 63U) != 0U};
    const unsigned z {top ? sh - 1U : sh};
    const boost::int128::uint128_t high {big[3], big[2]};
    const auto h {top ? high : (high << 1U) | (big[1] >> 63U)};
    int k {sig_digits(1) + fx_floor_log10_2(static_cast<int>(z))};
    std::uint64_t p[4] {};
    mul_words(h, pow10(static_cast<boost::int128::uint128_t>(k)), p);
    out.sig = z == 0U ? p[2] : (p[2] >> z) | (p[3] << (64U - z));
    out.tail = z == 0U ? p[1] : (p[1] >> z) | (p[2] << (64U - z));
    fill_digit(out, sig_digits(1), k);
    out.k = k;
    return out;
}

// Word j of pi/4 in pw fraction words, low word first
constexpr auto quarter_pi_word(const int pw, const int j) noexcept -> std::uint64_t
{
    return (fx_table::pio2_bin[pw - j] >> 1U) | (fx_table::pio2_bin[pw - j - 1] << 63U);
}

// Cody-Waite reduction in words of 64 bits for |x| < 10^19, with pi/2 to 2^-64(bin_frac + 1). The relative errors
// of r are below 2^-64, 2^-128 and 2^-190, less than those of trig_reduce.
template <typename T, typename Significand>
constexpr auto trig_reduce_binary(const Significand m, const int e, const bool xneg) noexcept -> trig_arg<trig_traits<T>::words>
{
    constexpr int words {trig_traits<T>::words};
    constexpr int frac {trig_traits<T>::bin_frac};
    constexpr int pw {frac + 1};
    constexpr int rw {words + 1};
    BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
    {
        return trig_reduce_word<words>(m, e, xneg);
    }

    std::uint64_t x[static_cast<std::size_t>(frac + 1)] {};
    binary_words<frac>(m, e, fx_make_span(x));
    trig_arg<words> out {};

    // n = round(x * 2/pi) from the integer word and the top fraction word, with an error below 2^-61
    std::uint64_t n {};
    std::uint64_t f {boost::int128::detail::umul(x[frac], fx_table::two_over_pi_bin[0], n)};
    std::uint64_t h[2] {};
    boost::int128::detail::umul(x[frac], fx_table::two_over_pi_bin[1], h[0]);
    boost::int128::detail::umul(x[frac - 1], fx_table::two_over_pi_bin[0], h[1]);
    for (const auto t : h)
    {
        f += t;
        n += f < t ? 1U : 0U;
    }
    n += f >> 63U;

    // r = x - n * pi/2 in pw fraction words and an integer word, two's complement
    std::uint64_t v[static_cast<std::size_t>(pw + 1)] {};
    std::uint64_t carry {};
    std::uint64_t borrow {};
    for (int i {}; i <= pw; ++i)
    {
        const std::uint64_t a {i >= 1 ? x[i - 1] : 0U};
        const std::uint64_t t {mul_add(n, fx_table::pio2_bin[pw - i], 0U, carry)};
        v[i] = a - t - borrow;
        borrow = (a < t || (a == t && borrow != 0U)) ? 1U : 0U;
    }
    bool rneg {(v[pw] >> 63U) != 0U};
    if (rneg)
    {
        carry = 1U;
        for (int i {}; i <= pw; ++i)
        {
            v[i] = ~v[i] + carry;
            carry = (carry != 0U && v[i] == 0U) ? 1U : 0U;
        }
    }

    // Near pi/4 the estimate of n can be one off: then |r| > pi/4, and pi/2 - |r| is r on the other side.
    int top {pw - 1};
    while (top > 0 && v[top] == quarter_pi_word(pw, top))
    {
        --top;
    }
    if (v[top] > quarter_pi_word(pw, top))
    {
        borrow = 0U;
        for (int i {}; i <= pw; ++i)
        {
            const std::uint64_t p {fx_table::pio2_bin[pw - i]};
            const std::uint64_t t {v[i]};
            v[i] = p - t - borrow;
            borrow = (p < t || (p == t && borrow != 0U)) ? 1U : 0U;
        }
        n += rneg ? UINT64_MAX : 1U;
        rneg = !rneg;
    }
    unsigned q {static_cast<unsigned>(n & 3U)};
    negate_arg(xneg, q, rneg);
    out.n = q;
    out.neg = rneg;

    // r is not zero: its worst case has 27, 58 and 117 leading zero bits (see trig_traits).
    int lead {};
    while (lead < pw - 1 && v[pw - 1 - lead] == 0U)
    {
        ++lead;
    }

    // r = s * 2^-64(lead + rw), with rw + 1 words of r from word lead. The four spare
    // words take the scale by 10^k below.
    std::uint64_t s[static_cast<std::size_t>(rw + 6)] {};
    for (int j {}; j <= rw; ++j)
    {
        const int i {pw - lead - rw + j};
        s[j] = i >= 0 ? v[i] : 0U;
    }
    const auto r {fx_make_span(s)};

    // rf = r * 2^(64N - 2): drop lead + 1 words and 2 bits.
    for (int i {}; i < words; ++i)
    {
        const int w {i + lead + 1};
        const std::uint64_t lo {w <= rw ? r[w] : 0U};
        const std::uint64_t hi {w + 1 <= rw ? r[w + 1] : 0U};
        out.rf.w[i] = (lo >> 2U) | (hi << 62U);
    }
    out.fixed = true;

    // For r in [2^-(z + 1), 2^-z), sig = floor(r * 10^k) with k = D + floor(z * log10(2)) has D - 1 or D digits.
    int count {rw + 1};
    while (r[count - 1] == 0U)
    {
        --count;
    }
    const int z {64 * (lead + rw - count) + countl_zero(r[count - 1])};
    constexpr int digits {sig_digits(words)};
    int k {digits + fx_floor_log10_2(z)};
    int left {k};
    while (left >= 19)
    {
        count = times_word(r, count, pow10(static_cast<std::uint64_t>(19)));
        left -= 19;
    }
    if (left != 0)
    {
        count = times_word(r, count, pow10(static_cast<std::uint64_t>(left)));
    }
    const int unit {lead + rw};
    out.sig = boost::int128::uint128_t {unit + 1 < count ? r[unit + 1] : 0U, unit < count ? r[unit] : 0U};
    out.tail = r[unit - 1];
    fill_digit(out, digits, k);
    out.k = k;
    return out;
}

// floor(pi/4 * 10^19)
constexpr std::uint64_t quarter_pi_19 {UINT64_C(7853981633974483096)};

// x = sig * 10^-k with sig of D = sig_digits(words) digits; for |x| <= pi/4 this is r itself.
template <typename T>
constexpr auto trig_prepare(const T x) noexcept -> trig_arg<trig_traits<T>::words>
{
    constexpr int words {trig_traits<T>::words};

    // floor(pi/4 * 10^D), where 10^38 gives 78539816339744830961566084581987572104
    constexpr int sd {sig_digits(words)};
    constexpr boost::int128::uint128_t quarter_pi {sd == 19 ? boost::int128::uint128_t {UINT64_C(0), quarter_pi_19}
                                                            : boost::int128::uint128_t {UINT64_C(4257651975108193236), UINT64_C(12340784068543502728)}};

    trig_arg<words> out {};
    const bool neg {signbit(x)};
    int e {};
    const auto m {frexp10(neg ? -x : x, &e)};
    out.neg = neg;
    if (m == 0U)
    {
        out.zero = true;
        return out;
    }

    const auto digits {num_digits(m)};
    const auto sig {static_cast<boost::int128::uint128_t>(m) * pow10(static_cast<boost::int128::uint128_t>(sd - digits))};
    const int k {(sd - digits) - e};
    if (k > sd || (k == sd && sig <= quarter_pi))
    {
        out.sig = sig;
        out.k = k;
        return out;
    }
    // |x| < 10^digits * 10^e
    if (digits + e <= 19)
    {
        return trig_reduce_binary<T>(m, e, neg);
    }
    return trig_reduce<T>(m, e, neg);
}

// trig_prepare in one word for a first try, called only for decimal64_t: false for x = 0 and |x| >= 10^19. frexp10
// gives m of precision_v<T> digits.
template <typename T>
constexpr auto trig_prepare_word(const T x, trig_arg<1>& out) noexcept -> bool
{
    constexpr int sd {sig_digits(1)};
    constexpr int prec {precision_v<T>};
    const bool neg {signbit(x)};
    int e {};
    const auto m {frexp10(neg ? -x : x, &e)};
    if (m == 0U || prec + e > 19)
    {
        return false;
    }
    const auto sig {static_cast<std::uint64_t>(m) * pow10(static_cast<std::uint64_t>(sd - prec))};
    const int k {(sd - prec) - e};
    if (k > sd || (k == sd && sig <= quarter_pi_19))
    {
        out.sig = sig;
        out.k = k;
        out.neg = neg;
        return true;
    }
    out = trig_reduce_word<1>(m, e, neg);
    return true;
}

} // namespace trig
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_TRIG_REDUCE_HPP
