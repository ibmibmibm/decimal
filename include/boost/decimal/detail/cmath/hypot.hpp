// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_HYPOT_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_HYPOT_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/promotion.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/u256.hpp>
#include <boost/decimal/detail/cmath/frexp10.hpp>
#include <boost/decimal/detail/cmath/fabs.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <cmath>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace hypot_detail {

// Bits [pos, pos + 64) of the K words of t, with zeros below bit 0
template <std::size_t K>
constexpr auto hypot_bits(const std::uint64_t (&t)[K], const int pos) noexcept -> std::uint64_t
{
    if (pos <= -64)
    {
        return 0U;
    }
    if (pos < 0)
    {
        return t[0] << static_cast<unsigned>(-pos);
    }
    const auto i {static_cast<std::size_t>(pos / 64)};
    const auto b {static_cast<unsigned>(pos % 64)};
    const std::uint64_t lo {i < K ? t[i] : 0U};
    const std::uint64_t hi {i + 1U < K ? t[i + 1U] : 0U};
    return b == 0U ? lo : (lo >> b) | (hi << (64U - b));
}

// floor(sqrt(t)) within a few units, for t > 0 in K words: t = m * 4^s with m in [1, 4), and sqrt(t) = m / sqrt(m)
// * 2^s, with the Newton steps of 1 / sqrt(m) in one word, then in N words
template <int N, int Steps, std::size_t K>
BOOST_DECIMAL_FORCE_INLINE constexpr auto hypot_sqrt(const std::uint64_t (&t)[K]) noexcept -> boost::int128::uint128_t
{
    std::size_t top {K - 1U};
    while (t[top] == 0U)
    {
        --top;
    }
    const int s {(64 * static_cast<int>(top) + 63 - countl_zero(t[top])) / 2};
    fx<N> m {};
    for (int i {}; i < N; ++i)
    {
        m.w[i] = hypot_bits(t, 2 * s - (64 * N - 2) + 64 * i);
    }
    const auto r {fx_mul(m, fx_rsqrt<N, Steps>(m))};
    return {hypot_bits(r.w, 64 * N - 2 - s + 64), hypot_bits(r.w, 64 * N - 2 - s)};
}

constexpr auto hypot_words(const std::uint64_t& t, std::uint64_t (&w)[1]) noexcept -> void
{
    w[0] = t;
}

constexpr auto hypot_words(const boost::int128::uint128_t& t, std::uint64_t (&w)[2]) noexcept -> void
{
    w[0] = t.low;
    w[1] = t.high;
}

constexpr auto hypot_words(const u256& t, std::uint64_t (&w)[4]) noexcept -> void
{
    for (int i {}; i < 4; ++i)
    {
        w[i] = t[static_cast<std::size_t>(i)];
    }
}

// q^2 in the wide type W of the sum
template <typename W>
constexpr auto hypot_square(const std::uint64_t q) noexcept -> W
{
    return static_cast<W>(q) * q;
}

template <typename W>
constexpr auto hypot_square(const boost::int128::uint128_t& q) noexcept -> W
{
    std::uint64_t r[4] {};
    mul_words(q, q, r);
    return {r[3], r[2], r[1], r[0]};
}

// sqrt(x^2 + y^2 + z^2) for finite values, rounded once: with x = sx * 10^ex and the largest exponent e, the integer
// t = sum sx^2 * 10^(4 - 2(e - ex)) < 3 * 10^(2p + 4) gives floor(sqrt(t)) * 10^(e - 2) and a sticky digit.
template <typename T>
constexpr auto hypot_fixed(const T x, const T y, const T z) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    using wide_type = std::conditional_t<(p <= 7), std::uint64_t, std::conditional_t<(p <= 16), boost::int128::uint128_t, u256>>;
    using root_type = std::conditional_t<(p <= 16), std::uint64_t, boost::int128::uint128_t>;
    using digits_type = fx_digits_t<fx_words(p)>;
    using sig_type = typename T::significand_type;
    constexpr std::size_t K {sizeof(wide_type) / sizeof(std::uint64_t)};

    const T v[3] {x, y, z};
    typename T::significand_type sig[3] {};
    int e[3] {};
    int top {(std::numeric_limits<int>::min)()};
    int nonzero {};
    int last {};
    for (int i {}; i < 3; ++i)
    {
        const auto s {frexp10(v[i], &e[i])};
        sig[i] = s;
        if (s != 0U)
        {
            top = e[i] > top ? e[i] : top;
            ++nonzero;
            last = i;
        }
    }
    // All zero gives 0, and one argument that is not zero gives its absolute value exactly
    if (nonzero <= 1)
    {
        return nonzero == 0 ? T {0} : fabs(v[last]);
    }

    wide_type t {};
    // The parts rem / 10^cut that t drops, at most two since the largest term is exact
    wide_type rem[2] {};
    int cut[2] {};
    int n {};
    for (int i {}; i < 3; ++i)
    {
        if (sig[i] == 0U)
        {
            continue;
        }
        const int shift {2 * (top - e[i]) - 4};
        // A shift of -2 or -4 scales the significand by 10 or 100 before the square, below 10^(p + 2)
        auto sq {hypot_square<wide_type>(shift < 0 ? sig[i] * pow10(static_cast<sig_type>(-shift / 2)) : sig[i])};
        if (shift > 2 * p)
        {
            rem[n] = sq;
            cut[n++] = shift;
            continue;
        }
        else if (shift > 0)
        {
            const auto dr {impl::divmod_pow10_dispatch(sq, shift, pow10(static_cast<wide_type>(static_cast<std::uint64_t>(shift))))};
            sq = dr.quotient;
            rem[n] = dr.remainder;
            cut[n++] = shift;
        }
        t = t + sq;
    }

    std::uint64_t w[K] {};
    hypot_words(t, w);
    auto q {static_cast<root_type>(hypot_sqrt<p <= 16 ? 1 : 2, p <= 7 ? 2 : 3>(w))};
    while (hypot_square<wide_type>(q) > t)
    {
        --q;
    }
    while (!(t < hypot_square<wide_type>(q + 1U)))
    {
        ++q;
    }
    bool sticky {rem[0] != 0U || rem[1] != 0U || hypot_square<wide_type>(q) != t};

    // Two dropped parts can add up to 1, which changes q when t + 1 = (q + 1)^2. That holds when rem[b] / 10^(cut[b] -
    // cut[a]) >= d for d = 10^cut[a] - rem[a].
    const int a {cut[0] < cut[1] ? 0 : 1};
    const int b {1 - a};
    const int k {cut[b] - cut[a]};
    if (rem[0] != 0U && rem[1] != 0U && cut[a] <= 2 * p && k <= 2 * p && hypot_square<wide_type>(q + 1U) - t == 1U)
    {
        const auto d {pow10(static_cast<wide_type>(static_cast<std::uint64_t>(cut[a]))) - rem[a]};
        auto r {rem[b]};
        bool exact {true};
        if (k > 0)
        {
            const auto dr {impl::divmod_pow10_dispatch(r, k, pow10(static_cast<wide_type>(static_cast<std::uint64_t>(k))))};
            r = dr.quotient;
            exact = dr.remainder == 0U;
        }
        if (!(r < d))
        {
            ++q;
            sticky = !(exact && r == d);
        }
    }

    // q has p + 2 or p + 3 digits. fx_round takes 37 or more digits (17 for one word) that truncate a value which is
    // not exact, here in (q, q + 1), and an exact q rounds in the constructor.
    if (!sticky)
    {
        return T {q, top - 2};
    }
    constexpr int scale {p <= 7 ? 8 : 37 - p - 2};
    return fx_round<T>(static_cast<digits_type>(q) * pow10(static_cast<std::uint64_t>(scale)), top - 2 - scale, false);
}

} // namespace hypot_detail

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
constexpr auto hypot_impl(const T x, const T y, const T z) noexcept
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isinf(x) || isinf(y) || isinf(z))
    {
        return std::numeric_limits<T>::infinity();
    }
    else if (isnan(x))
    {
        return issignaling(x) ? nan_conversion(x) : x;
    }
    else if (isnan(y))
    {
        return issignaling(y) ? nan_conversion(y) : y;
    }
    else if (isnan(z))
    {
        return issignaling(z) ? nan_conversion(z) : z;
    }
    #endif

    return hypot_detail::hypot_fixed(x, y, z);
}

template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T>
constexpr auto hypot_impl(const T x, const T y) noexcept
{
    return hypot_impl(x, y, T {0});
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T1, typename T2>
constexpr auto hypot(const T1 x, const T2 y) noexcept
    BOOST_DECIMAL_REQUIRES_TWO(detail::is_decimal_floating_point_v, T1, detail::is_decimal_floating_point_v, T2)
{
    #if BOOST_DECIMAL_DEC_EVAL_METHOD == 0

    using evaluation_type = detail::promote_args_t<T1, T2>;

    #elif BOOST_DECIMAL_DEC_EVAL_METHOD == 1

    using evaluation_type = detail::promote_args_t<T1, T2, decimal64_t>;

    #else // BOOST_DECIMAL_DEC_EVAL_METHOD == 2

    using evaluation_type = detail::promote_args_t<T1, T2, decimal128_t>;

    #endif

    using return_type = detail::promote_args_t<T1, T2>;

    return static_cast<return_type>(detail::hypot_impl(static_cast<evaluation_type>(x), static_cast<evaluation_type>(y)));
}

BOOST_DECIMAL_EXPORT template <typename T1, typename T2, typename T3>
constexpr auto hypot(const T1 x, const T2 y, const T3 z) noexcept
    BOOST_DECIMAL_REQUIRES_THREE(detail::is_decimal_floating_point_v, T1, detail::is_decimal_floating_point_v, T2, detail::is_decimal_floating_point_v, T3)
{
    #if BOOST_DECIMAL_DEC_EVAL_METHOD == 0

    using evaluation_type = detail::promote_args_t<T1, T2, T3>;

    #elif BOOST_DECIMAL_DEC_EVAL_METHOD == 1

    using evaluation_type = detail::promote_args_t<T1, T2, T3, decimal64_t>;

    #else // BOOST_DECIMAL_DEC_EVAL_METHOD == 2

    using evaluation_type = detail::promote_args_t<T1, T2, T3, decimal128_t>;

    #endif

    using return_type = detail::promote_args_t<T1, T2, T3>;

    return static_cast<return_type>(detail::hypot_impl(static_cast<evaluation_type>(x),
                                                       static_cast<evaluation_type>(y),
                                                       static_cast<evaluation_type>(z)));
}

} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_HYPOT_HPP
