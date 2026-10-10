// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_BIG_INT_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_BIG_INT_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {

// A signed integer of up to W words, low word first, where ok turns false once a result does not fit. The exact
// compares of pow, beta and the polynomials use it.
template <int W>
struct big_int
{
    std::uint64_t w[static_cast<std::size_t>(W)];
    int len;
    bool neg;
    bool ok;
};

template <int W>
constexpr auto big_trim(big_int<W>& r) noexcept -> void
{
    while (r.len > 0 && r.w[r.len - 1] == 0U)
    {
        --r.len;
    }
}

// a in r, by the words that a uses
template <int W>
constexpr auto big_copy(const big_int<W>& a, big_int<W>& r) noexcept -> void
{
    for (int i {}; i < a.len; ++i)
    {
        r.w[i] = a.w[i];
    }
    r.len = a.len;
    r.neg = a.neg;
    r.ok = a.ok;
}

template <int W>
constexpr auto big_mul_u(big_int<W>& a, const std::uint64_t f) noexcept -> void
{
    std::uint64_t carry {};
    for (int i {}; i < a.len; ++i)
    {
        const auto t {static_cast<boost::int128::uint128_t>(a.w[i]) * f + carry};
        a.w[i] = t.low;
        carry = t.high;
    }
    a.ok = a.ok && (carry == 0U || a.len < W);
    if (carry != 0U && a.ok)
    {
        a.w[a.len++] = carry;
    }
}

template <int W>
constexpr auto big_mul_pow10(big_int<W>& a, int j) noexcept -> void
{
    for (; j > 0 && a.ok; j -= 19)
    {
        big_mul_u(a, pow10(static_cast<std::uint64_t>(j < 19 ? j : 19)));
    }
}

// v * 10^j in r
template <int W>
constexpr auto big_set(big_int<W>& r, const boost::int128::uint128_t& v, const int j) noexcept -> void
{
    r.w[0] = v.low;
    r.w[1] = v.high;
    r.len = 2;
    r.neg = false;
    r.ok = true;
    big_trim(r);
    big_mul_pow10(r, j);
}

// a * b in r, which is neither
template <int W>
constexpr auto big_mul(const big_int<W>& a, const big_int<W>& b, big_int<W>& r) noexcept -> void
{
    r.ok = a.ok && b.ok && a.len + b.len <= W;
    r.len = r.ok ? a.len + b.len : 0;
    r.neg = a.neg != b.neg;
    for (int i {}; i < r.len; ++i)
    {
        r.w[i] = 0U;
    }
    for (int i {}; i < (r.ok ? a.len : 0); ++i)
    {
        std::uint64_t carry {};
        for (int j {}; j < b.len; ++j)
        {
            const auto t {static_cast<boost::int128::uint128_t>(a.w[i]) * b.w[j] + r.w[i + j] + carry};
            r.w[i + j] = t.low;
            carry = t.high;
        }
        r.w[i + b.len] = carry;
    }
    big_trim(r);
}

// The order of |a| and |b|: -1, 0 or 1
template <int W>
constexpr auto big_cmp(const big_int<W>& a, const big_int<W>& b) noexcept -> int
{
    if (a.len != b.len)
    {
        return a.len < b.len ? -1 : 1;
    }
    for (int i {a.len - 1}; i >= 0; --i)
    {
        if (a.w[i] != b.w[i])
        {
            return a.w[i] < b.w[i] ? -1 : 1;
        }
    }
    return 0;
}

// a + b, or a - b for sub, in r, which may be a
template <int W>
constexpr auto big_add(const big_int<W>& a, const big_int<W>& b, const bool sub, big_int<W>& r) noexcept -> void
{
    const bool bneg {b.neg != sub};
    const bool fits {a.len < W && b.len < W};
    const int len {!fits ? 0 : a.len > b.len ? a.len : b.len};
    const bool ok {a.ok && b.ok && fits};
    const bool swap {a.neg != bneg && big_cmp(a, b) < 0};
    const auto& x {swap ? b : a};
    const auto& y {swap ? a : b};
    const bool neg {swap ? bneg : a.neg};
    std::uint64_t carry {};
    for (int i {}; i < len; ++i)
    {
        const std::uint64_t xi {i < x.len ? x.w[i] : 0U};
        const std::uint64_t yi {i < y.len ? y.w[i] : 0U};
        if (a.neg == bneg)
        {
            const auto t {static_cast<boost::int128::uint128_t>(xi) + yi + carry};
            r.w[i] = t.low;
            carry = t.high;
        }
        else
        {
            r.w[i] = xi - yi - carry;
            carry = (xi < yi || (xi == yi && carry != 0U)) ? 1U : 0U;
        }
    }
    r.w[len] = a.neg == bneg ? carry : 0U;
    r.len = len + 1;
    big_trim(r);
    r.neg = neg;
    r.ok = ok;
}

// v^n in r, with t for the products
template <int W>
constexpr auto big_pow(const boost::int128::uint128_t& v, unsigned n, big_int<W>& r, big_int<W>& t) noexcept -> void
{
    big_int<W> base {};
    big_set(base, v, 0);
    big_set(r, 1U, 0);
    while (n != 0U && r.ok && base.ok)
    {
        if ((n & 1U) != 0U)
        {
            big_mul(r, base, t);
            big_copy(t, r);
        }
        n >>= 1U;
        if (n != 0U)
        {
            big_mul(base, base, t);
            big_copy(t, base);
        }
    }
    r.ok = r.ok && base.ok;
}

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_BIG_INT_HPP
