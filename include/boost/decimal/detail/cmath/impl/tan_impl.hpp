// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_TAN_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_TAN_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/u256.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/impl/trig_reduce.hpp>
#include <boost/decimal/detail/cmath/impl/sin_impl.hpp>
#include <boost/decimal/detail/cmath/impl/cos_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace trig {

// The digits of tan(r) for an even quadrant, else of -cot(r), from s = sin(r) / r and c = cos(r).
template <int N>
constexpr auto tan_digits(const trig_arg<N>& r, const fx<N>& s, const fx<N>& c) noexcept -> boost::int128::uint128_t
{
    if ((r.n & 1U) == 0U)
    {
        // tan(r) = r * (sin(r) / r) / cos(r), and the factor is above 1 for any r != 0.
        BOOST_DECIMAL_IF_CONSTEXPR (N == 1)
        {
            // One word: floor(sig * s / c) in one division, below 1.28 * sig < 2^64. Without the tail, all
            // decimal32_t results stay the same.
            const auto num {static_cast<boost::int128::uint128_t>(s.w[0]) * r.sig.low};
            std::uint64_t rem {};
            const auto q {boost::int128::detail::udiv_2by1(num.high, num.low, c.w[0], rem)};
            return q > r.sig.low ? q : r.sig.low;
        }
        auto m {fx_div(s, c)};
        if (!fx_less(fx_one<N>(), m))
        {
            m = fx_one<N>();
            m.w[0] += 1U;
        }
        return fx_scale(r.sig, r.tail, m);
    }

    // -cot(r) = -(cos(r) / (sin(r) / r)) / r. One word: c * 10^37 / (s * (sig + tail / 2^64)), with sig of 19 digits,
    // has 18 or 19 digits. Both shift left by sh = countl_zero(den.high), 2 <= sh <= 6, to set the top bit of den.
    BOOST_DECIMAL_IF_CONSTEXPR (N == 1)
    {
        constexpr auto ten37 {pow10(static_cast<boost::int128::uint128_t>(37))};
        const auto lo {static_cast<boost::int128::uint128_t>(c.w[0]) * ten37.low};
        const auto hi {static_cast<boost::int128::uint128_t>(c.w[0]) * ten37.high + lo.high};
        const boost::int128::uint128_t s0 {s.w[0]};
        const auto den {s0 * r.sig.low + ((s0 * r.tail) >> 64U)};
        const auto sh {static_cast<unsigned>(countl_zero(den.high))};
        std::uint64_t rem {};
        return boost::int128::detail::udiv_2by1((hi.high << sh) | (hi.low >> (64U - sh)), (hi.low << sh) | (lo.low >> (64U - sh)), (den << sh).high, rem);
    }

    // cos(r) / (sin(r) / r) is below 1 for any r != 0.
    const auto m {fx_clamp_below_one(fx_mul(c, fx_recip(s)))};

    // floor(m * 10^75 / 2^(64N - 2)), then divide by sig: for m >= 0.7, the quotient has 37 or 38 digits.
    std::uint64_t p[static_cast<std::size_t>(N + 4)] {};
    for (int i {}; i < N; ++i)
    {
        std::uint64_t carry {};
        for (int j {}; j < 4; ++j)
        {
            const auto t {static_cast<boost::int128::uint128_t>(m.w[i]) * fx_table::pow10_75[j] + p[i + j] + carry};
            p[i + j] = t.low;
            carry = t.high;
        }
        p[i + 4] = carry;
    }
    std::uint64_t num[4] {};
    for (int i {}; i < 4; ++i)
    {
        num[i] = (p[N - 1 + i] >> 62U) | (p[N + i] << 2U);
    }
    const auto dm {impl::div_mod(u256 {num[3], num[2], num[1], num[0]}, r.sig)};
    auto q {boost::int128::uint128_t {dm.quotient.bytes[1], dm.quotient.bytes[0]}};

    // num / (sig + tail / 2^64) = q - (q * tail / 2^64 - rem) / sig, to far below one unit, and
    // q * tail / 2^64 is below 10 * sig: subtract 1 for each sig in it, rounded up. Only three words need it.
    BOOST_DECIMAL_IF_CONSTEXPR (N == 3)
    {
        const boost::int128::uint128_t rem {dm.remainder.bytes[1], dm.remainder.bytes[0]};
        std::uint64_t carry {};
        boost::int128::detail::umul(q.low, r.tail, carry);
        auto excess {static_cast<boost::int128::uint128_t>(q.high) * r.tail + carry};
        while (excess > rem)
        {
            --q;
            excess -= excess - rem > r.sig ? r.sig : excess - rem;
        }
    }
    return q;
}

// tan(r) for an even quadrant, else -cot(r), for r != 0. Near a rounding boundary the Taylor series
// decide, as in sin_cos_of.
template <typename T>
constexpr auto tan_of(const trig_arg<trig_traits<T>::words>& r) noexcept -> T
{
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    auto v {tan_digits(r, sinc_fx(z), cos_fx(z))};
    BOOST_DECIMAL_IF_CONSTEXPR (trig_traits<T>::words == 3)
    {
        if (z.w[2] != 0U && fx_near_boundary(v))
        {
            const auto rt {fixed_r_with_tail(r)};
            const auto zt {fx_mul(rt, rt)};
            v = tan_digits(r, trig_taylor(zt, false), trig_taylor(zt, true));
        }
    }
    constexpr int digits {sig_digits(trig_traits<T>::words)};
    return (r.n & 1U) == 0U ? trig_round<T>(v, -r.k, r.neg) : trig_round<T>(v, r.k - (2 * digits - 1), !r.neg);
}

// tan_of in one word for decimal64_t with its polynomials: the error of v is in (-16.6, 13.9) units for tan(r) and in
// (-9.6, 9.9) for -cot(r), thus w = 17 and 10. It rounds once in out where v cannot be next to a rounding boundary.
template <typename T>
constexpr auto tan_try(const trig_arg<1>& r, T& out) noexcept -> bool
{
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    const auto v {tan_digits(r, sinc_fx(z, sin_poly_top(z)), cos_fx(z, cos_poly_top(z)))};
    constexpr int digits {sig_digits(1)};
    return (r.n & 1U) == 0U ? trig_try_round(v.low, -r.k, r.neg, 17U, out) : trig_try_round(v.low, r.k - (2 * digits - 1), !r.neg, 10U, out);
}

} // namespace trig
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_TAN_IMPL_HPP
