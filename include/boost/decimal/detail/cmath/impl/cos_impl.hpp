// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>
#include <boost/decimal/detail/cmath/impl/trig_reduce.hpp>
#include <boost/decimal/detail/cmath/impl/sin_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace trig {

constexpr auto cos_poly(const fx<1>& z) noexcept -> fx<1>
{
    return fx_alternating(z, fx_table::cos_d32);
}

constexpr auto cos_poly(const fx<2>& z) noexcept -> fx<2>
{
    return fx_alternating_split<fx_table::cos_d64_shift>(z, fx_table::cos_d64_head, fx_table::cos_d64_tail);
}

constexpr auto cos_poly(const fx<3>& z) noexcept -> fx<3>
{
    return fx_alternating_split<fx_table::cos_d128_mid_shift, fx_table::cos_d128_tail_shift>(z, fx_table::cos_d128_head, fx_table::cos_d128_mid, fx_table::cos_d128_tail);
}

// cos(r) for z = r^2, which is below 1 for any r != 0.
template <int N>
constexpr auto cos_fx(const fx<N>& z) noexcept -> fx<N>
{
    fx<N> half {};
    half.w[N - 1] = UINT64_C(1) << 61U;
    const auto inner {fx_sub(half, fx_mul(z, cos_poly(z)))};
    return fx_clamp_below_one(fx_sub(fx_one<N>(), fx_mul(z, inner)));
}

// sin(r) / r = 1 - z/(2*3) * (1 - z/(4*5) * (1 - ...)), or with use_cos cos(r) = 1 - z/(1*2) * (1 - ...),
// to 22 terms with an error below 2^-186 for z <= (pi/4)^2. It is slow, but it needs no table.
template <int N>
constexpr auto trig_taylor(const fx<N>& z, const bool use_cos) noexcept -> fx<N>
{
    const std::uint64_t c {use_cos ? 1U : 0U};
    auto s {fx_one<N>()};
    for (std::uint64_t k {22}; k > 0U; --k)
    {
        s = fx_sub(fx_one<N>(), fx_div_small(fx_mul(z, s), (2U * k - c) * (2U * k + 1U - c)));
    }
    return fx_clamp_below_one(s);
}

// sin(|r|), or with use_cos cos(r), with the sign neg, for r != 0. For three words the polynomials are good
// to 2 units of the 38 digits, thus near a rounding boundary the Taylor series decides.
template <typename T>
constexpr auto sin_cos_of(const trig_arg<trig_traits<T>::words>& r, const bool use_cos, const bool neg) noexcept -> T
{
    constexpr int digits {sig_digits(trig_traits<T>::words)};
    constexpr auto one {pow10(static_cast<boost::int128::uint128_t>(digits))};
    const auto rf {fixed_r(r)};
    const auto z {fx_mul(rf, rf)};
    auto v {use_cos ? fx_scale(one, cos_fx(z)) : fx_scale(r.sig, r.tail, sinc_fx(z))};

    // Below z = 2^-62 the polynomials are as good as the series.
    BOOST_DECIMAL_IF_CONSTEXPR (trig_traits<T>::words == 3)
    {
        if (z.w[2] != 0U && fx_near_boundary(v))
        {
            const auto rt {fixed_r_with_tail(r)};
            const auto f {trig_taylor(fx_mul(rt, rt), use_cos)};
            v = use_cos ? fx_scale(one, f) : fx_scale(r.sig, r.tail, f);
        }
    }
    return trig_round<T>(v, use_cos ? -digits : -r.k, neg);
}

} // namespace trig
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_COS_IMPL_HPP
