// Copyright 2023-2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_SIN_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_SIN_IMPL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/cmath/impl/fixed_point.hpp>

namespace boost {
namespace decimal {
namespace detail {
namespace trig {

constexpr auto sin_poly(const fx<1>& z) noexcept -> fx<1>
{
    return fx_alternating(z, fx_table::sin_d32);
}

constexpr auto sin_poly(const fx<2>& z) noexcept -> fx<2>
{
    return fx_alternating_split<fx_table::sin_d64_shift>(z, fx_table::sin_d64_head, fx_table::sin_d64_tail);
}

constexpr auto sin_poly(const fx<3>& z) noexcept -> fx<3>
{
    return fx_alternating_split<fx_table::sin_d128_mid_shift, fx_table::sin_d128_tail_shift>(z, fx_table::sin_d128_head, fx_table::sin_d128_mid, fx_table::sin_d128_tail);
}

// The polynomial of decimal64_t in one word
constexpr auto sin_poly_top(const fx<1>& z) noexcept -> fx<1>
{
    return fx_alternating_top<fx_table::sin_d64_shift>(z, fx_table::sin_d64_head, fx_table::sin_d64_tail);
}

// The polynomial of decimal128_t in two words
constexpr auto sin_poly_top(const fx<2>& z) noexcept -> fx<2>
{
    return fx_alternating_top<fx_table::sin_d128_mid_shift, fx_table::sin_d128_tail_shift>(z, fx_table::sin_d128_head, fx_table::sin_d128_mid, fx_table::sin_d128_tail);
}

// 16 (1 - sin(r) / r) for z = r^2 in one word, with the decimal64_t polynomial at 16 times its scale. For z below r^2
// by less than 2.6 units, the error is in (-9.2, 1) units.
constexpr auto sin_rest_top(const fx<1>& z) noexcept -> fx<1>
{
    fx<1> s {{fx_tail_sum(z.w[0], fx_table::sin_d64_tail) >> (fx_table::sin_d64_shift - 4)}};
    s = fx_sub(fx_top_shl(fx_table::sin_d64_head[1], 4U), fx_mul(z, s));
    s = fx_sub(fx_top_shl(fx_table::sin_d64_head[0], 4U), fx_mul(z, s));
    return fx_mul(z, s);
}

// sin(r) / r for z = r^2 and the polynomial p, which is below 1 for any r != 0.
template <int N>
constexpr auto sinc_fx(const fx<N>& z, const fx<N>& p) noexcept -> fx<N>
{
    return fx_clamp_below_one(fx_sub(fx_one<N>(), fx_mul(z, p)));
}

template <int N>
constexpr auto sinc_fx(const fx<N>& z) noexcept -> fx<N>
{
    return sinc_fx(z, sin_poly(z));
}

} // namespace trig
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_SIN_IMPL_HPP
