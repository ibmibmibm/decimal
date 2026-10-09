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
