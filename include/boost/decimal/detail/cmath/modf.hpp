// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_MODF_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_MODF_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/write_payload.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/decompose.hpp>
#include <boost/decimal/detail/add_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <type_traits>
#include <cmath>
#endif

namespace boost {
namespace decimal {

namespace detail {

// Divides sig by 10^shift. A 128-bit sig that fits in 64 bits divides fastest as it is,
// a wider one by the reciprocal of 10^shift
template <typename U>
BOOST_DECIMAL_CUDA_CONSTEXPR auto split_pow10(const U sig, const int shift) noexcept
{
    return impl::divmod_pow10_dispatch(sig, shift, pow10(static_cast<U>(shift)));
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto split_pow10(const int128::uint128_t sig, const int shift) noexcept
{
    const auto p {pow10(static_cast<int128::uint128_t>(shift))};
    return sig.high == 0U ? impl::divmod(sig, p) : impl::divmod_pow10_dispatch(sig, shift, p);
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto modf(const T x, T* iptr) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using sig_type = typename T::significand_type;

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (isnan(x))
    {
        *iptr = detail::nan_result(x);
        return *iptr;
    }
    if (isinf(x))
    {
        *iptr = x;
        return T {0U, 0, signbit(x)};
    }
    #endif

    // Both parts have the sign of x. In the IEEE types, the integer part has the exponent max(q(x), 0),
    // and the fractional part keeps q(x). The fast types can normalize both parts.
    const auto c {decompose(x)};
    if (c.exp >= 0)
    {
        *iptr = x;
        return detail::pack_in_range<T>(static_cast<sig_type>(0U), c.exp, c.sign);
    }

    const auto shift {-static_cast<int>(c.exp)};
    auto integer_part {static_cast<sig_type>(0U)};
    sig_type fraction {c.sig};
    // A non canonical significand from from_bits can have precision + 1 digits
    if (shift < detail::precision_v<T> || (shift == detail::precision_v<T> && c.sig > detail::max_significand_v<T>))
    {
        const auto dr {detail::split_pow10(c.sig, shift)};
        integer_part = dr.quotient;
        fraction = dr.remainder;
    }

    *iptr = detail::pack_in_range<T>(integer_part, 0, c.sign);
    return detail::pack_in_range<T>(fraction, c.exp, c.sign);
}

} // namespace decimal
} // namespace boost

#endif //BOOST_DECIMAL_DETAIL_CMATH_MODF_HPP
