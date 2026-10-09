// Copyright 2023 Matt Borland
// Copyright 2023 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_ATANH_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_ATANH_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace log_detail {

// atanh(x) = ln(a / b) / 2 for 1/256 <= |x| < 1 with a = 1 + |x| and b = 1 - |x| exact. With a = ua * 10^na and
// b = ub * 10^nb, and va = ua * 2^sa and vb = ub * 2^sb in [1, 2), a / b = (va / vb) * 2^(sb - sa) * 10^(na - nb).
template <int N>
constexpr auto atanh_r(const boost::int128::uint128_t& sig, const int e) noexcept -> fx_wide<N>
{
    const auto one {pow10(static_cast<boost::int128::uint128_t>(-e))};
    const auto a {one + sig};
    const auto b {one - sig};
    const int na {num_digits(a)};
    const int nb {num_digits(b)};
    unsigned sa {};
    unsigned sb {};
    const auto va {fx_normalize(fx_from_decimal<N>(a, na), sa)};
    const auto vb {fx_normalize(fx_from_decimal<N>(b, nb), sb)};
    auto v {fx_div(va, vb)};
    int n2 {static_cast<int>(sb) - static_cast<int>(sa)};
    if (v.w[N - 1] < (UINT64_C(1) << 62U))
    {
        v = fx_shl(v, 1U);
        --n2;
    }
    auto r {log_wide_of(v, n2, na - nb, log_base::e)};

    // r / 2, for r > 0
    for (int i {}; i < N; ++i)
    {
        r.w[i] = (r.w[i] >> 1U) | (r.w[i + 1] << 63U);
    }
    r.w[N] >>= 1U;
    return r;
}

// atanh of a finite x != 0, rounded once in the current mode
template <typename T>
constexpr auto atanh_fixed(const T x) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    int e {};
    const auto sig {static_cast<boost::int128::uint128_t>(frexp10(x, &e))};
    const bool neg {signbit(x)};

    // |x| >= 1 has e + p > 0: atanh(+-1) = +-infinity, and NaN above, also for a non-canonical x such as
    // 10000001 * 10^-7 that frexp10 rounds to 1
    if (e + p > 0)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        if (e + p == 1 && sig == pow10(static_cast<boost::int128::uint128_t>(p - 1)) && abs(x) == T {1})
        {
            return neg ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
        }
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }

    // For |x| < 1/256, atanh(x) = x * atanh(x) / x, and atanh(x) / x >= 1 gives digits >= m
    if (decimal_tiny(sig, e))
    {
        const auto u {fx_from_decimal<words>(sig, -e)};
        const auto m {sig * pow10(static_cast<boost::int128::uint128_t>(scale10 - p))};
        const auto z {fx_mul(u, u)};
        const auto d {fx_scale(m, atanh_series(z))};
        BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
        {
            // Next to a rounding boundary, atanh(x) = m (1 + x^2 R) with R = 1/3 + x^2/5 + ..., and the rest t^2 R 10^-2a
            // takes its own precision
            if (fx_near_boundary(d, 2U))
            {
                int a {};
                const auto t {fx_lead_t<words>(sig, -e, a)};
                const auto g {fx_mul(fx_mul(t, t), fx_horner_range<false>(z, fx_table::inv_odd, 1U, 11U, fx<words> {}))};
                return fx_round<T>(static_cast<digits_type>(fx_lead_digits(m, g, 2 * a, false)), e - (scale10 - p), neg);
            }
        }
        return fx_round<T>(static_cast<digits_type>(d), e - (scale10 - p), neg);
    }

    BOOST_DECIMAL_IF_CONSTEXPR (words > 1)
    {
        // A first try in one word fewer for |x| >= 0.1, with an error below 128 units of its last bit, or in
        // (-10.1, 5.1) units of 2^-62 in one word
        T r {};
        if (decimal_at_least(sig, e, p, 10U, 0) && fx_wide_try_round(atanh_r<words == 3 ? 2 : 1>(sig, e), neg, r, words == 3 ? 7U : 4U))
        {
            return r;
        }
    }
    return fx_wide_round<T, words>(atanh_r<words>(sig, e), neg);
}

} // namespace log_detail

template <typename T>
constexpr auto atanh_impl(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    T result { };

    const auto fpc {fpclassify(x)};

    if (fpc != FP_NORMAL && fpc != FP_SUBNORMAL)
    {
        // An infinity is outside [-1, 1]
        result = isinf(x) ? std::numeric_limits<T>::quiet_NaN() : x;
    }
    else
    {
        result = log_detail::atanh_fixed(x);
    }

    return result;
}

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto atanh(const T x) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    using evaluation_type = detail::evaluation_type_t<T>;

    return static_cast<T>(detail::atanh_impl(static_cast<evaluation_type>(x)));
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_ATANH_HPP
