// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_IMPL_ROUND_INTEGRAL_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_IMPL_ROUND_INTEGRAL_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/cmath/fpclassify.hpp>
#include <boost/decimal/detail/cmath/decompose.hpp>
#include <boost/decimal/detail/add_impl.hpp>
#include <boost/decimal/cfenv.hpp>
#include <boost/decimal/detail/concepts.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {
namespace detail {

// Rounds c.sig * 10^c.exp, with c.exp < 0, to an integer in the mode round:
// one division by 10^-c.exp
template <typename T, typename Components>
BOOST_DECIMAL_FORCE_INLINE constexpr auto round_integral_sig(const Components& c, const rounding_mode round) noexcept -> typename T::significand_type
{
    using sig_type = typename T::significand_type;

    const auto shift {-static_cast<int>(c.exp)};
    sig_type q {0U};
    bool inc {};

    if (shift > precision_v<T>)
    {
        // 0 < |x| < 0.1
        inc = round == (c.sign ? rounding_mode::fe_dec_downward : rounding_mode::fe_dec_upward);
    }
    else
    {
        const auto p {pow10(static_cast<sig_type>(shift))};
        q = static_cast<sig_type>(c.sig / p);
        const auto r {static_cast<sig_type>(c.sig - q * p)};
        if (r == 0U)
        {
            return q;
        }

        switch (round)
        {
            case rounding_mode::fe_dec_upward:
                inc = !c.sign;
                break;
            case rounding_mode::fe_dec_downward:
                inc = c.sign;
                break;
            case rounding_mode::fe_dec_to_nearest_from_zero:
                inc = r >= p - r;
                break;
            case rounding_mode::fe_dec_to_nearest:
                inc = r > p - r || (r == p - r && (static_cast<std::uint32_t>(q) & 1U) == 1U);
                break;
            case rounding_mode::fe_dec_toward_zero:
            default:
                break;
        }
    }

    if (inc)
    {
        ++q;
    }

    return q;
}

// A fast type has a full significand, so q * 10^-exp is full or a carry to 10^precision, as in ceil(9.999999).
// A q of 1 from a smaller exponent takes min_exp. So direct_pack_biased stores an exponent in [min_exp, 0].
template <typename T, typename Components, std::enable_if_t<is_fast_type_v<T>, bool> = true>
BOOST_DECIMAL_FORCE_INLINE constexpr auto round_integral_pack(const typename T::significand_type q, const Components& c) noexcept -> T
{
    using sig_type = typename T::significand_type;
    constexpr int min_exp {1 - precision_v<T>};

    sig_type sig {0U};
    int exp {0};
    if (q != 0U)
    {
        exp = static_cast<int>(c.exp) < min_exp ? min_exp : static_cast<int>(c.exp);
        sig = static_cast<sig_type>(q * pow10(static_cast<sig_type>(-exp)));
        if (BOOST_DECIMAL_UNLIKELY(sig == pow10(static_cast<sig_type>(precision_v<T>))))
        {
            sig = pow10(static_cast<sig_type>(precision_v<T> - 1));
            ++exp;
        }
    }
    return direct_pack_biased(sig, exp + bias_v<T>, c.sign);
}

template <typename T, typename Components, std::enable_if_t<!is_fast_type_v<T>, bool> = true>
BOOST_DECIMAL_FORCE_INLINE constexpr auto round_integral_pack(const typename T::significand_type q, const Components& c) noexcept -> T
{
    return pack_in_range<T>(q, 0, c.sign);
}

// Rounds x to an integer in the mode round. An IEEE result has exponent 0, and a zero result
// keeps the sign of x.
template <typename T>
constexpr auto round_integral_impl(const T x, const rounding_mode round) noexcept -> T
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(x))
    {
        return x;
    }
    #endif

    const auto c {decompose(x)};
    if (c.exp >= 0 || c.sig == 0U)
    {
        return x;
    }

    return round_integral_pack<T>(round_integral_sig<T>(c, round), c);
}

template <BOOST_DECIMAL_INTEGRAL Int>
constexpr auto round_integral_saturate(const bool sign) noexcept -> Int
{
    return sign ? (std::numeric_limits<Int>::min)() : (std::numeric_limits<Int>::max)();
}

// The full product of a and b, as the high and the low 64 bits
constexpr auto umul_wide(const std::uint64_t a, const std::uint64_t b, std::uint64_t& high) noexcept -> std::uint64_t
{
    return int128::detail::umul(a, b, high);
}

constexpr auto umul_wide(const std::uint32_t a, const std::uint32_t b, std::uint64_t& high) noexcept -> std::uint64_t
{
    high = 0U;
    return static_cast<std::uint64_t>(a) * b;
}

// MSVC 14.1 warns of unary minus being applied to unsigned type from numeric_limits::min
// 14.2 and on get it right
#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable: 4146)
#endif

// Rounds x to an integer of type Int in the mode round, and saturates out of range
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE T, BOOST_DECIMAL_INTEGRAL Int>
constexpr auto round_integral_int_impl(const T num, const rounding_mode round) noexcept -> Int
{
    using sig_type = typename T::significand_type;
    using uint_type = std::make_unsigned_t<Int>;
    using work_type = std::conditional_t<(sizeof(sig_type) > sizeof(uint_type)), sig_type, uint_type>;
    static_assert(sizeof(uint_type) <= sizeof(std::uint64_t), "The product below has 64 bits for the result");

    #ifndef BOOST_DECIMAL_FAST_MATH
    if (!isfinite(num))
    {
        // Implementation defined what to return here
        return std::numeric_limits<Int>::min();
    }
    #endif

    const auto c {decompose(num)};
    if (c.sig == 0U)
    {
        return 0;
    }

    // Largest magnitude of the result: max for x > 0, and -min for x < 0
    const auto limit {static_cast<work_type>(static_cast<uint_type>((std::numeric_limits<Int>::max)()) + (c.sign ? 1U : 0U))};
    work_type mag {};
    if (c.exp >= 0)
    {
        // 10^exp or the significand above the limit overflows: saturate. Else both factors fit in factor_type.
        if (c.exp > std::numeric_limits<uint_type>::digits10 || c.sig > limit)
        {
            return round_integral_saturate<Int>(c.sign);
        }
        using factor_type = std::conditional_t<(sizeof(uint_type) > sizeof(std::uint32_t)), std::uint64_t, std::uint32_t>;
        std::uint64_t high {};
        const auto low {umul_wide(static_cast<factor_type>(c.sig), pow10(static_cast<factor_type>(c.exp)), high)};
        if (high != 0U || low > limit)
        {
            return round_integral_saturate<Int>(c.sign);
        }
        mag = static_cast<work_type>(low);
    }
    else
    {
        // |x| >= 10^(digits10 + 1) is out of range: saturate before the division
        constexpr int int_digits {std::numeric_limits<uint_type>::digits10 + 1};
        const auto shift {-static_cast<int>(c.exp)};
        if (shift < precision_v<T> - int_digits && c.sig >= pow10(static_cast<sig_type>(shift + int_digits)))
        {
            return round_integral_saturate<Int>(c.sign);
        }
        mag = static_cast<work_type>(round_integral_sig<T>(c, round));
    }

    if (mag > limit)
    {
        return round_integral_saturate<Int>(c.sign);
    }

    const auto umag {static_cast<uint_type>(mag)};
    return static_cast<Int>(c.sign ? static_cast<uint_type>(0U - umag) : umag);
}

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_IMPL_ROUND_INTEGRAL_HPP
