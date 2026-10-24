// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DPD_CONVERSION_HPP
#define BOOST_DECIMAL_DPD_CONVERSION_HPP

#include <boost/decimal/decimal32_t.hpp>
#include <boost/decimal/decimal64_t.hpp>
#include <boost/decimal/decimal128_t.hpp>
#include <boost/decimal/decimal_fast32_t.hpp>
#include <boost/decimal/decimal_fast64_t.hpp>
#include <boost/decimal/decimal_fast128_t.hpp>
#include <boost/decimal/bid_conversion.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/uint128_t.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <limits>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The declet of the three digits of n < 1000 (IEEE 754 Table 3.4)
BOOST_DECIMAL_CUDA_CONSTEXPR auto encode_dpd(const std::uint32_t n) noexcept -> std::uint32_t
{
    const auto d1 {n / 100U};
    const auto d2 {n / 10U % 10U};
    const auto d3 {n % 10U};

    // The low bits of the digits that are 8 or 9
    const auto d {(d1 & 1U) << 7U};
    const auto h {(d2 & 1U) << 4U};
    const auto m {d3 & 1U};

    switch (((d1 >> 3U) << 2U) | ((d2 >> 3U) << 1U) | (d3 >> 3U))
    {
        case 0U:
            return (d1 << 7U) | (d2 << 4U) | d3;
        case 1U:
            return (d1 << 7U) | (d2 << 4U) | 0b1000U | m;
        case 2U:
            return (d1 << 7U) | ((d3 & 6U) << 4U) | h | 0b1010U | m;
        case 3U:
            return (d1 << 7U) | 0b1000000U | h | 0b1110U | m;
        case 4U:
            return ((d3 & 6U) << 7U) | d | (d2 << 4U) | 0b1100U | m;
        case 5U:
            return ((d2 & 6U) << 7U) | d | 0b0100000U | h | 0b1110U | m;
        case 6U:
            return ((d3 & 6U) << 7U) | d | h | 0b1110U | m;
        default:
            return d | 0b1100000U | h | 0b1110U | m;
    }
}

// The number n < 1000 of a declet (IEEE 754 Table 3.3), which also reads the noncanonical declets
BOOST_DECIMAL_CUDA_CONSTEXPR auto decode_dpd(const std::uint32_t declet) noexcept -> std::uint32_t
{
    // The bit pairs pq, st and wx as even numbers, and the low bits r, u and y
    const auto pq {(declet >> 7U) & 6U};
    const auto st {(declet >> 4U) & 6U};
    const auto r {(declet >> 7U) & 1U};
    const auto u {(declet >> 4U) & 1U};
    const auto y {declet & 1U};

    std::uint32_t d1 {pq | r};
    std::uint32_t d2 {st | u};
    std::uint32_t d3 {(declet & 6U) | y};
    if ((declet & 8U) != 0U)
    {
        switch ((declet >> 1U) & 3U)
        {
            case 0U:
                d3 = 8U | y;
                break;
            case 1U:
                d2 = 8U | u;
                d3 = st | y;
                break;
            case 2U:
                d1 = 8U | r;
                d3 = pq | y;
                break;
            default:
                switch (st >> 1U)
                {
                    case 0U:
                        d1 = 8U | r;
                        d2 = 8U | u;
                        d3 = pq | y;
                        break;
                    case 1U:
                        d1 = 8U | r;
                        d2 = pq | u;
                        d3 = 8U | y;
                        break;
                    case 2U:
                        d2 = 8U | u;
                        d3 = 8U | y;
                        break;
                    default:
                        d1 = 8U | r;
                        d2 = 8U | u;
                        d3 = 8U | y;
                        break;
                }
                break;
        }
    }

    return d1 * 100U + d2 * 10U + d3;
}

// The declets of 0 to 999 and the numbers of the 1024 declets, made at compile time
struct dpd_declets
{
    std::uint16_t value[1000];
};

struct dpd_numbers
{
    std::uint16_t value[1024];
};

constexpr auto make_dpd_declets() noexcept -> dpd_declets
{
    dpd_declets t {};
    for (std::uint32_t n {}; n < 1000U; ++n)
    {
        t.value[n] = static_cast<std::uint16_t>(encode_dpd(n));
    }
    return t;
}

constexpr auto make_dpd_numbers() noexcept -> dpd_numbers
{
    dpd_numbers t {};
    for (std::uint32_t n {}; n < 1024U; ++n)
    {
        t.value[n] = static_cast<std::uint16_t>(decode_dpd(n));
    }
    return t;
}

// Members of a class template, so that a C++14 program has one copy of each table that it uses
template <bool>
struct dpd_table
{
    static constexpr dpd_declets declets {make_dpd_declets()};
    static constexpr dpd_numbers numbers {make_dpd_numbers()};
};

#if !defined(__cpp_inline_variables) || __cpp_inline_variables < 201606L

template <bool b> constexpr dpd_declets dpd_table<b>::declets;
template <bool b> constexpr dpd_numbers dpd_table<b>::numbers;

#endif

// Device code computes each declet, and host code from every compiler reads the same tables

BOOST_DECIMAL_CUDA_CONSTEXPR auto declet_of(const std::uint32_t n) noexcept -> std::uint32_t
{
    #ifdef __CUDA_ARCH__
    return encode_dpd(n);
    #else
    return dpd_table<true>::declets.value[n];
    #endif
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto number_of(const std::uint32_t declet) noexcept -> std::uint32_t
{
    #ifdef __CUDA_ARCH__
    return decode_dpd(declet);
    #else
    return dpd_table<true>::numbers.value[declet];
    #endif
}

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int dpd_declet_bits {10};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint32_t dpd_declet_mask {UINT32_C(0x3FF)};

// The K declets of the low 3K digits of n
template <int K, typename T>
BOOST_DECIMAL_CUDA_CONSTEXPR auto encode_declets(T n) noexcept -> T
{
    T bits {};
    for (int i {}; i < K; ++i)
    {
        bits |= static_cast<T>(declet_of(static_cast<std::uint32_t>(n % 1000U))) << (dpd_declet_bits * i);
        n /= 1000U;
    }
    return bits;
}

// The number of K declets
template <int K, typename T>
BOOST_DECIMAL_CUDA_CONSTEXPR auto decode_declets(const T bits) noexcept -> T
{
    T n {};
    for (int i {K}; i-- > 0;)
    {
        n = n * 1000U + number_of(static_cast<std::uint32_t>(bits >> (dpd_declet_bits * i)) & dpd_declet_mask);
    }
    return n;
}

// decimal128_t has 11 declets: the low 6 hold 18 digits in 60 bits, and the high 5 hold 15 digits
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t dpd_ten_15 {UINT64_C(1000000000000000)};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t dpd_ten_18 {UINT64_C(1000000000000000000)};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned dpd_low_bits {6U * 10U};

struct dpd_declets_128
{
    int128::uint128_t declets;
    std::uint32_t lead;
};

// The 11 declets of the low 33 digits of n < 10^34, and the leading digit
BOOST_DECIMAL_CUDA_CONSTEXPR auto encode_declets_128(const int128::uint128_t& n) noexcept -> dpd_declets_128
{
    // Both products are modulo 2^64, and the true remainder is below 10^18
    const auto top {static_cast<std::uint64_t>(n / dpd_ten_18)};
    const auto low {n.low - top * dpd_ten_18};
    return {(static_cast<int128::uint128_t>(encode_declets<5>(top % dpd_ten_15)) << dpd_low_bits) | encode_declets<6>(low),
            static_cast<std::uint32_t>(top / dpd_ten_15)};
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto decode_declets_128(const int128::uint128_t& bits) noexcept -> int128::uint128_t
{
    const auto high {decode_declets<5>(static_cast<std::uint64_t>(bits >> dpd_low_bits))};
    const auto low {decode_declets<6>(bits.low & ((UINT64_C(1) << dpd_low_bits) - 1U))};
    return static_cast<int128::uint128_t>(high) * dpd_ten_18 + low;
}

// The combination field of the leading digit d0 and the top two bits of the biased exponent (IEEE 754 3.5.2.c.1)
BOOST_DECIMAL_CUDA_CONSTEXPR auto dpd_combination(const std::uint32_t d0, const std::uint32_t lead) noexcept -> std::uint32_t
{
    return d0 >= 8U ? 0b11000U | (lead << 1U) | (d0 & 1U) : (lead << 3U) | d0;
}

// The leading digit and the top two bits of the biased exponent of a finite combination field
BOOST_DECIMAL_CUDA_CONSTEXPR auto dpd_leading_digit(const std::uint32_t combination) noexcept -> std::uint32_t
{
    return combination >= 0b11000U ? 8U | (combination & 1U) : combination & 7U;
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto dpd_leading_exponent(const std::uint32_t combination) noexcept -> std::uint32_t
{
    return combination >= 0b11000U ? (combination >> 1U) & 3U : combination >> 3U;
}

// The bits that an infinity keeps are the sign and the class. A NaN also keeps the signaling bit.
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint32_t d32_dpd_inf_bits {d32_sign_mask | d32_inf_mask};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint32_t d32_dpd_nan_bits {d32_sign_mask | d32_snan_mask};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t d64_dpd_inf_bits {d64_sign_mask | d64_inf_mask};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t d64_dpd_nan_bits {d64_sign_mask | d64_snan_mask};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t d128_dpd_inf_bits {d128_sign_mask | d128_inf_mask.high};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t d128_dpd_nan_bits {d128_sign_mask | d128_snan_mask.high};

// The widths of the exponent continuation and the trailing significand, which is in declets
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d32_dpd_exp_bits {6U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d32_dpd_trailing_bits {2U * 10U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d64_dpd_exp_bits {8U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d64_dpd_trailing_bits {5U * 10U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d128_dpd_exp_bits {12U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d128_dpd_trailing_high_bits {11U * 10U - 64U};

// The 5-bit combination field is below the sign bit
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d32_dpd_combination_shift {d32_dpd_exp_bits + d32_dpd_trailing_bits};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d64_dpd_combination_shift {d64_dpd_exp_bits + d64_dpd_trailing_bits};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned d128_dpd_combination_shift {d128_dpd_exp_bits + d128_dpd_trailing_high_bits};

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint32_t d32_dpd_trailing_mask {(UINT32_C(1) << d32_dpd_trailing_bits) - 1U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t d64_dpd_trailing_mask {(UINT64_C(1) << d64_dpd_trailing_bits) - 1U};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t d128_dpd_trailing_mask {(UINT64_C(1) << d128_dpd_trailing_high_bits) - 1U, UINT64_MAX};

// 10^(p-1), which is the weight of the leading digit
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint32_t d32_dpd_lead_weight {UINT32_C(1000000)};
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t d128_dpd_lead_weight {static_cast<int128::uint128_t>(dpd_ten_15) * dpd_ten_18};

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename DecimalType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd_d32(const DecimalType val) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_decimal_floating_point_v, DecimalType, std::uint32_t)
{
    static_assert(std::is_same<DecimalType, decimal32_t>::value ||
                  std::is_same<DecimalType, decimal_fast32_t>::value, "The input must be a 32-bit decimal type");

    if (!isfinite(val))
    {
        const auto bid {to_bid(val)};
        if ((bid & detail::d32_nan_mask) != detail::d32_nan_mask)
        {
            return bid & detail::d32_dpd_inf_bits;
        }
        // The payload of a NaN goes in declets, and a noncanonical one is 0
        const auto payload {bid & detail::d32_dpd_trailing_mask};
        return (bid & detail::d32_dpd_nan_bits) |
               (payload >= detail::d32_dpd_lead_weight ? 0U : detail::encode_declets<2>(payload));
    }

    const auto significand {static_cast<std::uint32_t>(val.full_significand())};
    BOOST_DECIMAL_ASSERT(significand < 10U * detail::d32_dpd_lead_weight);
    const auto exp {static_cast<std::uint32_t>(val.unbiased_exponent())};
    const auto combination {detail::dpd_combination(significand / detail::d32_dpd_lead_weight, exp >> detail::d32_dpd_exp_bits)};

    return (val.isneg() ? detail::d32_sign_mask : 0U) | (combination << detail::d32_dpd_combination_shift) |
           ((exp & ((1U << detail::d32_dpd_exp_bits) - 1U)) << detail::d32_dpd_trailing_bits) |
           detail::encode_declets<2>(significand % detail::d32_dpd_lead_weight);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal_fast32_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd_d32(const std::uint32_t dpd) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    static_assert(std::is_same<DecimalType, decimal32_t>::value || std::is_same<DecimalType, decimal_fast32_t>::value,
                  "Target decimal type must be 32-bits");

    // A NaN keeps its payload, and the other bits of an infinity or a NaN are as in BID
    if ((dpd & detail::d32_inf_mask) == detail::d32_inf_mask)
    {
        return from_bid<DecimalType>((dpd & detail::d32_nan_mask) == detail::d32_nan_mask ?
                                     (dpd & detail::d32_dpd_nan_bits) | detail::decode_declets<2>(dpd & detail::d32_dpd_trailing_mask) :
                                     dpd & detail::d32_dpd_inf_bits);
    }

    const auto combination {dpd >> detail::d32_dpd_combination_shift & 0x1FU};
    const auto exp {(detail::dpd_leading_exponent(combination) << detail::d32_dpd_exp_bits) |
                    ((dpd >> detail::d32_dpd_trailing_bits) & ((1U << detail::d32_dpd_exp_bits) - 1U))};
    const auto significand {detail::dpd_leading_digit(combination) * detail::d32_dpd_lead_weight +
                            detail::decode_declets<2>(dpd & detail::d32_dpd_trailing_mask)};

    return DecimalType{significand, static_cast<std::int32_t>(exp) - detail::bias_v<DecimalType>, (dpd & detail::d32_sign_mask) != 0U};
}

BOOST_DECIMAL_EXPORT template <typename DecimalType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd_d64(const DecimalType val) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_decimal_floating_point_v, DecimalType, std::uint64_t)
{
    static_assert(std::is_same<DecimalType, decimal64_t>::value ||
                  std::is_same<DecimalType, decimal_fast64_t>::value, "The input must be a 64-bit decimal type");

    if (!isfinite(val))
    {
        const auto bid {to_bid(val)};
        if ((bid & detail::d64_nan_mask) != detail::d64_nan_mask)
        {
            return bid & detail::d64_dpd_inf_bits;
        }
        // The payload of a NaN goes in declets, and a noncanonical one is 0
        const auto payload {bid & detail::d64_dpd_trailing_mask};
        return (bid & detail::d64_dpd_nan_bits) | (payload >= detail::dpd_ten_15 ? 0U : detail::encode_declets<5>(payload));
    }

    const auto significand {static_cast<std::uint64_t>(val.full_significand())};
    BOOST_DECIMAL_ASSERT(significand < 10U * detail::dpd_ten_15);
    const auto exp {static_cast<std::uint64_t>(val.unbiased_exponent())};
    const auto combination {detail::dpd_combination(static_cast<std::uint32_t>(significand / detail::dpd_ten_15),
                                                    static_cast<std::uint32_t>(exp >> detail::d64_dpd_exp_bits))};

    return (val.isneg() ? detail::d64_sign_mask : 0U) | (static_cast<std::uint64_t>(combination) << detail::d64_dpd_combination_shift) |
           ((exp & ((UINT64_C(1) << detail::d64_dpd_exp_bits) - 1U)) << detail::d64_dpd_trailing_bits) |
           detail::encode_declets<5>(significand % detail::dpd_ten_15);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal_fast64_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd_d64(const std::uint64_t dpd) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    static_assert(std::is_same<DecimalType, decimal64_t>::value || std::is_same<DecimalType, decimal_fast64_t>::value,
                  "Target decimal type must be 64-bits");

    // A NaN keeps its payload, and the other bits of an infinity or a NaN are as in BID
    if ((dpd & detail::d64_inf_mask) == detail::d64_inf_mask)
    {
        return from_bid<DecimalType>((dpd & detail::d64_nan_mask) == detail::d64_nan_mask ?
                                     (dpd & detail::d64_dpd_nan_bits) | detail::decode_declets<5>(dpd & detail::d64_dpd_trailing_mask) :
                                     dpd & detail::d64_dpd_inf_bits);
    }

    const auto combination {static_cast<std::uint32_t>(dpd >> detail::d64_dpd_combination_shift) & 0x1FU};
    const auto exp {(static_cast<std::uint64_t>(detail::dpd_leading_exponent(combination)) << detail::d64_dpd_exp_bits) |
                    ((dpd >> detail::d64_dpd_trailing_bits) & ((UINT64_C(1) << detail::d64_dpd_exp_bits) - 1U))};
    const auto significand {detail::dpd_leading_digit(combination) * detail::dpd_ten_15 +
                            detail::decode_declets<5>(dpd & detail::d64_dpd_trailing_mask)};

    return DecimalType{significand, static_cast<std::int32_t>(exp) - detail::bias_v<DecimalType>, (dpd & detail::d64_sign_mask) != 0U};
}

BOOST_DECIMAL_EXPORT template <typename DecimalType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd_d128(const DecimalType val) noexcept
    BOOST_DECIMAL_REQUIRES_RETURN(detail::is_decimal_floating_point_v, DecimalType, boost::int128::uint128_t)
{
    static_assert(std::is_same<DecimalType, decimal128_t>::value ||
                  std::is_same<DecimalType, decimal_fast128_t>::value, "The input must be a 128-bit decimal type");

    if (!isfinite(val))
    {
        const auto bid {to_bid(val)};
        if ((bid & detail::d128_nan_mask) != detail::d128_nan_mask)
        {
            return int128::uint128_t {bid.high & detail::d128_dpd_inf_bits, UINT64_C(0)};
        }
        // The payload of a NaN goes in declets, and a noncanonical one is 0
        const auto payload {bid & detail::d128_dpd_trailing_mask};
        int128::uint128_t dpd {payload >= detail::d128_dpd_lead_weight ? int128::uint128_t{} :
                               detail::encode_declets_128(payload).declets};
        dpd.high |= bid.high & detail::d128_dpd_nan_bits;
        return dpd;
    }

    const int128::uint128_t significand {val.full_significand()};
    BOOST_DECIMAL_ASSERT(significand < detail::d128_dpd_lead_weight * 10U);
    const auto exp {static_cast<std::uint64_t>(val.unbiased_exponent())};

    auto dpd {detail::encode_declets_128(significand)};
    const auto combination {detail::dpd_combination(dpd.lead, static_cast<std::uint32_t>(exp >> detail::d128_dpd_exp_bits))};
    dpd.declets.high |= (val.isneg() ? detail::d128_sign_mask : 0U) |
                        (static_cast<std::uint64_t>(combination) << detail::d128_dpd_combination_shift) |
                        ((exp & ((UINT64_C(1) << detail::d128_dpd_exp_bits) - 1U)) << detail::d128_dpd_trailing_high_bits);
    return dpd.declets;
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal_fast128_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd_d128(const int128::uint128_t dpd) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    static_assert(std::is_same<DecimalType, decimal128_t>::value || std::is_same<DecimalType, decimal_fast128_t>::value,
                  "Target decimal type must be 128-bits");

    // A NaN keeps its payload, and the other bits of an infinity or a NaN are as in BID
    if ((dpd & detail::d128_inf_mask) == detail::d128_inf_mask)
    {
        int128::uint128_t bid {};
        if ((dpd & detail::d128_nan_mask) == detail::d128_nan_mask)
        {
            bid = detail::decode_declets_128(dpd & detail::d128_dpd_trailing_mask);
            bid.high |= dpd.high & detail::d128_dpd_nan_bits;
        }
        else
        {
            bid.high = dpd.high & detail::d128_dpd_inf_bits;
        }
        return from_bid<DecimalType>(bid);
    }

    const auto combination {static_cast<std::uint32_t>(dpd.high >> detail::d128_dpd_combination_shift) & 0x1FU};
    const auto exp {(static_cast<std::uint64_t>(detail::dpd_leading_exponent(combination)) << detail::d128_dpd_exp_bits) |
                    ((dpd.high >> detail::d128_dpd_trailing_high_bits) & ((UINT64_C(1) << detail::d128_dpd_exp_bits) - 1U))};
    const auto significand {detail::dpd_leading_digit(combination) * detail::d128_dpd_lead_weight +
                            detail::decode_declets_128(dpd & detail::d128_dpd_trailing_mask)};

    return DecimalType{significand, static_cast<std::int32_t>(exp) - detail::bias_v<DecimalType>, (dpd.high & detail::d128_sign_mask) != 0U};
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal32_t val) noexcept -> std::uint32_t
{
    return to_dpd_d32(val);
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal_fast32_t val) noexcept -> std::uint32_t
{
    return to_dpd_d32(val);
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal64_t val) noexcept -> std::uint64_t
{
    return to_dpd_d64(val);
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal_fast64_t val) noexcept -> std::uint64_t
{
    return to_dpd_d64(val);
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal128_t val) noexcept -> int128::uint128_t
{
    return to_dpd_d128(val);
}

BOOST_DECIMAL_EXPORT BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const decimal_fast128_t& val) noexcept -> int128::uint128_t
{
    return to_dpd_d128(val);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_dpd(const DecimalType val) noexcept
{
    static_assert(detail::is_decimal_floating_point_v<DecimalType>, "Must be a decimal floating point type.");
    return to_dpd(val);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal32_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd(const std::uint32_t bits) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    return from_dpd_d32<DecimalType>(bits);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal64_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd(const std::uint64_t bits) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    return from_dpd_d64<DecimalType>(bits);
}

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal128_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd(const int128::uint128_t bits) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    return from_dpd_d128<DecimalType>(bits);
}

#ifdef BOOST_DECIMAL_HAS_INT128

BOOST_DECIMAL_EXPORT template <typename DecimalType = decimal128_t>
BOOST_DECIMAL_CUDA_CONSTEXPR auto from_dpd(const detail::builtin_uint128_t bits) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, DecimalType)
{
    const int128::uint128_t converted_bits {bits};
    return from_dpd_d128<DecimalType>(converted_bits);
}

#endif // BOOST_DECIMAL_HAS_INT128

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DPD_CONVERSION_HPP
