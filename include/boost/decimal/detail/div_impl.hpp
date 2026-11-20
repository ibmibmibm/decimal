// Copyright 2023 - 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_DIV_IMPL_HPP
#define BOOST_DECIMAL_DETAIL_DIV_IMPL_HPP

#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/add_impl.hpp>
#include <boost/decimal/detail/select_bits.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>
#include <boost/decimal/detail/normalize.hpp>
#include <boost/decimal/detail/power_tables.hpp>
#include <boost/decimal/detail/remove_trailing_zeros.hpp>
#include <boost/decimal/detail/components.hpp>
#include <boost/decimal/detail/u256.hpp>
#include "int128.hpp"

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4127) // Conditional expression is constant (pre-C++17 if-constexpr fallback)
#endif

namespace impl {

// Division finalization with inline round-half-to-even. The finalizers
// round in the mode that the caller gives, with this rule for
// fe_dec_to_nearest. The driver computes
// q = (lhs.sig * 10^k) / rhs.sig and r = (lhs.sig * 10^k) % rhs.sig where
// both operand significands are at full precision p, and k is p - 1 when
// lhs.sig >= rhs.sig, else p. Thus q has exactly p digits, and the natural
// remainder r encodes the exact sticky bit needed for correct RHTE.
//
// Round up iff 2r > divisor, or 2r == divisor and q is odd. With
// half_div = divisor - r, this is r + (q & 1) > half_div, without overflow.
//
// The post-rounding carry q == 10^p shifts down to 10^(p-1) and bumps the
// exponent by one digit.

// IEEE 754-2019 5.4.1: an exact quotient has the exponent nearest to the difference of the
// exponents of the operands. Zeros at the end of q go in steps of FirstStep, ..., 2 and 1 digits,
// and q has fewer than 2 * FirstStep of them.
template <int FirstStep, typename SigType, typename ExpType>
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto exact_quotient_exponent(SigType& q, ExpType& exp, const ExpType ideal_exp) noexcept -> void
{
    for (int step {FirstStep}; step != 0; step /= 2)
    {
        if (exp + step <= ideal_exp)
        {
            const auto power {pow10(static_cast<SigType>(static_cast<std::uint64_t>(step)))};
            if (q % power == 0U)
            {
                q /= power;
                exp += step;
            }
        }
    }
}

#if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)

// Inverses of 5^k modulo 2^128
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int128::uint128_t inverse_pow5_u128[34] = {
    {UINT64_C(0x0000000000000000), UINT64_C(0x0000000000000001)},
    {UINT64_C(0xCCCCCCCCCCCCCCCC), UINT64_C(0xCCCCCCCCCCCCCCCD)},
    {UINT64_C(0x28F5C28F5C28F5C2), UINT64_C(0x8F5C28F5C28F5C29)},
    {UINT64_C(0x6E978D4FDF3B645A), UINT64_C(0x1CAC083126E978D5)},
    {UINT64_C(0x495182A9930BE0DE), UINT64_C(0xD288CE703AFB7E91)},
    {UINT64_C(0xDB76B3BB83CF2CF9), UINT64_C(0x5D4E8FB00BCBE61D)},
    {UINT64_C(0xC57E23F24D8FD5CB), UINT64_C(0x790FB65668C26139)},
    {UINT64_C(0xC1193A63A91CC45B), UINT64_C(0xE5032477AE8D46A5)},
    {UINT64_C(0xF36B7213EE9F5A78), UINT64_C(0xC767074B22E90E21)},
    {UINT64_C(0x97157D372FB9787E), UINT64_C(0x8E47CE423A2E9C6D)},
    {UINT64_C(0x516AB2A4A3251819), UINT64_C(0x4FA7F60D3ED61F49)},
    {UINT64_C(0x76AEF08753D43805), UINT64_C(0x0FEE64690C913975)},
    {UINT64_C(0xB156301B10C40B34), UINT64_C(0x3662E0E1CF503EB1)},
    {UINT64_C(0x2377A3389CF4023D), UINT64_C(0xA47A2CF9F6433FBD)},
    {UINT64_C(0x0717ED71B8FD9A0C), UINT64_C(0x54186F653140A659)},
    {UINT64_C(0xCE37FC49F1CC5202), UINT64_C(0x7738164770402145)},
    {UINT64_C(0xF60B3275305C1066), UINT64_C(0xE4A4D1417CD9A041)},
    {UINT64_C(0x6468A3B109AC0347), UINT64_C(0xC75429D9E5C5200D)},
    {UINT64_C(0xE0E1BA569B88CD74), UINT64_C(0xC1773B91FAC10669)},
    {UINT64_C(0x93605877B8B4F5E4), UINT64_C(0x26B172506559CE15)},
    {UINT64_C(0x83E011B18B576460), UINT64_C(0xD489E3A9ADDEC2D1)},
    {UINT64_C(0x4D9336BD1BDE4746), UINT64_C(0x90E860BB892C8D5D)},
    {UINT64_C(0xDC50A48C38C60E41), UINT64_C(0x502E79BF1B6F4F79)},
    {UINT64_C(0x9276874F3E8E02D9), UINT64_C(0xDCD618596BE30FE5)},
    {UINT64_C(0xEA17B4A972E933C5), UINT64_C(0x2C2AD1AB7BFA3661)},
    {UINT64_C(0xFB9E575516FB70C1), UINT64_C(0x08D55D224BFED7AD)},
    {UINT64_C(0x658611776AFF168D), UINT64_C(0x01C445D3A8CC9189)},
    {UINT64_C(0xE11AD04B156637B5), UINT64_C(0xCD27412A54F5B6B5)},
    {UINT64_C(0xF9D229A89DE13E57), UINT64_C(0x8F6E403BAA978AF1)},
    {UINT64_C(0x31F6D521B92D0C77), UINT64_C(0xE97C733F221E4EFD)},
    {UINT64_C(0xA397C439F1D5CF4B), UINT64_C(0x2EB27D7306D2DC99)},
    {UINT64_C(0xED84C0D863912975), UINT64_C(0x6FBD4C4A34909285)},
    {UINT64_C(0x62B42691AD836EB1), UINT64_C(0x16590F420A835081)},
    {UINT64_C(0x46F0D483891A4956), UINT64_C(0x9E11CFDA021A434D)}
};

#endif

// A division of two words is slow, so this overload removes all the zeros at one time, and puts
// back those that are past the preferred exponent
template <typename ExpType>
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto exact_quotient_exponent(int128::uint128_t& q, ExpType& exp, const ExpType ideal_exp) noexcept -> void
{
    if (exp < ideal_exp)
    {
        #if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)
        // Most quotients reach the preferred exponent. For a multiple of 10^gap, (q >> gap) times
        // the inverse of 5^gap is q / 10^gap, and any other q gives a product of 10^(34 - gap) or more.
        const auto gap {ideal_exp - exp < 33 ? static_cast<unsigned>(ideal_exp - exp) : 33U};
        const auto m {(q >> gap) * inverse_pow5_u128[gap]};
        if ((q.low & ((UINT64_C(1) << gap) - 1U)) == 0U && m < pow10(int128::uint128_t{static_cast<std::uint64_t>(34U - gap)}))
        {
            q = m;
            exp += static_cast<ExpType>(gap);
            return;
        }
        #endif

        const auto zeros {remove_trailing_zeros(q)};
        const auto removed {static_cast<ExpType>(zeros.number_of_removed_zeros)};
        const auto steps {removed < ideal_exp - exp ? removed : ideal_exp - exp};
        q = zeros.trimmed_number * pow10(int128::uint128_t{static_cast<std::uint64_t>(removed - steps)});
        exp += steps;
    }
}

// d32/fast32 finalizer. Dividend = lhs.sig * 10^7 fits in uint64; quotient is
// at most 10^7 - 1 which fits in uint32; remainder is bounded by the divisor
// which fits in uint32 (rhs.sig < 10^7).
template <typename ReturnType, typename ExpType>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto div_finalize_u64(
    std::uint64_t q, std::uint64_t r, std::uint32_t divisor,
    ExpType result_exp, ExpType ideal_exp, bool sign, const rounding_mode round) noexcept -> ReturnType
{
    constexpr auto ten_to_p {pow10(static_cast<std::uint64_t>(detail::precision_v<ReturnType>))};
    constexpr auto ten_to_p_minus_1 {pow10(static_cast<std::uint64_t>(detail::precision_v<ReturnType> - 1))};

    // An exact quotient takes the exponent nearest to the preferred one
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(r == 0U))
    {
        exact_quotient_exponent<4>(q, result_exp, ideal_exp);
        return pack_in_range<ReturnType>(static_cast<typename ReturnType::significand_type>(q), result_exp, sign);
    }

    // Below the subnormal exponent the constructor rounds once, with the remainder as a sticky digit.
    // Rounding to the precision here first would round the result twice.
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(result_exp < detail::etiny_v<ReturnType>))
    {
        return pack_in_range<ReturnType>(static_cast<typename ReturnType::significand_type>(q * 10U + static_cast<std::uint64_t>(r != 0U)), result_exp - 1, sign);
    }

    int extra {0};

    const auto half_div {static_cast<std::uint64_t>(divisor) - r};
    const bool round_up {round == rounding_mode::fe_dec_to_nearest ?
                         r + (q & UINT64_C(1)) > half_div :
                         detail::steps_up_in_mode(round, sign, r != 0U, r >= half_div)};
    q += static_cast<unsigned>(round_up);
    if (BOOST_DECIMAL_UNLIKELY(q == ten_to_p))
    {
        q = ten_to_p_minus_1;
        ++extra;
    }

    using sig_type = typename ReturnType::significand_type;
    return pack_in_range<ReturnType>(static_cast<sig_type>(q),
                                     result_exp + static_cast<ExpType>(extra),
                                     sign);
}

// Out-of-line helper for the directed modes, without FORCE_INLINE: an inlined
// second copy of the finalizer makes the operators slower in the nearest mode.
template <typename ReturnType, typename ExpType>
BOOST_DECIMAL_CUDA_CONSTEXPR auto div_finalize_u64_out_of_line(
    std::uint64_t q, std::uint64_t r, std::uint32_t divisor,
    ExpType result_exp, ExpType ideal_exp, bool sign, const rounding_mode round) noexcept -> ReturnType
{
    return div_finalize_u64<ReturnType>(q, r, divisor, result_exp, ideal_exp, sign, round);
}

// d64/fast64 finalizer. Quotient and remainder both fit in uint64 because
// quotient is at most 10^16 - 1 < 2^54 and remainder is bounded by the
// uint64 divisor.
template <typename ReturnType, typename ExpType>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto div_finalize_u128(
    std::uint64_t q, std::uint64_t r, std::uint64_t divisor,
    ExpType result_exp, ExpType ideal_exp, bool sign, const rounding_mode round) noexcept -> ReturnType
{
    constexpr auto ten_to_p {pow10(static_cast<std::uint64_t>(detail::precision_v<ReturnType>))};
    constexpr auto ten_to_p_minus_1 {pow10(static_cast<std::uint64_t>(detail::precision_v<ReturnType> - 1))};

    // An exact quotient takes the exponent nearest to the preferred one
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(r == 0U))
    {
        exact_quotient_exponent<8>(q, result_exp, ideal_exp);
        return pack_in_range<ReturnType>(static_cast<typename ReturnType::significand_type>(q), result_exp, sign);
    }

    // Below the subnormal exponent the constructor rounds once, with the remainder as a sticky digit.
    // Rounding to the precision here first would round the result twice.
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(result_exp < detail::etiny_v<ReturnType>))
    {
        return pack_in_range<ReturnType>(static_cast<typename ReturnType::significand_type>(q * 10U + static_cast<std::uint64_t>(r != 0U)), result_exp - 1, sign);
    }

    int extra {0};

    const auto half_div {divisor - r};
    const bool round_up {round == rounding_mode::fe_dec_to_nearest ?
                         r + (q & UINT64_C(1)) > half_div :
                         detail::steps_up_in_mode(round, sign, r != 0U, r >= half_div)};
    q += static_cast<unsigned>(round_up);
    if (BOOST_DECIMAL_UNLIKELY(q == ten_to_p))
    {
        q = ten_to_p_minus_1;
        ++extra;
    }

    using sig_type = typename ReturnType::significand_type;
    return pack_in_range<ReturnType>(static_cast<sig_type>(q),
                                     result_exp + static_cast<ExpType>(extra),
                                     sign);
}

// d128/fast128 finalizer. Quotient and remainder both fit in uint128 because
// quotient is at most 10^34 - 1 (< 2^113) and remainder is bounded by the
// uint128 divisor (< 10^34).
template <typename ReturnType, typename ExpType>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto div_finalize_u256(
    int128::uint128_t q, int128::uint128_t r, int128::uint128_t divisor,
    ExpType result_exp, ExpType ideal_exp, bool sign, const rounding_mode round) noexcept -> ReturnType
{
    constexpr auto ten_to_p {pow10(int128::uint128_t{static_cast<std::uint64_t>(detail::precision_v<ReturnType>)})};
    constexpr auto ten_to_p_minus_1 {pow10(int128::uint128_t{static_cast<std::uint64_t>(detail::precision_v<ReturnType> - 1)})};

    // An exact quotient takes the exponent nearest to the preferred one
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(r == 0U))
    {
        exact_quotient_exponent(q, result_exp, ideal_exp);
        return pack_in_range<ReturnType>(q, result_exp, sign);
    }

    // Below the subnormal exponent the constructor rounds once, with the remainder as a sticky digit.
    // Rounding to the precision here first would round the result twice.
    if (!detail::is_fast_type_v<ReturnType> && BOOST_DECIMAL_UNLIKELY(result_exp < detail::etiny_v<ReturnType>))
    {
        return pack_in_range<ReturnType>(static_cast<typename ReturnType::significand_type>(q * 10U + static_cast<std::uint64_t>(r != int128::uint128_t{0U})), result_exp - 1, sign);
    }

    int extra {0};

    const auto half_div {divisor - r};
    const bool round_up {round == rounding_mode::fe_dec_to_nearest ?
                         r + (q.low & UINT64_C(1)) > half_div :
                         detail::steps_up_in_mode(round, sign, r != 0U, r >= half_div)};
    q += static_cast<unsigned>(round_up);
    if (BOOST_DECIMAL_UNLIKELY(q == ten_to_p))
    {
        q = ten_to_p_minus_1;
        ++extra;
    }

    return pack_in_range<ReturnType>(q,
                                     result_exp + static_cast<ExpType>(extra),
                                     sign);
}

} // namespace impl

// d32/fast32 division driver. Accepts a decimal type or a components struct
// (both expose to_components()) and rounds in the current mode.
template <typename DecimalType, typename T>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto generic_div_impl(const T& lhs, const T& rhs) noexcept -> DecimalType
{
    auto lhs_c {lhs.to_components()};
    auto rhs_c {rhs.to_components()};
    const auto ideal_exp {static_cast<int>(lhs_c.exp) - static_cast<int>(rhs_c.exp)};

    detail::expand_significand<DecimalType>(lhs_c.sig, lhs_c.exp);
    detail::expand_significand<DecimalType>(rhs_c.sig, rhs_c.exp);

    const bool sign {lhs_c.sign != rhs_c.sign};

    if (BOOST_DECIMAL_UNLIKELY(lhs_c.sig == 0U))
    {
        using sig_type = typename decltype(lhs_c)::significand_type;
        return DecimalType{sig_type{0U}, ideal_exp, sign};
    }

    // Scale by 10^(p-1) when lhs >= rhs, thus the quotient has exactly p digits
    constexpr auto ten_to_p {pow10(static_cast<std::uint64_t>(precision_v<DecimalType>))};
    const bool shorter {lhs_c.sig >= rhs_c.sig};
    const auto big_sig {static_cast<std::uint64_t>(lhs_c.sig) * select_bits(shorter, ten_to_p / 10U, ten_to_p)};
    const auto divisor {static_cast<std::uint64_t>(rhs_c.sig)};
    const auto q {big_sig / divisor};
    const auto r {big_sig - q * divisor};
    const auto res_exp {(lhs_c.exp - static_cast<int>(precision_v<DecimalType>)) - rhs_c.exp + static_cast<int>(shorter)};

    // Only this driver keeps the other modes out of line: inlined, they made
    // the d32 operators slower in the nearest mode, but not d64 or d128.
    // Not const, see current_rounding_mode
    auto round {detail::current_rounding_mode()};
    if (BOOST_DECIMAL_UNLIKELY(round != rounding_mode::fe_dec_to_nearest))
    {
        return impl::div_finalize_u64_out_of_line<DecimalType>(q, r, static_cast<std::uint32_t>(divisor), res_exp, ideal_exp, sign, round);
    }

    return impl::div_finalize_u64<DecimalType>(q, r, static_cast<std::uint32_t>(divisor), res_exp, ideal_exp, sign,
                                               rounding_mode::fe_dec_to_nearest);
}

// d64/fast64 division driver. Same structure as the d32 driver but uses
// uint128 for the pre-scaled dividend (lhs.sig * 10^16 overflows uint64).
// The quotient is at most 10^16 - 1 < 2^54 so it narrows safely to uint64
// for the finalizer.
template <typename DecimalType, typename T>
BOOST_DECIMAL_FORCE_INLINE BOOST_DECIMAL_CUDA_CONSTEXPR auto d64_generic_div_impl(const T& lhs, const T& rhs, const bool sign) noexcept -> DecimalType
{
    using unsigned_int128_type = boost::int128::uint128_t;

    auto lhs_c {lhs.to_components()};
    auto rhs_c {rhs.to_components()};
    const auto ideal_exp {static_cast<int>(lhs_c.exp) - static_cast<int>(rhs_c.exp)};

    detail::expand_significand<DecimalType>(lhs_c.sig, lhs_c.exp);
    detail::expand_significand<DecimalType>(rhs_c.sig, rhs_c.exp);

    if (BOOST_DECIMAL_UNLIKELY(lhs_c.sig == 0U))
    {
        using sig_type = typename decltype(lhs_c)::significand_type;
        return DecimalType{sig_type{0U}, ideal_exp, sign};
    }

    // Scale by 10^(p-1) when lhs >= rhs, thus the quotient has exactly p digits
    constexpr auto ten_to_p {pow10(static_cast<std::uint64_t>(precision_v<DecimalType>))};
    const bool shorter {lhs_c.sig >= rhs_c.sig};
    const auto big_sig {detail::umul128(static_cast<std::uint64_t>(lhs_c.sig), select_bits(shorter, ten_to_p / 10U, ten_to_p))};
    const auto divisor {static_cast<std::uint64_t>(rhs_c.sig)};
    const auto q_wide {big_sig / static_cast<unsigned_int128_type>(divisor)};
    const std::uint64_t remainder {big_sig.low - q_wide.low * divisor};
    const auto res_exp {(lhs_c.exp - static_cast<int>(precision_v<DecimalType>)) - rhs_c.exp + static_cast<int>(shorter)};

    return impl::div_finalize_u128<DecimalType>(static_cast<std::uint64_t>(q_wide.low),
                                                remainder,
                                                divisor,
                                                res_exp,
                                                ideal_exp,
                                                sign,
                                                detail::current_rounding_mode());
}

// d128/fast128 division driver. The pre-scaled dividend lhs.sig * 10^34
// needs u256 because the product may reach ~10^68 (well above uint128 max).
// The quotient itself fits in uint128 (< 10^34) so we narrow before the
// finalizer.
template <typename DecimalType, typename T>
BOOST_DECIMAL_CUDA_CONSTEXPR auto d128_generic_div_impl(const T& lhs, const T& rhs, const bool sign) noexcept -> DecimalType
{
    auto lhs_c {lhs.to_components()};
    auto rhs_c {rhs.to_components()};
    const auto ideal_exp {static_cast<int>(lhs_c.exp) - static_cast<int>(rhs_c.exp)};

    detail::expand_significand<DecimalType>(lhs_c.sig, lhs_c.exp);
    detail::expand_significand<DecimalType>(rhs_c.sig, rhs_c.exp);

    if (BOOST_DECIMAL_UNLIKELY(lhs_c.sig == int128::uint128_t{0U}))
    {
        return DecimalType{int128::uint128_t{0U}, ideal_exp, sign};
    }

    // Scale by 10^(p-1) when lhs >= rhs, thus the quotient has exactly p digits
    constexpr auto ten_to_p {pow10(int128::uint128_t{static_cast<std::uint64_t>(precision_v<DecimalType>)})};
    constexpr auto ten_to_p_minus_1 {pow10(int128::uint128_t{static_cast<std::uint64_t>(precision_v<DecimalType> - 1)})};
    const bool shorter {lhs_c.sig >= rhs_c.sig};
    const auto big_sig {detail::umul256(lhs_c.sig, select_bits(shorter, ten_to_p_minus_1, ten_to_p))};
    const auto divisor {rhs_c.sig};
    const auto dr {impl::div_mod(big_sig, divisor)};

    const int128::uint128_t q {dr.quotient.bytes[1], dr.quotient.bytes[0]};
    const int128::uint128_t r {dr.remainder.bytes[1], dr.remainder.bytes[0]};
    const auto res_exp {lhs_c.exp - rhs_c.exp - static_cast<int>(precision_v<DecimalType>) + static_cast<int>(shorter)};

    return impl::div_finalize_u256<DecimalType>(q, r, divisor, res_exp, ideal_exp, sign, detail::current_rounding_mode());
}

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

} // namespace detail
} // namespace decimal
} // namespace boost

#endif //BOOST_DECIMAL_DETAIL_DIV_IMPL_HPP
