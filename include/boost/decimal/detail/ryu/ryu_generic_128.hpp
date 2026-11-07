// Copyright 2018 - 2023 Ulf Adams
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_RYU_RYU_GENERIC_128_HPP
#define BOOST_DECIMAL_DETAIL_RYU_RYU_GENERIC_128_HPP

#include <boost/decimal/detail/ryu/generic_128.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/binary_bits.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/bit_cast.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cinttypes>
#include <cstdint>
#endif

#if !defined(BOOST_DECIMAL_DISABLE_IOSTREAM)
#ifdef BOOST_DECIMAL_DEBUG_RYU
#  include <iostream>
#endif
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace ryu {

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int32_t fd128_exceptional_exponent = 0x7FFFFFFF;
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE unsigned_128_type one = 1;

#ifdef _MSC_VER
#  pragma warning(push)
#  pragma warning(disable : 4324) // Structure was padded due to alignment specifier
#endif

// Internal use only
#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpadded"
#endif

struct floating_decimal_128
{
    unsigned_128_type mantissa;
    int32_t exponent;
    bool sign;
    // Whether the exact value is above (1) or below (-1) the digits, or 0 when they are exact, so that the digits
    // can round correctly to a narrower precision
    int32_t direction;
};

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#ifdef BOOST_DECIMAL_DEBUG_RYU
static char* s(unsigned_128_type v) {
  int len = num_digits(v);
  char* b = (char*) malloc((len + 1) * sizeof(char));
  for (int i = 0; i < len; i++) {
    const uint32_t c = (uint32_t) (v % 10);
    v /= 10;
    b[len - 1 - i] = (char) ('0' + c);
  }
  b[len] = 0;
  return b;
}
#endif

// x / 10 is the top 128 bits of x * ceil(2^131 / 10) shifted right by 3, exact for all x as in 64 bits.
// The division of uint128_t by 10 calls a long division.
constexpr auto div10(const unsigned_128_type x) noexcept -> unsigned_128_type
{
    constexpr std::uint64_t high {UINT64_C(0xCCCCCCCCCCCCCCCC)};
    constexpr std::uint64_t low {UINT64_C(0xCCCCCCCCCCCCCCCD)};
    const auto ll {unsigned_128_type {x.low} * low};
    const auto lh {unsigned_128_type {x.low} * high};
    const auto hl {unsigned_128_type {x.high} * low};
    const auto mid {unsigned_128_type {ll.high} + unsigned_128_type {lh.low} + unsigned_128_type {hl.low}};
    const auto top {unsigned_128_type {x.high} * high + unsigned_128_type {lh.high} + unsigned_128_type {hl.high} + unsigned_128_type {mid.high}};
    return top >> 3U;
}

// x mod 10^8 from x mod 2^8 and (x >> 8) mod 5^8, where 2^64 mod 5^8 is 176616.
// It avoids the 128-bit division.
constexpr auto mod_pow10_8(const unsigned_128_type x) noexcept -> std::uint64_t
{
    constexpr std::uint64_t pow5_8 {390625U};
    const auto y {x >> 8U};
    const auto r {((y.high % pow5_8) * UINT64_C(176616) + y.low % pow5_8) % pow5_8};
    return (r << 8U) | (x.low & 255U);
}

// A struct of arrays has one definition of each table, also in C++14
template <bool b>
struct inv_pow5_table_imp
{
    // 5^-k mod 2^128 for k from 0 to 8, which divides a multiple of 5^k exactly
    static constexpr unsigned_128_type values[9] = {
        {UINT64_C(0x0000000000000000), UINT64_C(0x0000000000000001)},
        {UINT64_C(0xCCCCCCCCCCCCCCCC), UINT64_C(0xCCCCCCCCCCCCCCCD)},
        {UINT64_C(0x28F5C28F5C28F5C2), UINT64_C(0x8F5C28F5C28F5C29)},
        {UINT64_C(0x6E978D4FDF3B645A), UINT64_C(0x1CAC083126E978D5)},
        {UINT64_C(0x495182A9930BE0DE), UINT64_C(0xD288CE703AFB7E91)},
        {UINT64_C(0xDB76B3BB83CF2CF9), UINT64_C(0x5D4E8FB00BCBE61D)},
        {UINT64_C(0xC57E23F24D8FD5CB), UINT64_C(0x790FB65668C26139)},
        {UINT64_C(0xC1193A63A91CC45B), UINT64_C(0xE5032477AE8D46A5)},
        {UINT64_C(0xF36B7213EE9F5A78), UINT64_C(0xC767074B22E90E21)},
    };
};

#if !(defined(__cpp_inline_variables) && __cpp_inline_variables >= 201606L) && (!defined(_MSC_VER) || _MSC_VER != 1900)

template <bool b> constexpr unsigned_128_type inv_pow5_table_imp<b>::values[];

#endif

using inv_pow5_table = inv_pow5_table_imp<true>;

inline constexpr auto generic_binary_to_decimal(
        const unsigned_128_type bits,
        const std::uint32_t mantissaBits, const std::uint32_t exponentBits, const bool explicitLeadingBit) noexcept -> struct floating_decimal_128
{
    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("IN=");
    for (int32_t bit = 127; bit >= 0; --bit)
    {
        printf("%u", (uint32_t) ((bits >> bit) & 1));
    }
    printf("\n");
    #endif

    const std::uint32_t local_bias = (1u << (exponentBits - 1)) - 1;
    const bool ieeeSign = ((bits >> (mantissaBits + exponentBits)) & 1U) != 0U;
    const unsigned_128_type ieeeMantissa = bits & ((one << mantissaBits) - 1U);
    const auto ieeeExponent = static_cast<std::uint32_t>((bits >> mantissaBits) & ((one << exponentBits) - 1u));

    if (ieeeExponent == 0 && ieeeMantissa == 0U)
    {
        struct floating_decimal_128 fd {0, 0, ieeeSign, 0};
        return fd;
    }
    if (ieeeExponent == ((1u << exponentBits) - 1u))
    {
        struct floating_decimal_128 fd {};
        fd.mantissa = explicitLeadingBit ? ieeeMantissa & ((one << (mantissaBits - 1)) - 1U) : ieeeMantissa;
        fd.exponent = fd128_exceptional_exponent;
        fd.sign = ieeeSign;
        return fd;
    }

    std::int32_t e2 {};
    unsigned_128_type m2 {};
    // We subtract 2 in all cases so that the bounds computation has 2 additional bits.
    if (explicitLeadingBit)
    {
        // mantissaBits includes the explicit leading bit, so we need to correct for that here.
        if (BOOST_DECIMAL_UNLIKELY(ieeeExponent == 0))
        {
            e2 = static_cast<std::int32_t>(1 - local_bias - mantissaBits + 1 - 2); // LCOV_EXCL_LINE
        }
        else
        {
            e2 = static_cast<std::int32_t>(ieeeExponent - local_bias - mantissaBits + 1 - 2);
        }
        m2 = ieeeMantissa;
    }
    else
    {
        if (BOOST_DECIMAL_UNLIKELY(ieeeExponent == 0))
        {
            e2 = static_cast<std::int32_t>(1 - local_bias - mantissaBits - 2);  // LCOV_EXCL_LINE
            m2 = ieeeMantissa;                                                  // LCOV_EXCL_LINE
        } else
        {
            e2 = static_cast<std::int32_t>(ieeeExponent - local_bias - mantissaBits - 2U);
            m2 = (one << mantissaBits) | ieeeMantissa;
        }
    }
    const bool even = (m2 & 1U) == 0U;
    const bool acceptBounds = even;

    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("-> %s %s * 2^%d\n", ieeeSign ? "-" : "+", s(m2), e2 + 2);
    #endif

    // Step 2: Determine the interval of legal decimal representations.
    const unsigned_128_type mv = 4U * m2;
    // Implicit bool -> int conversion. True is 1, false is 0.
    const std::uint32_t mmShift =
            (ieeeMantissa != (explicitLeadingBit ? one << (mantissaBits - 1) : static_cast<unsigned_128_type>(0)))
            || (ieeeExponent == 0);

    // Step 3: Convert to a decimal power base using 128-bit arithmetic.
    unsigned_128_type vr {};
    unsigned_128_type vp {};
    unsigned_128_type vm {};
    int32_t e10 {};
    bool vmIsTrailingZeros = false;
    bool vrIsTrailingZeros = false;
    if (e2 >= 0)
    {
        // I tried special-casing q == 0, but there was no effect on performance.
        // This expression is slightly faster than max(0, log10Pow2(e2) - 1).
        const std::uint32_t q = log10Pow2(e2) - (e2 > 3);
        e10 = static_cast<int32_t>(q);
        const std::int32_t k = BOOST_DECIMAL_POW5_INV_BITCOUNT + static_cast<std::int32_t>(pow5bits(q)) - 1;
        const std::int32_t i = -e2 + static_cast<std::int32_t>(q) + k;
        std::uint64_t pow5[4] {};
        generic_computeInvPow5(q, pow5);
        vr = mulShiftAll(m2, pow5, static_cast<std::uint32_t>(i), mmShift, vp, vm);

        #ifdef BOOST_DECIMAL_DEBUG_RYU
        printf("%s * 2^%d / 10^%d\n", s(mv), e2, q);
        printf("V+=%s\nV =%s\nV-=%s\n", s(vp), s(vr), s(vm));
        #endif

        // floor(log_5(2^128)) = 55, this is very conservative
        if (q <= 55)
        {
            // Only one of mp, mv, and mm can be a multiple of 5, if any. With q == 0, vr is the exact value.
            vrIsTrailingZeros = q == 0U;
            if (multipleOfPowerOf5(mv, 1))
            {
                vrIsTrailingZeros = multipleOfPowerOf5(mv, q);
            }
            else if (acceptBounds)
            {
                // Same as min(e2 + (~mm & 1), pow5Factor(mm)) >= q
                // <=> e2 + (~mm & 1) >= q && pow5Factor(mm) >= q
                // <=> true && pow5Factor(mm) >= q, since e2 >= q.
                vmIsTrailingZeros = multipleOfPowerOf5(mv - 1U - mmShift, q);
            }
            else
            {
                // Same as min(e2 + 1, pow5Factor(mp)) >= q.
                vp -= multipleOfPowerOf5(mv + 2U, q);
            }
        }
    }
    else
    {
        // This expression is slightly faster than max(0, log10Pow5(-e2) - 1).
        const uint32_t q = log10Pow5(-e2) - static_cast<uint32_t>(-e2 > 1);
        e10 = static_cast<int32_t>(q) + e2;
        const int32_t i = -e2 - static_cast<int32_t>(q);
        const int32_t k = static_cast<int32_t>(pow5bits(static_cast<uint32_t>(i))) - BOOST_DECIMAL_POW5_BITCOUNT;
        const int32_t j = static_cast<int32_t>(q) - k;
        uint64_t pow5[4U] {};
        generic_computePow5(static_cast<uint32_t>(i), pow5);
        vr = mulShiftAll(m2, pow5, static_cast<std::uint32_t>(j), mmShift, vp, vm);

        #ifdef BOOST_DECIMAL_DEBUG_RYU
        printf("%s * 5^%d / 10^%d\n", s(mv), -e2, q);
        printf("%d %d %d %d\n", q, i, k, j);
        printf("V+=%s\nV =%s\nV-=%s\n", s(vp), s(vr), s(vm));
        #endif

        if (q <= 1)
        {
            // {vr,vp,vm} is trailing zeros if {mv,mp,mm} has at least q trailing 0 bits.
            // mv = 4 m2, so it always has at least two trailing 0 bits.
            vrIsTrailingZeros = true;
            if (acceptBounds)
            {
                // mm = mv - 1 - mmShift, so it has 1 trailing 0 bit iff mmShift == 1.
                vmIsTrailingZeros = mmShift == 1;
            }
            else
            {
                // mp = mv + 2, so it always has at least one trailing 0 bit.
                --vp;
            }
        }
        else if (q < 127)
        {
            // We need to compute min(ntz(mv), pow5Factor(mv) - e2) >= q
            // <=> ntz(mv) >= q  &&  pow5Factor(mv) - e2 >= q
            // <=> ntz(mv) >= q    (e2 is negative and -e2 >= q)
            // <=> (mv & ((1 << q) - 1)) == 0
            // We also need to make sure that the left shift does not overflow.
            vrIsTrailingZeros = multipleOfPowerOf2(mv, q);

            #ifdef BOOST_DECIMAL_DEBUG_RYU
            printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
            #endif
        }
    }

    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("e10=%d\n", e10);
    printf("V+=%s\nV =%s\nV-=%s\n", s(vp), s(vr), s(vm));
    printf("vm is trailing zeros=%s\n", vmIsTrailingZeros ? "true" : "false");
    printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
    #endif

    // Step 4: Find the shortest decimal representation in the interval of legal representations.
    uint32_t removed = 0;
    uint8_t lastRemovedDigit = 0;
    unsigned_128_type output {};

    // A multiple of 10^k lies in (vm, vp] if vp mod 10^k < vp - vm. With the remainders
    // mod 10^8, the digits go 8 at a time and then one at a time in 64 bits.
    auto rp {mod_pow10_8(vp)};
    while (vp - vm > rp)
    {
        const auto rr {mod_pow10_8(vr)};
        const auto rm {mod_pow10_8(vm)};
        vmIsTrailingZeros &= rm == 0U;
        vrIsTrailingZeros &= lastRemovedDigit == 0 && rr % 10000000U == 0U;
        lastRemovedDigit = static_cast<uint8_t>(rr / 10000000U);
        vr = ((vr - rr) >> 8U) * inv_pow5_table::values[8];
        vp = ((vp - rp) >> 8U) * inv_pow5_table::values[8];
        vm = ((vm - rm) >> 8U) * inv_pow5_table::values[8];
        rp = mod_pow10_8(vp);
        removed += 8U;
    }

    // vm <= vr <= vp and vp - vm <= vp mod 10^8, so vr and vm mod 10^8 follow
    const auto delta {(vp - vm).low};
    auto qp {rp};
    auto qr {rp - (vp - vr).low};
    auto qm {rp - delta};
    std::uint64_t ap {};
    std::uint64_t ar {};
    std::uint64_t am {};
    std::uint64_t pow10 {1U};
    std::uint32_t k {};
    while (ap + qp % 10U * pow10 < delta)
    {
        vrIsTrailingZeros &= lastRemovedDigit == 0;
        lastRemovedDigit = static_cast<uint8_t>(qr % 10U);
        ap += qp % 10U * pow10;
        ar += static_cast<std::uint64_t>(lastRemovedDigit) * pow10;
        am += qm % 10U * pow10;
        qp /= 10U;
        qr /= 10U;
        qm /= 10U;
        pow10 *= 10U;
        ++k;
    }
    // floor(vm / 10^k) is the new vr if vr - vm is at most vr mod 10^k
    vmIsTrailingZeros &= am == 0U;
    auto vrIsVm {(vr - vm).low <= ar};
    vr = ((vr - ar) >> k) * inv_pow5_table::values[k];
    removed += k;

    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("V+=%s\nV =%s\nV-=%s\n", s(vp), s(vr), s(vm));
    printf("d-10=%s\n", vmIsTrailingZeros ? "true" : "false");
    #endif

    if (vmIsTrailingZeros)
    {
        vm = ((vm - am) >> k) * inv_pow5_table::values[k];
        auto vm_div10 {div10(vm)};
        while (vm == vm_div10 * 10U)
        {
            const auto vr_div10 {div10(vr)};
            vrIsTrailingZeros &= lastRemovedDigit == 0;
            lastRemovedDigit = static_cast<uint8_t>(vr - vr_div10 * 10U);
            vr = vr_div10;
            vm = vm_div10;
            vm_div10 = div10(vm);
            ++removed;
        }
        vrIsVm = vr == vm;
    }

    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("%s %d\n", s(vr), lastRemovedDigit);
    printf("vr is trailing zeros=%s\n", vrIsTrailingZeros ? "true" : "false");
    #endif

    if (vrIsTrailingZeros && (lastRemovedDigit == 5) && (vr % 2U == 0U))
    {
        // Round even if the exact numbers is .....50..0.
        lastRemovedDigit = 4;
    }
    // We need to take vr+1 if vr is outside bounds, or we need to round up.
    output = vr + static_cast<unsigned_128_type>((vrIsVm && (!acceptBounds || !vmIsTrailingZeros)) || (lastRemovedDigit >= 5));
    const auto exp = e10 + static_cast<std::int32_t>(removed);

    #ifdef BOOST_DECIMAL_DEBUG_RYU
    printf("V+=%s\nV =%s\nV-=%s\n", s(vp), s(vr), s(vm));
    printf("O=%s\n", s(output));
    printf("EXP=%d\n", exp);
    #endif

    // The digits are above the exact value when they rounded up, and below it when a removed digit is not zero
    const std::int32_t direction {output != vr ? -1 :
                                  static_cast<std::int32_t>(!vrIsTrailingZeros || lastRemovedDigit != 0)};
    return {output, exp, ieeeSign, direction};
}

template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128(T val) noexcept -> floating_decimal_128;

#if BOOST_DECIMAL_LDBL_BITS == 80

template <>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128<long double>(long double val) noexcept -> floating_decimal_128
{
    return generic_binary_to_decimal(wide_bits(val), 64, 15, true);
}

#elif BOOST_DECIMAL_LDBL_BITS == 128

template <>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128<long double>(long double val) noexcept -> floating_decimal_128
{
    auto bits = bit_cast<unsigned_128_type>(val);

    return generic_binary_to_decimal(bits, 112, 15, false);
}

#elif BOOST_DECIMAL_LDBL_BITS == 0

template <>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128<long double>(long double val) noexcept -> floating_decimal_128
{
    BOOST_DECIMAL_ASSERT_MSG(1==0, "Unsupported configuration");

    auto bits = bit_cast<unsigned_128_type>(val);
    return generic_binary_to_decimal(bits, 112, 15, false);
}

#endif

#if defined(BOOST_DECIMAL_HAS_FLOAT128) && !defined(BOOST_DECIMAL_LDBL_IS_FLOAT128)

template <>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128<__float128>(__float128 val) noexcept -> floating_decimal_128
{
    auto bits = bit_cast<unsigned_128_type>(val);

    return generic_binary_to_decimal(bits, 112, 15, false);
}

#endif

#ifdef BOOST_DECIMAL_HAS_STDFLOAT128

template <>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd128<std::float128_t>(std::float128_t d) noexcept -> floating_decimal_128
{
    auto bits = bit_cast<unsigned_128_type>(d);

    return generic_binary_to_decimal(bits, 112, 15, false);
}

#endif

}}}} // Namespaces

#endif //BOOST_DECIMAL_RYU_GENERIC_128_HPP
