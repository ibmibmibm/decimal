// Copyright 2020-2023 Daniel Lemire
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT80_128_HPP
#define BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT80_128_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/apply_sign.hpp>
#include <boost/decimal/detail/fast_float/compute_float64.hpp>
#include <boost/decimal/detail/integer_search_trees.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <algorithm>
#include <array>
#include <climits>
#include <limits>
#include <cmath>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace fast_float {

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::array<long double, 56> powers_of_ten_ld = {
        1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L,
        1e7L, 1e8L, 1e9L, 1e10L, 1e11L, 1e12L, 1e13L,
        1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L, 1e20L,
        1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L,
        1e28L, 1e29L, 1e30L, 1e31L, 1e32L, 1e33L, 1e34L,
        1e35L, 1e36L, 1e37L, 1e38L, 1e39L, 1e40L, 1e41L,
        1e42L, 1e43L, 1e44L, 1e45L, 1e46L, 1e47L, 1e48L,
        1e49L, 1e50L, 1e51L, 1e52L, 1e53L, 1e54L, 1e55L
};

// The top 128 bits of 5^(300 * k) for k = -17 to 16, rounded down, with p such that
// 10^(300 * k) is near the table value times 2^p
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t coarse_power_high[] = {
        UINT64_C(0x8fae1896aa9cbc6b), UINT64_C(0xd68bd3c92066a797), UINT64_C(0xa02eba538ac25a32),
        UINT64_C(0xef3023b80a732d93), UINT64_C(0xb294951d99cd3295), UINT64_C(0x85547fa14b14e46f),
        UINT64_C(0xc71763fcc709b18f), UINT64_C(0x94a4cf2019d7ba15), UINT64_C(0xddf5489f3473ea88),
        UINT64_C(0xa5b763b319d7f1dc), UINT64_C(0xf773878e7ec7dd45), UINT64_C(0xb8bff4a88f1fb463),
        UINT64_C(0x89efad01750d64bd), UINT64_C(0xcdf829eaaf012977), UINT64_C(0x99c76be80626e8cc),
        UINT64_C(0xe5a04aa62b553e99), UINT64_C(0xab70fe17c79ac6ca), UINT64_C(0x8000000000000000),
        UINT64_C(0xbf21e44003acdd2c), UINT64_C(0x8eb39714297efb27), UINT64_C(0xd515c4344c1e8ef2),
        UINT64_C(0x9f17732dfac7617b), UINT64_C(0xed8f1d9ba3850fa9), UINT64_C(0xb15d3a58cd9f7620),
        UINT64_C(0x846c09b028ae0395), UINT64_C(0xc5bc4672073224f7), UINT64_C(0x93a1a62148b73c1e),
        UINT64_C(0xdc724cd44411af53), UINT64_C(0xa496769b32506bee), UINT64_C(0xf5c419447c50865f),
        UINT64_C(0xb77dd84175a2aadb), UINT64_C(0x88ff2f2bade74531), UINT64_C(0xcc910e78664e5136),
        UINT64_C(0x98bb4ee309f04d45)
};

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::uint64_t coarse_power_low[] = {
        UINT64_C(0x9f9df6607da0466d), UINT64_C(0x326cb526b3747638), UINT64_C(0xe85f0d271ebf22c6),
        UINT64_C(0xf5a7800f23ef67b8), UINT64_C(0x1a2cb1cc24e829ad), UINT64_C(0x1381a641a3d35b14),
        UINT64_C(0x509c2b58aa2fcd73), UINT64_C(0x711ae6040de62000), UINT64_C(0xe615da33ac7d627e),
        UINT64_C(0x0a0f429d93058121), UINT64_C(0x2b566ef4caf507b0), UINT64_C(0xb02d9653eb387400),
        UINT64_C(0x19ddee9ee91b7777), UINT64_C(0x36f89dc4f34187c3), UINT64_C(0x549dd110d4c13756),
        UINT64_C(0xf74ce198d654ca94), UINT64_C(0x6dbd630a48aaf406), UINT64_C(0x0000000000000000),
        UINT64_C(0xe0470a63e6bd56c3), UINT64_C(0xbcafdba96ebcb609), UINT64_C(0x915861c376f9d550),
        UINT64_C(0x3b7dffb204b6d932), UINT64_C(0xab3792186b8256cb), UINT64_C(0xe79894cf86a0d163),
        UINT64_C(0x04f609974dd3ffe9), UINT64_C(0xb2c46d6d298a0658), UINT64_C(0xe4df06eb48627ae8),
        UINT64_C(0x37b24bf1b205e08c), UINT64_C(0xe15dbb1c1b7a2785), UINT64_C(0x30f451703325e354),
        UINT64_C(0x49e6afe3d2b682d8), UINT64_C(0xc9ac50475e25293a), UINT64_C(0x10a1bd37a7340862),
        UINT64_C(0x5a050b215eebc516)
};

BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE std::int32_t coarse_power_exponent[] = {
        -17069, -16073, -15076,
        -14080, -13083, -12086,
        -11090, -10093, -9097,
        -8100, -7104, -6107,
        -5110, -4114, -3117,
        -2121, -1124, -127,
        869, 1866, 2862,
        3859, 4855, 5852,
        6849, 7845, 8842,
        9838, 10835, 11831,
        12828, 13825, 14821,
        15818
};

template <typename Unsigned_Integer>
constexpr auto fast_path(const std::int64_t q, const Unsigned_Integer &w, bool negative) noexcept -> long double
{
    // The general idea is as follows.
    // if 0 <= s <= 2^64 and if 10^0 <= p <= 10^27
    // Both s and p can be represented exactly
    // because of this s*p and s/p will produce
    // correctly rounded values

    auto ld = static_cast<long double>(w);

    if (q < 0)
    {
        ld /= powers_of_ten_ld[static_cast<std::size_t>(-q)];
    }
    else
    {
        ld *= powers_of_ten_ld[static_cast<std::size_t>(q)];
    }

    if (negative)
    {
        ld = -ld;
    }

    return ld;
}

// The correctly rounded value of w * 10^q with digits significand bits, round half to even.
// The value is mantissa * 2^exponent, or an infinity.
struct long_double_parts
{
    int128::uint128_t mantissa;
    std::int64_t exponent;
    bool infinite;
};

// bits * 2^exponent rounded to digits bits, or to fewer bits for a subnormal value.
// bits has more than digits bits, and sticky is set if the exact value is above bits * 2^exponent.
template <int digits, int min_exponent, int max_exponent>
constexpr auto round_long_double_parts(const int128::uint128_t bits, const std::int64_t exponent, const bool sticky) noexcept -> long_double_parts
{
    constexpr std::int64_t min_bit {min_exponent - digits};

    const auto bit_count {static_cast<std::int64_t>(128 - countl_zero(bits))};
    const auto last_bit {(std::max)(exponent + bit_count - digits, min_bit)};
    const auto drop {last_bit - exponent};
    if (drop > bit_count)
    {
        return {0U, 0, false};
    }

    auto mantissa {bits >> static_cast<unsigned>(drop)};
    const auto half {int128::uint128_t {1U} << static_cast<unsigned>(drop - 1)};
    const auto rest {bits & ((half << 1U) - 1U)};
    if (rest > half || (rest == half && (sticky || (mantissa & 1U) != 0U)))
    {
        ++mantissa;
    }

    if (last_bit + 128 - countl_zero(mantissa) > max_exponent)
    {
        return {0U, 0, true};
    }

    return {mantissa, last_bit, false};
}

// x * 2^exponent rounded as round_long_double_parts does, for x with the top bit set.
// For a normal result, the kept bits are the top digits bits of x.
template <int digits, int min_exponent, int max_exponent>
constexpr auto round_top_bits(const int128::uint128_t x, const std::int64_t exponent, const bool sticky) noexcept -> long_double_parts
{
    constexpr auto drop {128 - digits};
    if (exponent + drop < min_exponent - digits)
    {
        constexpr auto shift {drop - 2};
        constexpr auto mask {(int128::uint128_t {1U} << static_cast<unsigned>(shift)) - 1U};
        return round_long_double_parts<digits, min_exponent, max_exponent>(x >> static_cast<unsigned>(shift), exponent + shift, sticky || (x & mask) != 0U);
    }

    constexpr auto half {int128::uint128_t {1U} << static_cast<unsigned>(drop - 1)};
    auto mantissa {x >> static_cast<unsigned>(drop)};
    const auto rest {x & ((half << 1U) - 1U)};
    const auto up {static_cast<unsigned>(rest > half) | (static_cast<unsigned>(rest == half) & (static_cast<unsigned>(sticky) | static_cast<unsigned>(mantissa.low & 1U)))};
    mantissa += up;

    if (exponent + drop + 128 - countl_zero(mantissa) > max_exponent)
    {
        return {0U, 0, true};
    }

    return {mantissa, exponent + drop, false};
}

// The top 128 bits of a product, which is bits * 2^drop and the bits below
struct top_product
{
    int128::uint128_t bits;
    int drop;
    bool sticky;
};

// x * (high, low) for x and high with the top bit set, so the product has 255 or 256 bits.
// The products of 64-bit words keep the values in registers, unlike the u256 of umul256.
constexpr auto multiply_top(const int128::uint128_t x, const std::uint64_t high, const std::uint64_t low) noexcept -> top_product
{
    const auto ll {umul128(x.low, low)};
    const auto lh {umul128(x.low, high)};
    const auto hl {umul128(x.high, low)};
    const auto hh {umul128(x.high, high)};
    const auto mid {int128::uint128_t {ll.high} + int128::uint128_t {lh.low} + int128::uint128_t {hl.low}};
    const auto top {hh + int128::uint128_t {lh.high} + int128::uint128_t {hl.high} + int128::uint128_t {mid.high}};

    // Shift left by one if the top bit is clear, without a branch: the bit is random
    const auto n {static_cast<unsigned>(1U - (top.high >> 63U))};
    const int128::uint128_t bits {(top.high << n) | ((top.low >> 63U) & n), (top.low << n) | ((mid.low >> 63U) & n)};
    return {bits, 128 - static_cast<int>(n), ((mid.low << n) | ll.low) != 0U};
}

// x * 5^step from the table
constexpr auto multiply_top(const int128::uint128_t x, const std::int64_t step) noexcept -> top_product
{
    const auto index {static_cast<std::size_t>(step - smallest_power_of_five)};
    return multiply_top(x, significand_64[index], significand_128[index]);
}

// 10^p is at most half of the smallest subnormal value for p <= long_double_min_power10
constexpr auto long_double_min_power10(const int digits, const int min_exponent) noexcept -> std::int64_t
{
    return -((1 - (min_exponent - digits)) * INT64_C(30103) + 99999) / 100000;
}

// 10^p is above the largest value for p > long_double_max_power10
constexpr auto long_double_max_power10(const int max_exponent) noexcept -> std::int64_t
{
    return max_exponent * INT64_C(30103) / 100000;
}

// For binary formats with min_exponent and max_exponent as in std::numeric_limits
template <int digits, int min_exponent, int max_exponent>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compute_long_double_parts(const std::int64_t q, const int128::uint128_t w) noexcept -> long_double_parts
{
    // w has at most 39 digits
    if (w == 0U || q + 39 <= long_double_min_power10(digits, min_exponent))
    {
        return {0U, 0, false};
    }
    if (q > long_double_max_power10(max_exponent))
    {
        return {0U, 0, true};
    }

    // w * 10^q is x * 10^q * 2^exponent, with the top bit of x set
    const auto lz {countl_zero(w)};
    auto x {w << static_cast<unsigned>(lz)};
    std::int64_t exponent {-lz};

    constexpr auto shift {128 - (digits + 2)};
    constexpr auto mask {(int128::uint128_t {1U} << static_cast<unsigned>(shift)) - 1U};

    // The power table has 5^q exactly for 0 <= q <= 55, so the product is exact
    if (q >= 0 && q <= 55)
    {
        const auto top {multiply_top(x, q)};
        exponent += binary_power(static_cast<std::int32_t>(q)) - 190 + top.drop;
        return round_top_bits<digits, min_exponent, max_exponent>(top.bits, exponent, top.sticky);
    }

    // Else multiply by 10^q from the table, with the top 128 bits of each product.
    // Each step adds a relative error below 2^-126, so the result is off by less than 4 * steps.
    std::int64_t steps {1};
    auto rest {q};
    if (q < smallest_power_of_five || q > 308)
    {
        // 10^q is 10^(300 * k) from the coarse table times 10^rest with 0 <= rest < 300
        const auto k {q >= 0 ? q / 300 : -((299 - q) / 300)};
        rest -= 300 * k;
        const auto index {static_cast<std::size_t>(k + 17)};
        const auto top {multiply_top(x, coarse_power_high[index], coarse_power_low[index])};
        x = top.bits;
        exponent += coarse_power_exponent[index] + top.drop;
        ++steps;
    }

    const auto top {multiply_top(x, rest)};
    x = top.bits;
    exponent += binary_power(static_cast<std::int32_t>(rest)) - 190 + top.drop;

    const auto dropped {x & mask};
    const auto error {static_cast<std::uint64_t>(4 * steps)};
    if (dropped > error && dropped < mask - error)
    {
        return round_top_bits<digits, min_exponent, max_exponent>(x, exponent, true);
    }

    exponent += shift;

    // Else w * 10^q is next to boundary * 2^exponent, and an exact compare decides
    const auto boundary {(x >> static_cast<unsigned>(shift)) + (dropped > (mask >> 1U) ? 1U : 0U)};
    compare_bigint<200> value {w};
    compare_bigint<200> limit {boundary};
    if (q >= 0)
    {
        value.mul_pow5(q);
    }
    else
    {
        limit.mul_pow5(-q);
    }

    const auto two_shift {q - exponent};
    if (two_shift >= 0)
    {
        value.shl(two_shift);
    }
    else
    {
        limit.shl(-two_shift);
    }

    const auto order {value.compare(limit)};
    return round_long_double_parts<digits, min_exponent, max_exponent>(order < 0 ? boundary - 1U : boundary, exponent, order != 0);
}

// mantissa * 2^exponent, exact if the value fits in T
template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_long_double(const long_double_parts& parts, const bool negative) noexcept -> T
{
    if (parts.infinite)
    {
        return negative ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
    }

    #if !defined(BOOST_DECIMAL_NO_CONSTEVAL_DETECTION) && (BOOST_DECIMAL_LDBL_BITS == 128 || (BOOST_DECIMAL_LDBL_BITS == 80 && defined(__SIZEOF_LONG_DOUBLE__) && __SIZEOF_LONG_DOUBLE__ == 16 && BOOST_DECIMAL_ENDIAN_LITTLE_BYTE))
    #ifdef BOOST_DECIMAL_HAS_CONSTEXPR_BITCAST
    if (std::is_same<T, long double>::value && !BOOST_DECIMAL_IS_CONSTANT_EVALUATED(parts))
    #else
    if (std::is_same<T, long double>::value)
    #endif
    {
        #if BOOST_DECIMAL_LDBL_BITS == 80
        // The x87 format keeps the leading bit, so the bits come from the parts: a 64-bit mantissa below 2^63
        // is subnormal with the exponent field 0, and a carry to 2^64 is 2^63 one binade up
        const auto carry {static_cast<unsigned>(parts.mantissa.high)};
        const auto mantissa {(parts.mantissa.low >> carry) | (static_cast<std::uint64_t>(carry) << 63U)};
        const auto field {static_cast<std::uint64_t>(parts.exponent + 16446 + carry) & (UINT64_C(0) - (mantissa >> 63U))};
        const std::uint64_t words[2] {mantissa, field | (static_cast<std::uint64_t>(negative) << 15U)};
        return bit_cast<T>(words);
        #else
        // The same for the binary128 format, which drops the leading bit 2^112
        const auto carry {static_cast<unsigned>(parts.mantissa.high >> 49U)};
        const auto mantissa {parts.mantissa >> carry};
        const auto field {static_cast<std::uint64_t>(parts.exponent + 16495 + carry) & (UINT64_C(0) - (mantissa.high >> 48U))};
        const auto high {(mantissa.high & ((UINT64_C(1) << 48U) - 1U)) | (field << 48U) | (static_cast<std::uint64_t>(negative) << 63U)};
        return bit_cast<T>(int128::uint128_t {high, mantissa.low});
        #endif
    }
    #endif

    auto value {static_cast<T>(parts.mantissa.high) * static_cast<T>(UINT64_C(1) << 63U) * 2 + static_cast<T>(parts.mantissa.low)};

    // Each partial product has the same significand bits, so all are exact.
    // A double holds 2^e exactly for |e| <= 1000.
    constexpr std::int64_t step {1000};
    auto exponent {parts.exponent};
    for (; exponent > step; exponent -= step)
    {
        value *= static_cast<T>(bit_cast<double>(static_cast<std::uint64_t>(1023 + step) << 52U));
    }
    for (; exponent < -step; exponent += step)
    {
        value *= static_cast<T>(bit_cast<double>(static_cast<std::uint64_t>(1023 - step) << 52U));
    }
    value *= static_cast<T>(bit_cast<double>(static_cast<std::uint64_t>(1023 + exponent) << 52U));

    return negative ? -value : value;
}

template <typename Unsigned_Integer>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compute_float80_128(const std::int64_t q, const Unsigned_Integer &w,
                                                       const bool negative, bool &success) noexcept -> long double
{
    using limits = std::numeric_limits<long double>;

    success = true;

    #if BOOST_DECIMAL_LDBL_BITS == 80 || BOOST_DECIMAL_LDBL_BITS == 128
    // Clinger: if w and 10^|q| are exact, one multiply or divide rounds w * 10^q correctly
    constexpr std::int64_t max_exponent {BOOST_DECIMAL_LDBL_BITS == 80 ? 27 : 48};
    if (q >= -max_exponent && q <= max_exponent && static_cast<int128::uint128_t>(w) <= (int128::uint128_t {1U} << static_cast<unsigned>(limits::digits)) - 1U)
    {
        return fast_path(q, w, negative);
    }
    #endif

    const auto parts {compute_long_double_parts<limits::digits, limits::min_exponent, limits::max_exponent>(q, static_cast<int128::uint128_t>(w))};
    return to_long_double<long double>(parts, negative);
}

} //namespace fast_float
} //namespace detail
} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT80_128_HPP
