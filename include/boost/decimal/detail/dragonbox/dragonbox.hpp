// Copyright 2024-2025 Junekey Jeon, Toby Bell
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Dragonbox's subproject/simple/include/simple_dragonbox.h (github.com/jk-jeon/dragonbox) without its string output.

#ifndef BOOST_DECIMAL_DETAIL_DRAGONBOX_DRAGONBOX_HPP
#define BOOST_DECIMAL_DETAIL_DRAGONBOX_DRAGONBOX_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/bit_layouts.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/remove_trailing_zeros.hpp>
#include <boost/decimal/detail/ryu/ryu_generic_128.hpp>
#include <boost/decimal/detail/fast_float/significand_tables.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace dragonbox {

template <typename Unsigned>
struct floating_decimal
{
    Unsigned mantissa;
    std::int32_t exponent;
    bool sign;
};

struct cache_entry
{
    std::uint64_t high;
    std::uint64_t low;
};

template <typename Unsigned>
struct mul_result
{
    Unsigned integer_part;
    bool is_integer;
};

struct parity_result
{
    bool parity;
    bool is_integer;
};

constexpr auto umul128(const std::uint64_t x, const std::uint64_t y) noexcept -> int128::uint128_t
{
    return static_cast<int128::uint128_t>(x) * y;
}

constexpr auto umul192_upper128(const std::uint64_t x, const cache_entry y) noexcept -> cache_entry
{
    const auto high {umul128(x, y.high)};
    const auto middle {umul128(x, y.low).high};
    const auto low {high.low + middle};
    return {high.high + static_cast<std::uint64_t>(low < middle), low};
}

constexpr auto umul192_lower128(const std::uint64_t x, const cache_entry y) noexcept -> cache_entry
{
    const auto high_low {umul128(x, y.low)};
    return {x * y.high + high_low.high, high_low.low};
}

constexpr auto umul96_upper64(const std::uint32_t x, const std::uint64_t y) noexcept -> std::uint64_t
{
    return std::uint64_t {x} * (y >> 32U) + ((std::uint64_t {x} * static_cast<std::uint32_t>(y)) >> 32U);
}

constexpr auto floor_log10_pow2(const int e) noexcept -> int
{
    return (e * 315653) >> 20;
}

constexpr auto floor_log2_pow10(const int e) noexcept -> int
{
    return (e * 1741647) >> 19;
}

constexpr auto floor_log10_pow2_minus_log10_4_over_3(const int e) noexcept -> int
{
    return (e * 631305 - 261663) >> 21;
}

// A struct of arrays has one definition of each table, see https://github.com/boostorg/math/issues/923
template <bool b>
struct cache_table_imp
{
    static constexpr std::uint64_t binary32[] = {
        0x81ceb32c4b43fcf5, 0xa2425ff75e14fc32, 0xcad2f7f5359a3b3f, 0xfd87b5f28300ca0e,
        0x9e74d1b791e07e49, 0xc612062576589ddb, 0xf79687aed3eec552, 0x9abe14cd44753b53,
        0xc16d9a0095928a28, 0xf1c90080baf72cb2, 0x971da05074da7bef, 0xbce5086492111aeb,
        0xec1e4a7db69561a6, 0x9392ee8e921d5d08, 0xb877aa3236a4b44a, 0xe69594bec44de15c,
        0x901d7cf73ab0acda, 0xb424dc35095cd810, 0xe12e13424bb40e14, 0x8cbccc096f5088cc,
        0xafebff0bcb24aaff, 0xdbe6fecebdedd5bf, 0x89705f4136b4a598, 0xabcc77118461cefd,
        0xd6bf94d5e57a42bd, 0x8637bd05af6c69b6, 0xa7c5ac471b478424, 0xd1b71758e219652c,
        0x83126e978d4fdf3c, 0xa3d70a3d70a3d70b, 0xcccccccccccccccd, 0x8000000000000000,
        0xa000000000000000, 0xc800000000000000, 0xfa00000000000000, 0x9c40000000000000,
        0xc350000000000000, 0xf424000000000000, 0x9896800000000000, 0xbebc200000000000,
        0xee6b280000000000, 0x9502f90000000000, 0xba43b74000000000, 0xe8d4a51000000000,
        0x9184e72a00000000, 0xb5e620f480000000, 0xe35fa931a0000000, 0x8e1bc9bf04000000,
        0xb1a2bc2ec5000000, 0xde0b6b3a76400000, 0x8ac7230489e80000, 0xad78ebc5ac620000,
        0xd8d726b7177a8000, 0x878678326eac9000, 0xa968163f0a57b400, 0xd3c21bcecceda100,
        0x84595161401484a0, 0xa56fa5b99019a5c8, 0xcecb8f27f4200f3a, 0x813f3978f8940985,
        0xa18f07d736b90be6, 0xc9f2c9cd04674edf, 0xfc6f7c4045812297, 0x9dc5ada82b70b59e,
        0xc5371912364ce306, 0xf684df56c3e01bc7, 0x9a130b963a6c115d, 0xc097ce7bc90715b4,
        0xf0bdc21abb48db21, 0x96769950b50d88f5, 0xbc143fa4e250eb32, 0xeb194f8e1ae525fe,
        0x92efd1b8d0cf37bf, 0xb7abc627050305ae, 0xe596b7b0c643c71a, 0x8f7e32ce7bea5c70,
        0xb35dbf821ae4f38c, 0xe0352f62a19e306f};
};

#if !(defined(__cpp_inline_variables) && __cpp_inline_variables >= 201606L) && (!defined(_MSC_VER) || _MSC_VER != 1900)

template <bool b> constexpr std::uint64_t cache_table_imp<b>::binary32[];

#endif

using cache_table = cache_table_imp<true>;

static_assert(sizeof(cache_table::binary32) / sizeof(std::uint64_t) == 78, "one entry per k from -31 to 46");

// 10^k rounded up, from the fast_float table, which rounds down outside 10^-27 to 10^55.
// It starts at 10^-342.
constexpr auto cache_128(const int k) noexcept -> cache_entry
{
    const auto index {static_cast<std::size_t>(k + 342)};
    const auto up {static_cast<std::uint64_t>(k < -27 || k > 55)};
    const auto low {fast_float::significand_128[index] + up};
    return {fast_float::significand_64[index] + static_cast<std::uint64_t>(low < up), low};
}

struct binary32_format
{
    using carrier_uint = std::uint32_t;

    static constexpr int significand_bits {23};
    static constexpr int exponent_bias {-127};
    static constexpr int min_exponent {-126};
    static constexpr int min_k {-31};
    static constexpr int kappa {1};
    static constexpr int shorter_interval_tie_exponent {-35};
    static constexpr std::uint32_t big_divisor {100};
    static constexpr std::uint32_t small_divisor {10};
    static constexpr std::uint64_t divide_magic_number {6554};

    static constexpr auto get_cache(const int minus_k) noexcept -> std::uint64_t
    {
        return cache_table::binary32[static_cast<std::size_t>(-minus_k - min_k)];
    }

    static constexpr auto compute_mul(const std::uint32_t u, const std::uint64_t cache) noexcept
        -> mul_result<std::uint32_t>
    {
        const auto r {umul96_upper64(u, cache)};
        return {static_cast<std::uint32_t>(r >> 32U), static_cast<std::uint32_t>(r) == 0U};
    }

    static constexpr auto compute_delta(const std::uint64_t cache, const int beta) noexcept -> std::uint32_t
    {
        return static_cast<std::uint32_t>(cache >> (63 - beta));
    }

    static constexpr auto compute_mul_parity(const std::uint32_t two_f, const std::uint64_t cache,
                                             const int beta) noexcept -> parity_result
    {
        const auto r {two_f * cache};
        return {((r >> (64 - beta)) & 1U) != 0U, static_cast<std::uint32_t>(r >> (32 - beta)) == 0U};
    }

    static constexpr auto left_endpoint(const std::uint64_t cache, const int beta) noexcept -> std::uint32_t
    {
        return static_cast<std::uint32_t>((cache - (cache >> (significand_bits + 2))) >>
                                          (64 - significand_bits - 1 - beta));
    }

    static constexpr auto right_endpoint(const std::uint64_t cache, const int beta) noexcept -> std::uint32_t
    {
        return static_cast<std::uint32_t>((cache + (cache >> (significand_bits + 1))) >>
                                          (64 - significand_bits - 1 - beta));
    }

    static constexpr auto round_up(const std::uint64_t cache, const int beta) noexcept -> std::uint32_t
    {
        return (static_cast<std::uint32_t>(cache >> (64 - significand_bits - 2 - beta)) + 1U) / 2U;
    }

    static constexpr auto divide_by_big_divisor(const std::uint32_t n) noexcept -> std::uint32_t
    {
        return static_cast<std::uint32_t>((std::uint64_t {n} * 1374389535U) >> 37U);
    }

    static constexpr auto divide_by_10(const std::uint32_t n) noexcept -> std::uint32_t
    {
        return static_cast<std::uint32_t>((std::uint64_t {n} * 429496730U) >> 32U);
    }
};

struct binary64_format
{
    using carrier_uint = std::uint64_t;

    static constexpr int significand_bits {52};
    static constexpr int exponent_bias {-1023};
    static constexpr int min_exponent {-1022};
    static constexpr int kappa {2};
    static constexpr int shorter_interval_tie_exponent {-77};
    static constexpr std::uint64_t big_divisor {1000};
    static constexpr std::uint64_t small_divisor {100};
    static constexpr std::uint64_t divide_magic_number {656};

    static constexpr auto get_cache(const int minus_k) noexcept -> cache_entry
    {
        return cache_128(-minus_k);
    }

    static constexpr auto compute_mul(const std::uint64_t u, const cache_entry cache) noexcept
        -> mul_result<std::uint64_t>
    {
        const auto r {umul192_upper128(u, cache)};
        return {r.high, r.low == 0U};
    }

    static constexpr auto compute_delta(const cache_entry cache, const int beta) noexcept -> std::uint64_t
    {
        return cache.high >> (63 - beta);
    }

    static constexpr auto compute_mul_parity(const std::uint64_t two_f, const cache_entry cache,
                                             const int beta) noexcept -> parity_result
    {
        const auto r {umul192_lower128(two_f, cache)};
        return {((r.high >> (64 - beta)) & 1U) != 0U, ((r.high << beta) | (r.low >> (64 - beta))) == 0U};
    }

    static constexpr auto left_endpoint(const cache_entry cache, const int beta) noexcept -> std::uint64_t
    {
        return (cache.high - (cache.high >> (significand_bits + 2))) >> (64 - significand_bits - 1 - beta);
    }

    static constexpr auto right_endpoint(const cache_entry cache, const int beta) noexcept -> std::uint64_t
    {
        return (cache.high + (cache.high >> (significand_bits + 1))) >> (64 - significand_bits - 1 - beta);
    }

    static constexpr auto round_up(const cache_entry cache, const int beta) noexcept -> std::uint64_t
    {
        return ((cache.high >> (64 - significand_bits - 2 - beta)) + 1U) / 2U;
    }

    static constexpr auto divide_by_big_divisor(const std::uint64_t n) noexcept -> std::uint64_t
    {
        return umul128(n, UINT64_C(4722366482869645214)).high >> 8U;
    }

    static constexpr auto divide_by_10(const std::uint64_t n) noexcept -> std::uint64_t
    {
        return umul128(n, UINT64_C(1844674407370955162)).high;
    }
};

template <typename UInt>
constexpr auto without_trailing_zeros(const UInt significand, const int exponent, const bool sign) noexcept
    -> floating_decimal<UInt>
{
    // Few shortest digit strings end in 0, so one test of the last digit usually skips the removal
    if (significand % 10U != 0U)
    {
        return {significand, exponent, sign};
    }
    const auto removed {remove_trailing_zeros(significand)};
    return {removed.trimmed_number, exponent + static_cast<int>(removed.number_of_removed_zeros), sign};
}

// A zero significand of a normal binary float: the gap below the value is half the gap above it.
template <typename Format>
constexpr auto shorter_interval_case(const int binary_exponent, const bool sign) noexcept
    -> floating_decimal<typename Format::carrier_uint>
{
    const auto minus_k {floor_log10_pow2_minus_log10_4_over_3(binary_exponent)};
    const auto beta {binary_exponent + floor_log2_pow10(-minus_k)};
    const auto cache {Format::get_cache(minus_k)};

    auto xi {Format::left_endpoint(cache, beta)};
    const auto zi {Format::right_endpoint(cache, beta)};
    if (binary_exponent < 2 || binary_exponent > 3)
    {
        ++xi;
    }

    auto significand {Format::divide_by_10(zi)};
    if (significand * 10U >= xi)
    {
        return without_trailing_zeros(significand, minus_k + 1, sign);
    }

    significand = Format::round_up(cache, beta);
    if (significand % 2U != 0U && binary_exponent == Format::shorter_interval_tie_exponent)
    {
        --significand;
    }
    else if (significand < xi)
    {
        ++significand;
    }
    return {significand, minus_k, sign};
}

template <typename Format>
constexpr auto to_decimal(const typename Format::carrier_uint binary_significand, int binary_exponent,
                          const bool sign) noexcept -> floating_decimal<typename Format::carrier_uint>
{
    using carrier_uint = typename Format::carrier_uint;

    const auto is_even {binary_significand % 2U == 0U};
    auto two_fc {binary_significand * 2U};
    if (binary_exponent != 0)
    {
        binary_exponent += Format::exponent_bias - Format::significand_bits;
        if (two_fc == 0U)
        {
            return shorter_interval_case<Format>(binary_exponent, sign);
        }
        two_fc |= carrier_uint {1} << (Format::significand_bits + 1);
    }
    else
    {
        binary_exponent = Format::min_exponent - Format::significand_bits;
    }

    const auto minus_k {floor_log10_pow2(binary_exponent) - Format::kappa};
    const auto cache {Format::get_cache(minus_k)};
    const auto beta {binary_exponent + floor_log2_pow10(-minus_k)};
    const auto deltai {Format::compute_delta(cache, beta)};
    const auto z {Format::compute_mul((two_fc | 1U) << beta, cache)};

    auto significand {Format::divide_by_big_divisor(z.integer_part)};
    auto r {z.integer_part - Format::big_divisor * significand};
    auto found {false};
    if (r < deltai)
    {
        found = r != 0U || !z.is_integer || is_even;
        if (!found)
        {
            --significand;
            r = Format::big_divisor;
        }
    }
    else if (r == deltai)
    {
        const auto x {Format::compute_mul_parity(two_fc - 1U, cache, beta)};
        found = x.parity || (x.is_integer && is_even);
    }
    if (found)
    {
        return without_trailing_zeros(significand, minus_k + Format::kappa + 1, sign);
    }

    significand *= 10U;
    auto dist {r - (deltai / 2U) + (Format::small_divisor / 2U)};
    const auto approx_y_parity {((dist ^ (Format::small_divisor / 2U)) & 1U) != 0U};
    const auto product {static_cast<std::uint32_t>(dist * Format::divide_magic_number)};
    const auto divisible_by_small_divisor {(product & 0xFFFFU) < Format::divide_magic_number};
    dist = product >> 16U;
    significand += dist;
    if (divisible_by_small_divisor)
    {
        const auto y {Format::compute_mul_parity(two_fc, cache, beta)};
        if (y.parity != approx_y_parity || (significand % 2U != 0U && y.is_integer))
        {
            --significand;
        }
    }
    return {significand, minus_k + Format::kappa, sign};
}

// The shortest decimal reading back as the finite float, trailing zeros removed; zero gives 0E0 with its sign.
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd32(const float val) noexcept -> floating_decimal<std::uint32_t>
{
    const auto bits {bit_cast<std::uint32_t>(val)};
    const auto significand {bits & ((UINT32_C(1) << binary32_format::significand_bits) - 1U)};
    const auto exponent {static_cast<int>((bits >> binary32_format::significand_bits) & 0xFFU)};
    const auto sign {(bits >> 31U) != 0U};
    if (exponent == 0 && significand == 0U)
    {
        return {0U, 0, sign};
    }
    return to_decimal<binary32_format>(significand, exponent, sign);
}

// The same for double.
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd64(const double val) noexcept -> floating_decimal<std::uint64_t>
{
    const auto bits {bit_cast<std::uint64_t>(val)};
    const auto significand {bits & ((UINT64_C(1) << binary64_format::significand_bits) - 1U)};
    const auto exponent {static_cast<int>((bits >> binary64_format::significand_bits) & 0x7FFU)};
    const auto sign {(bits >> 63U) != 0U};
    if (exponent == 0 && significand == 0U)
    {
        return {0U, 0, sign};
    }
    return to_decimal<binary64_format>(significand, exponent, sign);
}

} // namespace dragonbox

// Dragonbox for float and double, Ryu for the wider types Dragonbox does not cover.
template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd(const T val) noexcept -> ryu::floating_decimal_128
{
    return ryu::floating_point_to_fd128(val);
}

BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd(const float val) noexcept
    -> dragonbox::floating_decimal<std::uint32_t>
{
    return dragonbox::floating_point_to_fd32(val);
}

BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd(const double val) noexcept
    -> dragonbox::floating_decimal<std::uint64_t>
{
    return dragonbox::floating_point_to_fd64(val);
}

#if BOOST_DECIMAL_LDBL_BITS == 64

BOOST_DECIMAL_CXX20_CONSTEXPR auto floating_point_to_fd(const long double val) noexcept
    -> dragonbox::floating_decimal<std::uint64_t>
{
    return dragonbox::floating_point_to_fd64(static_cast<double>(val));
}

#endif

} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_DRAGONBOX_DRAGONBOX_HPP
