// Copyright 2020-2023 Junekey Jeon
// Copyright 2022 Peter Dimov
// Copyright 2023-2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_TO_CHARS_INTEGER_IMPL_HPP
#define BOOST_TO_CHARS_INTEGER_IMPL_HPP

#include <boost/decimal/detail/apply_sign.hpp>
#include "int128.hpp"
#include <boost/decimal/detail/to_chars_result.hpp>
#include <boost/decimal/detail/memcpy.hpp>
#include <boost/decimal/detail/countl.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <limits>
#include <system_error>
#include <type_traits>
#include <cstdint>
#include <climits>
#endif

namespace boost {
namespace decimal {
namespace detail {

#if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE char digit_table[] = {
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
        'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
        'u', 'v', 'w', 'x', 'y', 'z'
};
#endif

#if defined(__GNUC__) && __GNUC__ == 7
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

// Use a simple lookup table to put together the Integer in character form
template <typename Integer, typename Unsigned_Integer>
BOOST_DECIMAL_CUDA_CONSTEXPR auto to_chars_integer_impl(char* first, char* last, Integer value, int) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_integral_v, Integer, detail::is_integral_v, Unsigned_Integer, to_chars_result)
{
    #if defined(__CUDACC__) && defined(BOOST_DECIMAL_ENABLE_CUDA)
    constexpr char digit_table[] = {
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
        'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
        'u', 'v', 'w', 'x', 'y', 'z'
    };
    #endif

    const std::ptrdiff_t output_length = last - first;

    Unsigned_Integer unsigned_value {};

    BOOST_DECIMAL_IF_CONSTEXPR (std::is_signed<Integer>::value)
    {
        if (value < static_cast<Integer>(0))
        {
            *first++ = '-';
            unsigned_value = static_cast<Unsigned_Integer>(detail::apply_sign(value));
        }
        else
        {
            unsigned_value = static_cast<Unsigned_Integer>(value);
        }
    }
    else
    {
        unsigned_value = static_cast<Unsigned_Integer>(value);
    }

    // The bit length, rounded up to whole 4-bit digits, gives the number of digits, which go straight
    // to the output from LSB to MSB
    constexpr int value_bits {static_cast<int>(sizeof(Unsigned_Integer) * CHAR_BIT)};
    const std::ptrdiff_t num_chars {(value_bits - countl_zero(unsigned_value) + 3) / 4};

    if (num_chars > output_length)
    {
        return {last, std::errc::value_too_large};
    }

    auto* cursor {first + num_chars};
    while (unsigned_value != 0U)
    {
        *--cursor = digit_table[static_cast<std::size_t>(unsigned_value & 15U)]; // 1<<4 - 1
        unsigned_value >>= 4U;
    }

    return {first + num_chars, std::errc()};
}

// Specialization for base-10

#if !defined(__CUDACC__) || !defined(BOOST_DECIMAL_ENABLE_CUDA)
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE char radix_table[] = {
    '0', '0', '0', '1', '0', '2', '0', '3', '0', '4',
    '0', '5', '0', '6', '0', '7', '0', '8', '0', '9',
    '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
    '1', '5', '1', '6', '1', '7', '1', '8', '1', '9',
    '2', '0', '2', '1', '2', '2', '2', '3', '2', '4',
    '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
    '3', '0', '3', '1', '3', '2', '3', '3', '3', '4',
    '3', '5', '3', '6', '3', '7', '3', '8', '3', '9',
    '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
    '4', '5', '4', '6', '4', '7', '4', '8', '4', '9',
    '5', '0', '5', '1', '5', '2', '5', '3', '5', '4',
    '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
    '6', '0', '6', '1', '6', '2', '6', '3', '6', '4',
    '6', '5', '6', '6', '6', '7', '6', '8', '6', '9',
    '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
    '7', '5', '7', '6', '7', '7', '7', '8', '7', '9',
    '8', '0', '8', '1', '8', '2', '8', '3', '8', '4',
    '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
    '9', '0', '9', '1', '9', '2', '9', '3', '9', '4',
    '9', '5', '9', '6', '9', '7', '9', '8', '9', '9'
};
#endif

// Writes the two digits of pair < 100. CUDA device code can not read the table at namespace scope,
// so it has a local copy
BOOST_DECIMAL_CUDA_CONSTEXPR void write_2_digits(char* dest, const std::uint32_t pair) noexcept
{
    #if defined(__CUDACC__) && defined(BOOST_DECIMAL_ENABLE_CUDA)
    constexpr char radix_table[] = {
        '0', '0', '0', '1', '0', '2', '0', '3', '0', '4',
        '0', '5', '0', '6', '0', '7', '0', '8', '0', '9',
        '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
        '1', '5', '1', '6', '1', '7', '1', '8', '1', '9',
        '2', '0', '2', '1', '2', '2', '2', '3', '2', '4',
        '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
        '3', '0', '3', '1', '3', '2', '3', '3', '3', '4',
        '3', '5', '3', '6', '3', '7', '3', '8', '3', '9',
        '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
        '4', '5', '4', '6', '4', '7', '4', '8', '4', '9',
        '5', '0', '5', '1', '5', '2', '5', '3', '5', '4',
        '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
        '6', '0', '6', '1', '6', '2', '6', '3', '6', '4',
        '6', '5', '6', '6', '6', '7', '6', '8', '6', '9',
        '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
        '7', '5', '7', '6', '7', '7', '7', '8', '7', '9',
        '8', '0', '8', '1', '8', '2', '8', '3', '8', '4',
        '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
        '9', '0', '9', '1', '9', '2', '9', '3', '9', '4',
        '9', '5', '9', '6', '9', '7', '9', '8', '9', '9'
    };
    #endif

    boost::decimal::detail::memcpy(dest, radix_table + pair * 2U, 2);
}

// See: https://jk-jeon.github.io/posts/2022/02/jeaiii-algorithm/
// https://arxiv.org/abs/2101.11408
BOOST_DECIMAL_CUDA_CONSTEXPR char* decompose32(std::uint32_t value, char* buffer) noexcept
{
    constexpr auto mask = (UINT64_C(1) << 57) - 1;
    auto y = value * UINT64_C(1441151881);

    for (std::size_t i {}; i < 10; i += 2)
    {
        write_2_digits(buffer + i, static_cast<std::uint32_t>(y >> 57));
        y &= mask;
        y *= 100U;
    }

    return buffer + 10;
}

constexpr std::uint32_t ten_8_u32 {100000000U};

// Writes the eight digits of group < 10^8, with the leading zeros
BOOST_DECIMAL_CUDA_CONSTEXPR void write_8_digits(char* dest, const std::uint32_t group) noexcept
{
    const auto high {group / 10000U};
    const auto low {group % 10000U};
    write_2_digits(dest, high / 100U);
    write_2_digits(dest + 2, high % 100U);
    write_2_digits(dest + 4, low / 100U);
    write_2_digits(dest + 6, low % 100U);
}

template <typename Unsigned_Integer>
using to_chars_work_t = std::conditional_t<(sizeof(Unsigned_Integer) > sizeof(std::uint32_t)), std::uint64_t, std::uint32_t>;

// Writes groups of eight digits before cursor until the value fits in 32 bits
BOOST_DECIMAL_CUDA_CONSTEXPR auto narrow_to_32_bits(std::uint64_t value, char*& cursor) noexcept -> std::uint32_t
{
    while (value > UINT32_MAX)
    {
        cursor -= 8;
        write_8_digits(cursor, static_cast<std::uint32_t>(value % ten_8_u32));
        value /= ten_8_u32;
    }
    return static_cast<std::uint32_t>(value);
}

BOOST_DECIMAL_CUDA_CONSTEXPR auto narrow_to_32_bits(const std::uint32_t value, char*&) noexcept -> std::uint32_t
{
    return value;
}

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable: 4127 4146)
#endif

template <typename Integer, std::enable_if_t<sizeof(Integer) != 16, bool> = true>
BOOST_DECIMAL_CUDA_CONSTEXPR to_chars_result to_chars_integer_impl(char* first, char* last, Integer value) noexcept
{
    using Unsigned_Integer = make_unsigned_t<Integer>;
    Unsigned_Integer unsigned_value {};

    bool is_negative = false;

    // Strip the sign from the value and apply at the end after parsing if the type is signed
    BOOST_DECIMAL_IF_CONSTEXPR (is_signed_v<Integer>)
    {
        if (value < static_cast<Integer>(0))
        {
            is_negative = true;
            unsigned_value = apply_sign(value);
        }
        else
        {
            unsigned_value = static_cast<Unsigned_Integer>(value);
        }
    }
    else
    {
        unsigned_value = static_cast<Unsigned_Integer>(value);
    }

    const std::ptrdiff_t user_buffer_size = last - first - static_cast<std::ptrdiff_t>(is_negative);

    const auto converted_value_digits {num_digits(unsigned_value)};

    if (converted_value_digits > user_buffer_size)
    {
        return {last, std::errc::value_too_large};
    }

    if (is_negative)
    {
        *first++ = '-';
    }

    // From the end: groups of eight digits while the value is longer, then two digits at a time.
    // A value that fits in 32 bits continues in 32 bits, so a 32-bit target does less 64-bit division.
    auto cursor {first + converted_value_digits};
    auto rest {narrow_to_32_bits(static_cast<to_chars_work_t<Unsigned_Integer>>(unsigned_value), cursor)};
    if (rest >= ten_8_u32)
    {
        cursor -= 8;
        write_8_digits(cursor, rest % ten_8_u32);
        rest /= ten_8_u32;
    }

    while (rest >= 100U)
    {
        cursor -= 2;
        write_2_digits(cursor, rest % 100U);
        rest /= 100U;
    }

    if (rest >= 10U)
    {
        write_2_digits(first, rest);
    }
    else
    {
        *first = static_cast<char>('0' + rest);
    }

    return {first + converted_value_digits, std::errc()};
}

// Writes the 19 digits of a value below 10^19, with the leading zeros
BOOST_DECIMAL_CUDA_CONSTEXPR char* to_chars_19_digits(char* first, const std::uint64_t value) noexcept
{
    constexpr std::uint64_t ten_9 {UINT64_C(1000000000)};
    const auto rest {value % (ten_9 * ten_9)};
    *first = static_cast<char>('0' + value / (ten_9 * ten_9));

    char buffer[10] {};
    decompose32(static_cast<std::uint32_t>(rest / ten_9), buffer); // Always returns 10 digits (to include leading 0s) which we want
    boost::decimal::detail::memcpy(first + 1, buffer + 1, 9);
    decompose32(static_cast<std::uint32_t>(rest % ten_9), buffer);
    boost::decimal::detail::memcpy(first + 10, buffer + 1, 9);

    return first + 19;
}

template <typename Integer, std::enable_if_t<sizeof(Integer) == 16, bool> = true>
BOOST_DECIMAL_CUDA_CONSTEXPR to_chars_result to_chars_integer_impl(char* first, char* last, Integer value) noexcept
{
    using Unsigned_Integer = boost::int128::uint128_t;

    Unsigned_Integer unsigned_value {};

    BOOST_DECIMAL_ATTRIBUTE_UNUSED bool is_negative = false;

    // Strip the sign from the value, and write it first if the type is signed
    BOOST_DECIMAL_IF_CONSTEXPR (std::numeric_limits<Integer>::is_signed
                        #ifdef BOOST_DECIMAL_HAS_INT128
                        || std::is_same<boost::decimal::detail::builtin_int128_t, Integer>::value
                        #endif
                        )
    {
        if (value < 0U)
        {
            is_negative = true;
            unsigned_value = -(static_cast<Unsigned_Integer>(value));
        }
        else
        {
            unsigned_value = static_cast<Unsigned_Integer>(value);
        }
    }
    else
    {
        unsigned_value = static_cast<Unsigned_Integer>(value);
    }

    auto converted_value = unsigned_value;

    if (is_negative)
    {
        if (first == last)
        {
            return {last, std::errc::value_too_large};
        }
        *first++ = '-';
    }

    // A magnitude that fits 64 bits goes to the 64-bit writer, which checks the buffer size
    if (converted_value.high == 0U)
    {
        return to_chars_integer_impl(first, last, static_cast<std::uint64_t>(converted_value));
    }

    const int converted_value_digits = num_digits(converted_value);

    if (converted_value_digits > last - first)
    {
        return {last, std::errc::value_too_large};
    }

    // One division by 10^19 leaves the last 19 digits in 64 bits. The digits in front are below
    // 2^128 / 10^19, so they are a value below 10^19, or one digit and 19 more digits
    const auto dr {impl::divmod_pow10_uint128(converted_value, 19)};
    if (dr.quotient < UINT64_C(10000000000000000000))
    {
        first = to_chars_integer_impl(first, last, static_cast<std::uint64_t>(dr.quotient)).ptr;
    }
    else
    {
        const auto top {impl::divmod_pow10_uint128(dr.quotient, 19)};
        *first++ = static_cast<char>('0' + static_cast<int>(top.quotient));
        first = to_chars_19_digits(first, static_cast<std::uint64_t>(top.remainder));
    }

    return {to_chars_19_digits(first, static_cast<std::uint64_t>(dr.remainder)), std::errc()};
}

#ifdef _MSC_VER
# pragma warning(pop)
#endif

#if defined(__GNUC__) && __GNUC__ == 7
#pragma GCC diagnostic pop
#endif

} // namespace detail
} // namespace decimal
} // namespace boost

#endif //BOOST_TO_CHARS_INTEGER_IMPL_HPP
