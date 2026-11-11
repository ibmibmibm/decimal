// Copyright 2023 - 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_CSTDLIB_HPP
#define BOOST_DECIMAL_CSTDLIB_HPP

#include <boost/decimal/decimal32_t.hpp>
#include <boost/decimal/decimal64_t.hpp>
#include <boost/decimal/decimal128_t.hpp>
#include <boost/decimal/decimal_fast32_t.hpp>
#include <boost/decimal/decimal_fast64_t.hpp>
#include <boost/decimal/decimal_fast128_t.hpp>
#include <boost/decimal/charconv.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/utilities.hpp>
#include <boost/decimal/detail/locale_conversion.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <memory>
#include <new>
#include <limits>
#include <locale>
#include <cstdint>
#include <clocale>
#include <cstring>
#endif

#if !defined(BOOST_DECIMAL_DISABLE_CLIB)

namespace boost {
namespace decimal {

namespace detail {

// 3.8.2
template <typename TargetDecimalType, typename charT>
inline auto strtod_calculation(const charT* str, charT** endptr, const char* number, const std::size_t str_length, const char separator) noexcept -> TargetDecimalType
{
    // Strtod allows leading whitespace characters, as isspace in the C locale
    auto first {number};
    std::size_t first_pos {0};
    while (first_pos < str_length && (*first == ' ' || (*first >= '\t' && *first <= '\r')))
    {
        ++first;
        ++first_pos;
    }

    // The plus sign is also allowed before the value, but not before a minus sign
    if (first_pos < str_length && *first == '+' && *(first + 1) != '-')
    {
        ++first;
    }

    TargetDecimalType d {};
    const auto r {from_chars_general_impl<TargetDecimalType, true>(first, number + str_length, d, chars_format::general)};
    const bool overflow {r.ec == std::errc::result_out_of_range};
    if (overflow)
    {
        errno = ERANGE;
    }
    else if (!r)
    {
        d = std::numeric_limits<TargetDecimalType>::signaling_NaN();
        errno = EINVAL;
    }

    if (endptr != nullptr)
    {
        // Without a conversion, endptr is str (C 7.24.1.5)
        *endptr = const_cast<charT*>(str + unconverted_length(str, r || overflow ? r.ptr - number : 0, separator));
    }

    return d;
}

// Copies the number as the C locale writes it: as char, without the separators between digits
// of the integer part, and with a '.' point after that part
template <typename TargetDecimalType, typename charT>
inline auto strtod_converted(const charT* str, charT** endptr, const std::size_t str_length, const std::size_t integer_end, const char separator, const char point) noexcept -> TargetDecimalType
{
    // A number that does not fit into the stack buffer gets a buffer from new(std::nothrow)
    char stack_buffer[1024U];
    std::unique_ptr<char[]> heap_buffer;
    if (str_length >= sizeof(stack_buffer))
    {
        heap_buffer.reset(new(std::nothrow) char[str_length + 1U]);
        if (heap_buffer == nullptr)
        {
            errno = ENOMEM;
            return std::numeric_limits<TargetDecimalType>::quiet_NaN();
        }
    }
    char* buffer {heap_buffer == nullptr ? stack_buffer : heap_buffer.get()};

    std::size_t length {str_length};
    if (separator == '\0')
    {
        for (std::size_t i {}; i < str_length; ++i)
        {
            buffer[i] = static_cast<char>(str[i]);
        }
    }
    else
    {
        length = 0U;
        for (std::size_t i {}; i < str_length; ++i)
        {
            const auto c {static_cast<char>(str[i])};
            buffer[length] = c;
            length += i >= integer_end || c != separator || i == 0U || !is_decimal_digit(static_cast<char>(str[i - 1U]));
        }
    }
    buffer[length] = '\0';

    // The separators are all before the point
    const auto point_pos {integer_end - (str_length - length)};
    if (point != '.')
    {
        // A '.' that stays is not the point of this locale, so it ends the number
        const auto dot {std::strchr(buffer, '.')};
        length = dot != nullptr ? static_cast<std::size_t>(dot - buffer) : length;
        buffer[length] = '\0';
    }
    if (point_pos < length && buffer[point_pos] == point)
    {
        buffer[point_pos] = '.';
    }

    return strtod_calculation<TargetDecimalType>(str, endptr, buffer, length, separator);
}

template <typename TargetDecimalType>
inline auto strtod_impl(const char* str, char** endptr) noexcept -> TargetDecimalType
{
    if (str == nullptr)
    {
        errno = EINVAL;
        return std::numeric_limits<TargetDecimalType>::quiet_NaN();
    }

    // Only the characters that can belong to the number are read, as a loop over a long text
    // with endptr would read the rest of the text again for each number
    const std::locale loc {};
    const std::numpunct<char>& np = std::use_facet<std::numpunct<char>>(loc);
    char separator {};
    std::size_t integer_end {};
    const auto str_length {number_length(str, np, separator, integer_end)};
    const auto point {np.decimal_point()};

    // Without a separator, a '.' point needs no conversion, so the number is read in place
    if (separator == '\0' && point == '.')
    {
        return strtod_calculation<TargetDecimalType>(str, endptr, str, str_length, separator);
    }

    return strtod_converted<TargetDecimalType>(str, endptr, str_length, integer_end, separator, point);
}

// 3.9.2
template <typename TargetDecimalType>
inline auto wcstod_impl(const wchar_t* str, wchar_t** endptr) noexcept -> TargetDecimalType
{
    if (str == nullptr)
    {
        errno = EINVAL;
        return std::numeric_limits<TargetDecimalType>::quiet_NaN();
    }

    const std::locale loc {};
    const std::numpunct<char>& np = std::use_facet<std::numpunct<char>>(loc);
    char separator {};
    std::size_t integer_end {};
    const auto str_length {number_length(str, np, separator, integer_end)};

    return strtod_converted<TargetDecimalType>(str, endptr, str_length, integer_end, separator, np.decimal_point());
}

} //namespace detail

BOOST_DECIMAL_EXPORT template <typename TargetDecimalType = decimal64_t>
inline auto strtod(const char* str, char** endptr) noexcept -> TargetDecimalType
{
    return detail::strtod_impl<TargetDecimalType>(str, endptr);
}

BOOST_DECIMAL_EXPORT template <typename TargetDecimalType = decimal64_t>
inline auto wcstod(const wchar_t* str, wchar_t** endptr) noexcept -> TargetDecimalType
{
    return detail::wcstod_impl<TargetDecimalType>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod32(const char* str, char** endptr) noexcept -> decimal32_t
{
    return detail::strtod_impl<decimal32_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod32(const wchar_t* str, wchar_t** endptr) noexcept -> decimal32_t
{
    return detail::wcstod_impl<decimal32_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod32f(const char* str, char** endptr) noexcept -> decimal_fast32_t
{
    return detail::strtod_impl<decimal_fast32_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod32f(const wchar_t* str, wchar_t** endptr) noexcept -> decimal_fast32_t
{
    return detail::wcstod_impl<decimal_fast32_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod64(const char* str, char** endptr) noexcept -> decimal64_t
{
    return detail::strtod_impl<decimal64_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod64(const wchar_t* str, wchar_t** endptr) noexcept -> decimal64_t
{
    return detail::wcstod_impl<decimal64_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod64f(const char* str, char** endptr) noexcept -> decimal_fast64_t
{
    return detail::strtod_impl<decimal_fast64_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod64f(const wchar_t* str, wchar_t** endptr) noexcept -> decimal_fast64_t
{
    return detail::wcstod_impl<decimal_fast64_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod128(const char* str, char** endptr) noexcept -> decimal128_t
{
    return detail::strtod_impl<decimal128_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod128(const wchar_t* str, wchar_t** endptr) noexcept -> decimal128_t
{
    return detail::wcstod_impl<decimal128_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto strtod128f(const char* str, char** endptr) noexcept -> decimal_fast128_t
{
    return detail::strtod_impl<decimal_fast128_t>(str, endptr);
}

BOOST_DECIMAL_EXPORT inline auto wcstod128f(const wchar_t* str, wchar_t** endptr) noexcept -> decimal_fast128_t
{
    return detail::wcstod_impl<decimal_fast128_t>(str, endptr);
}

} // namespace decimal
} // namespace boost

#endif

#endif // BOOST_DECIMAL_CSTDLIB_HPP
