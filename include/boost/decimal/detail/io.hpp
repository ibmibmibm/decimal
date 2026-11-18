// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_IO_HPP
#define BOOST_DECIMAL_DETAIL_IO_HPP

#include <boost/decimal/fwd.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/parser.hpp>
#include <boost/decimal/detail/attributes.hpp>
#include <boost/decimal/detail/fenv_rounding.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/locale_conversion.hpp>
#include <boost/decimal/charconv.hpp>

#if !defined(BOOST_DECIMAL_DISABLE_CLIB)

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cerrno>
#include <cstring>
#include <limits>
#include <ios>
#include <iostream>
#include <locale>
#include <system_error>
#include <type_traits>
#include <string>
#include <memory>
#endif

namespace boost {
namespace decimal {

namespace detail {

// The float field picks the format as for a double, so hexfloat is fixed | scientific
inline auto float_field_format(const std::ios_base::fmtflags flags) noexcept -> chars_format
{
    const auto float_field {flags & std::ios_base::floatfield};
    return float_field == (std::ios_base::fixed | std::ios_base::scientific) ? chars_format::hex :
           float_field == std::ios_base::scientific ? chars_format::scientific :
           float_field == std::ios_base::fixed ? chars_format::fixed : chars_format::general;
}

// Reads the characters of a number as num_get does, without a putback of what follows them, and sets
// eof at the end of the input. Returns the thousands separator if the characters have one, else '\0'.
// If the text has a separator or a point other than '.', converted gets the text in the C locale
template <typename charT, typename traits>
auto read_number(std::basic_istream<charT, traits>& is, std::string& text, std::string& converted, const chars_format fmt, bool& eof, bool& grouping_ok) -> char
{
    const auto& np = std::use_facet<std::numpunct<charT>>(is.getloc());
    const auto grouping {np.grouping()};
    const auto point {is.narrow(np.decimal_point(), '\0')};
    const auto separator {!grouping.empty() && grouping[0] > 0 ? thousands_separator(is.narrow(np.thousands_sep(), '\0')) : '\0'};

    number_chars number {point, separator, fmt == chars_format::hex, fmt != chars_format::fixed};
    auto* sb {is.rdbuf()};
    // A chunk of characters goes into text at once, which is faster than one character at a time
    char chunk[32];
    std::size_t size {};
    auto grouped {false};
    for (auto ch {sb->sgetc()};; ch = sb->snextc())
    {
        if (traits::eq_int_type(ch, traits::eof()))
        {
            eof = true;
            break;
        }
        const auto c {traits::to_char_type(ch)};
        // A space ends the number, also when the locale uses it as the separator
        const auto narrow {std::is_same<charT, char>::value ? static_cast<char>(c) : is.narrow(c, '\0')};
        if (narrow == '\0' || narrow == ' ' || (narrow >= '\t' && narrow <= '\r') || !number.take(narrow))
        {
            break;
        }
        if (size == sizeof(chunk))
        {
            text.append(chunk, size);
            size = 0U;
        }
        chunk[size++] = narrow;
        grouped = grouped || narrow == separator;
    }
    text.append(chunk, size);

    if (grouped || point != '.')
    {
        converted = text;
        convert_string_to_c_locale(&converted[0], grouping, grouped ? separator : '\0', point, grouping_ok);
    }

    return grouped ? separator : '\0';
}

} // namespace detail

// 3.2.10 Formatted input:
BOOST_DECIMAL_EXPORT template <typename charT, typename traits, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
auto operator>>(std::basic_istream<charT, traits>& is, DecimalType& d)
    -> std::enable_if_t<detail::is_decimal_floating_point_v<DecimalType>, std::basic_istream<charT, traits>&>
{
    const typename std::basic_istream<charT, traits>::sentry guard {is};
    if (!guard)
    {
        errno = EINVAL;
        return is;
    }

    // The float field picks the format as operator<< does
    const auto fmt {detail::float_field_format(is.flags())};

    // The facet and the copy for the conversion can throw too
    std::string text;
    std::string converted;
    auto eof {false};
    bool grouping_ok {true};
    #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
    char separator {};
    try
    {
        separator = detail::read_number(is, text, converted, fmt, eof, grouping_ok);
    }
    catch (...)
    {
        // An exception of the stream buffer sets badbit and keeps d, as in operator<<
        try
        {
            is.setstate(std::ios_base::badbit);
        }
        catch (const std::ios_base::failure&)
        {
        }
        if (is.exceptions() & std::ios_base::badbit)
        {
            throw;
        }
        return is;
    }
    #else
    const auto separator {detail::read_number(is, text, converted, fmt, eof, grouping_ok)};
    #endif
    // Without a separator, a '.' point needs no conversion, so the text needs no copy
    const auto buffer {converted.empty() ? &text[0] : &converted[0]};

    // Having a leading + sign is legal in iostream, but not allowed with charconv
    const auto first {buffer + (*buffer == '+')};

    auto r = from_chars(first, buffer + std::strlen(buffer), d, fmt);

    if (BOOST_DECIMAL_UNLIKELY(r.ec == std::errc::not_supported))
    {
        d = std::numeric_limits<DecimalType>::signaling_NaN();
    }
    else if (static_cast<int>(r.ec) == EINVAL)
    {
        errno = EINVAL;
    }

    // Put back the characters after the number, at most a short tail such as "e+"
    const auto consumed {detail::unconverted_length(text.c_str(), r.ptr - buffer, separator)};
    for (auto i {text.size()}; i > consumed; --i)
    {
        is.putback(is.widen(text[i - 1U]));
    }

    // The state is set last, because putback has no effect on a failed stream and setstate can throw.
    // A grouping that num_put does not give keeps the value and fails, as in num_get
    std::ios_base::iostate state {eof && consumed == text.size() ? std::ios_base::eofbit : std::ios_base::goodbit};
    if (!r || !grouping_ok)
    {
        state |= std::ios_base::failbit;
    }
    is.setstate(state);

    return is;
}

namespace detail {

// showpoint: the general format keeps P significant digits with the trailing zeros, as printf %#g,
// and every format keeps the decimal point. hexfloat ignores the precision
template <BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
auto to_chars_showpoint(char* first, char* last, const DecimalType& d, const chars_format fmt, const int local_precision) -> to_chars_result
{
    to_chars_result r {};
    if (fmt == chars_format::general && isfinite(d))
    {
        // The exponent after the rounding to P digits decides between the fixed and the scientific format
        const auto digits {local_precision == 0 ? 1 : local_precision};
        r = to_chars(first, last, d, chars_format::scientific, digits - 1);
        if (!r)
        {
            return r;
        }

        // The exponent is e+XX or e-XX at the end
        const auto exp_mark {static_cast<const char*>(std::memchr(first, 'e', static_cast<std::size_t>(r.ptr - first)))};
        int exp {};
        for (auto c {exp_mark + 2}; c != r.ptr; ++c)
        {
            exp = exp * 10 + (*c - '0');
        }
        exp = exp_mark[1] == '-' ? -exp : exp;
        if (exp >= -4 && exp < digits)
        {
            r = to_chars(first, last, d, chars_format::fixed, digits - 1 - exp);
        }
    }
    else if (fmt == chars_format::hex)
    {
        r = to_chars(first, last, d, fmt);
    }
    else
    {
        r = to_chars(first, last, d, fmt, local_precision);
    }

    if (r && isfinite(d) && std::memchr(first, '.', static_cast<std::size_t>(r.ptr - first)) == nullptr)
    {
        if (r.ptr == last)
        {
            return {last, std::errc::value_too_large};
        }
        auto point {static_cast<char*>(std::memchr(first, fmt == chars_format::hex ? 'p' : 'e', static_cast<std::size_t>(r.ptr - first)))};
        point = point == nullptr ? r.ptr : point;
        std::memmove(point + 1, point, static_cast<std::size_t>(r.ptr - point));
        *point = '.';
        ++r.ptr;
    }

    return r;
}

// Writes the C locale characters [first, last) as num_put does, for all the decimal types
template <typename charT, typename traits>
void put_chars(std::basic_ostream<charT, traits>& os, char* first, char* const last)
{
    const auto flags {os.flags()};
    if (flags & std::ios_base::uppercase)
    {
        detail::make_uppercase(first, last);
    }

    const char* str {first};
    const bool minus {str != last && *str == '-'};
    const bool plus {(flags & std::ios_base::showpos) && !minus};
    str += static_cast<int>(minus);

    // The numpunct of the stream groups the digits before the point as num_put does:
    // the last group size repeats, and a size of 0 or CHAR_MAX stops the grouping
    const auto loc {os.getloc()};
    const auto& np {std::use_facet<std::numpunct<charT>>(loc)};
    const auto& ct {std::use_facet<std::ctype<charT>>(loc)};
    const auto grouping {np.grouping()};
    std::size_t digits {};
    while (str + digits != last && str[digits] >= '0' && str[digits] <= '9')
    {
        ++digits;
    }
    std::size_t first_group {};
    auto separators {detail::count_separators(grouping, digits, first_group)};

    // Pad to the width as num_put does: internal puts the fill after the sign
    const auto total {static_cast<std::size_t>(last - str) + static_cast<std::size_t>(minus || plus) + separators};
    // As num_put does, the width resets before the output, so it also resets when the stream buffer throws
    const auto width {static_cast<std::size_t>(os.width() > 0 ? os.width() : 0)};
    os.width(0);
    auto padding {width > total ? width - total : 0U};
    const auto adjust {flags & std::ios_base::adjustfield};

    auto& sb {*os.rdbuf()};
    bool failed {false};
    const auto put = [&](const charT c) { failed = failed || traits::eq_int_type(sb.sputc(c), traits::eof()); };
    const auto put_padding = [&]
    {
        for (; padding > 0U && !failed; --padding)
        {
            put(os.fill());
        }
    };

    // A run is widened and written in blocks, and the output stops at the first failed write
    const auto put_run = [&](std::size_t n)
    {
        constexpr std::size_t block_size {64U};
        charT block[block_size];
        while (n > 0U && !failed)
        {
            const auto k {n < block_size ? n : block_size};
            ct.widen(str, str + k, block);
            failed = failed || sb.sputn(block, static_cast<std::streamsize>(k)) != static_cast<std::streamsize>(k);
            str += k;
            n -= k;
        }
    };

    if (adjust != std::ios_base::left && adjust != std::ios_base::internal)
    {
        put_padding();
    }
    if ((minus || plus) && !failed)
    {
        put(ct.widen(minus ? '-' : '+'));
    }
    if (adjust == std::ios_base::internal)
    {
        put_padding();
    }
    put_run(first_group);
    while (separators > 0U && !failed)
    {
        put(np.thousands_sep());
        put_run(detail::group_size(grouping, --separators, digits));
    }
    const auto point {failed ? nullptr : static_cast<const char*>(std::memchr(str, '.', static_cast<std::size_t>(last - str)))};
    if (point != nullptr)
    {
        put_run(static_cast<std::size_t>(point - str));
        if (!failed)
        {
            put(np.decimal_point());
        }
        ++str;
    }
    put_run(static_cast<std::size_t>(last - str));
    put_padding();

    if (failed)
    {
        os.setstate(std::ios_base::badbit);
    }
}

// Writes d as num_put writes a double
template <typename charT, typename traits, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
void put_decimal(std::basic_ostream<charT, traits>& os, const DecimalType& d)
{
    const auto flags {os.flags()};
    const auto fmt {float_field_format(flags)};

    // A negative precision means 6, as for printf. hexfloat, inf and NaN do not use the precision,
    // so it is 0 for them and their short text needs no large buffer
    constexpr std::streamsize max_precision {std::numeric_limits<int>::max() / 2};
    const auto stream_precision {os.precision() < 0 ? 6 : os.precision() < max_precision ? os.precision() : max_precision};
    const auto local_precision {fmt == chars_format::hex || !isfinite(d) ? 0 : static_cast<int>(stream_precision)};

    const auto print = [&](char* first, char* last)
    {
        if (flags & std::ios_base::showpoint)
        {
            return detail::to_chars_showpoint(first, last, d, fmt, local_precision);
        }

        // hexfloat ignores the precision
        return fmt == chars_format::hex ? to_chars(first, last, d, fmt) : to_chars(first, last, d, fmt, local_precision);
    };

    // Most values fit the local buffer. A longer one gets the length that to_chars needs at most,
    // where showpoint keeps the trailing zeros of the general format, as the scientific format does
    char local_buffer[512];
    std::unique_ptr<char[]> heap_buffer;
    char* buffer {local_buffer};
    auto r {print(local_buffer, local_buffer + sizeof(local_buffer))};
    if (r.ec == std::errc::value_too_large)
    {
        const auto sizing_fmt {fmt == chars_format::general && (flags & std::ios_base::showpoint) ? chars_format::scientific : fmt};
        const auto length {detail::to_chars_length(d, sizing_fmt, local_precision)};
        heap_buffer.reset(new char[length]);
        buffer = heap_buffer.get();
        r = print(buffer, buffer + length);
    }

    // The buffer sizes make this fail never, and a failure writes nothing and sets badbit
    if (BOOST_DECIMAL_UNLIKELY(!r))
    {
        os.setstate(std::ios_base::badbit);
        return;
    }

    detail::put_chars(os, buffer, r.ptr);
}

} // namespace detail

// 3.2.11 Formatted output
BOOST_DECIMAL_EXPORT template <typename charT, typename traits, BOOST_DECIMAL_DECIMAL_FLOATING_TYPE DecimalType>
auto operator<<(std::basic_ostream<charT, traits>& os, const DecimalType& d)
    -> std::enable_if_t<detail::is_decimal_floating_point_v<DecimalType>, std::basic_ostream<charT, traits>&>
{
    const typename std::basic_ostream<charT, traits>::sentry guard(os);
    if (guard)
    {
        #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
        try
        {
            detail::put_decimal(os, d);
        }
        catch (...)
        {
            // An exception during the output sets badbit, and goes to the caller only when badbit is in the
            // exception mask. setstate throws ios_base::failure for that mask, so the original exception is rethrown
            try
            {
                os.setstate(std::ios_base::badbit);
            }
            catch (const std::ios_base::failure&)
            {
            }
            if (os.exceptions() & std::ios_base::badbit)
            {
                throw;
            }
        }
        #else
        detail::put_decimal(os, d);
        #endif
    }

    return os;
}

} //namespace decimal
} //namespace boost

#endif // BOOST_DECIMAL_DISABLE_CLIB

#endif //BOOST_DECIMAL_DETAIL_IO_HPP
