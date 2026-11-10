// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_LOCALE_CONVERSION_HPP
#define BOOST_DECIMAL_DETAIL_LOCALE_CONVERSION_HPP

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <locale>
#include <climits>
#include <cstring>
#include <string>
#include <type_traits>
#endif

namespace boost {
namespace decimal {
namespace detail {

// A char numpunct can give the first byte of the UTF-8 U+202F, which is not a character alone,
// so it becomes a space
inline auto thousands_separator(const char separator) noexcept -> char
{
    return static_cast<unsigned char>(separator) == 0xE2U ? ' ' : separator;
}

// The grouping of num_put: group i counts from the right, the last size repeats, and a size of
// 0 or CHAR_MAX stops the grouping, so that the other digits form one group
inline auto group_size(const std::string& grouping, const std::size_t i, const std::size_t digits) noexcept -> std::size_t
{
    const auto g {grouping[i < grouping.size() ? i : grouping.size() - 1U]};
    return g > 0 && g != CHAR_MAX ? static_cast<std::size_t>(g) : digits;
}

// The number of separators between digits integer digits, and the size of the leftmost group
inline auto count_separators(const std::string& grouping, const std::size_t digits, std::size_t& first_group) noexcept -> std::size_t
{
    std::size_t separators {};
    first_group = digits;
    while (!grouping.empty() && first_group > group_size(grouping, separators, digits))
    {
        first_group -= group_size(grouping, separators, digits);
        ++separators;
    }

    return separators;
}

inline auto is_decimal_digit(const char c) noexcept -> bool
{
    return c >= '0' && c <= '9';
}

// Takes the characters of a number one at a time, so that the input ends with the number. A tail that the
// value does not use can stay in the input, such as "e+", the "," of "1," or a NaN payload without ")"
class number_chars
{
public:
    number_chars(const char point, const char separator, const bool hex, const bool exponent) noexcept
        : point_ {point}, separator_ {separator}, hex_ {hex}, exponent_allowed_ {exponent} {}

    auto take(const char c) noexcept -> bool
    {
        if (word_ == nullptr && is_decimal_digit(c))
        {
            digit_ = digit_ || !exponent_;
            last_ = c;
            return true;
        }

        const auto lower {static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c)};
        if (word_ != nullptr ? take_word(c, lower) : take_numeric(c, lower))
        {
            last_ = lower;
            return true;
        }

        return false;
    }

private:
    // inf, infinity, nan and snan, and a payload of digits or in parentheses after nan
    auto take_word(const char c, const char lower) noexcept -> bool
    {
        if (paren_)
        {
            return last_ != ')' && (is_decimal_digit(c) || (lower >= 'a' && lower <= 'z') || c == '_' || c == ')');
        }
        if (word_[word_size_] != '\0')
        {
            const auto ok {lower == word_[word_size_]};
            word_size_ += ok;
            return ok;
        }
        if (word_[0] == 'i')
        {
            return false;
        }
        paren_ = c == '(' && !digit_;
        digit_ = digit_ || is_decimal_digit(c);
        return paren_ || is_decimal_digit(c);
    }

    auto take_numeric(const char c, const char lower) noexcept -> bool
    {
        if (c == '+' || c == '-')
        {
            return last_ == '\0' || last_ == (hex_ ? 'p' : 'e');
        }
        if ((last_ == '\0' || last_ == '+' || last_ == '-') && !exponent_ && (lower == 'i' || lower == 'n' || lower == 's'))
        {
            word_ = lower == 'i' ? "infinity" : lower == 'n' ? "nan" : "snan";
            word_size_ = 1U;
            return true;
        }
        if (hex_ && !exponent_ && lower >= 'a' && lower <= 'f')
        {
            digit_ = true;
            return true;
        }
        if (c == point_ || (separator_ != '\0' && c == separator_))
        {
            const auto ok {!point_seen_ && !exponent_ && (c == point_ || is_decimal_digit(last_))};
            point_seen_ = point_seen_ || (ok && c == point_);
            return ok;
        }
        if (lower == (hex_ ? 'p' : 'e') && !exponent_ && digit_ && exponent_allowed_)
        {
            exponent_ = true;
            return true;
        }

        return false;
    }

    const char* word_ {};
    std::size_t word_size_ {};
    char point_;
    char separator_;
    bool hex_;
    bool exponent_allowed_;
    char last_ {};
    bool point_seen_ {};
    bool exponent_ {};
    bool digit_ {};
    bool paren_ {};
};

// Removes the separators removed_sep between the digits of the integer part, and sets grouping_ok
// to false if they are not where num_put puts them. A removed_sep of '\0' removes nothing
inline void convert_string_to_c_locale(char* buffer, const std::string& grouping, const char removed_sep, const char locale_decimal_point, bool& grouping_ok) noexcept
{
    grouping_ok = true;

    if (removed_sep != '\0')
    {
        auto start {buffer};
        while (*start == ' ' || (*start >= '\t' && *start <= '\r'))
        {
            ++start;
        }
        if (*start == '-' || *start == '+')
        {
            ++start;
        }

        auto int_end {start};
        std::size_t digits {};
        while (is_decimal_digit(*int_end) || (*int_end == removed_sep && int_end != start && is_decimal_digit(int_end[1])))
        {
            digits += *int_end++ != removed_sep;
        }

        // Check the groups from the right, as num_get does
        std::size_t first_group {};
        const auto expected {count_separators(grouping, digits, first_group)};
        std::size_t separators {};
        std::size_t group {};
        for (auto p {int_end}; p != start; --p)
        {
            if (p[-1] != removed_sep)
            {
                ++group;
                continue;
            }
            grouping_ok = grouping_ok && separators < expected && group == group_size(grouping, separators, digits);
            ++separators;
            group = 0;
        }
        grouping_ok = grouping_ok && (separators == 0U || (separators == expected && group == first_group));

        auto write {start};
        for (auto read {start}; read != int_end; ++read)
        {
            if (*read != removed_sep)
            {
                *write++ = *read;
            }
        }
        std::memmove(write, int_end, std::strlen(int_end) + 1U);
    }

    if (locale_decimal_point != '.')
    {
        // A '.' that stays is not the point of this locale, so it ends the number
        const auto dot {std::strchr(buffer, '.')};
        if (dot != nullptr)
        {
            *dot = '\0';
        }
        // A wide point without a char form gives '\0', and then the text has no point
        const auto p {locale_decimal_point != '\0' ? std::strchr(buffer, static_cast<int>(locale_decimal_point)) : nullptr};
        if (p != nullptr)
        {
            *p = '.';
        }
    }
}

// Returns the separator, or '\0' with no grouping
inline auto convert_string_to_c_locale(char* buffer) -> char
{
    const std::numpunct<char>& np = std::use_facet<std::numpunct<char>>(std::locale());
    const auto grouping {np.grouping()};
    const auto removed_sep {!grouping.empty() && grouping[0] > 0 ? thousands_separator(np.thousands_sep()) : '\0'};
    bool grouping_ok {};
    convert_string_to_c_locale(buffer, grouping, removed_sep, np.decimal_point(), grouping_ok);
    return removed_sep;
}

// The number of characters of str that give the first used characters of the converted buffer,
// where convert_string_to_c_locale removed each separator between digits of the integer part
template <typename charT>
auto unconverted_length(const charT* str, std::ptrdiff_t used, const char separator) noexcept -> std::size_t
{
    if (separator == '\0')
    {
        return static_cast<std::size_t>(used);
    }

    std::size_t start {};
    while (str[start] == ' ' || (str[start] >= '\t' && str[start] <= '\r'))
    {
        ++start;
    }
    start += str[start] == '-' || str[start] == '+';

    std::size_t length {};
    for (auto integer_part {true}; used > 0; ++length)
    {
        // A wide character that is not a char must not look like a digit or the separator
        const auto c {static_cast<std::make_unsigned_t<charT>>(str[length])};
        const auto digit {c >= '0' && c <= '9'};
        const auto removed {integer_part && c == static_cast<unsigned char>(separator) && length > start &&
                            str[length + 1U] >= '0' && str[length + 1U] <= '9'};
        integer_part = integer_part && (removed || length < start || digit);
        used -= !removed;
    }

    return length;
}

// The length of the characters of str that can belong to a number. Gives the separator if str has one
// between digits of the integer part, else '\0', and the end of the integer part
template <typename charT>
auto number_length(const charT* str, const std::numpunct<char>& np, char& separator, std::size_t& integer_end) noexcept -> std::size_t
{
    const auto grouping {np.grouping()};
    const auto sep {!grouping.empty() && grouping[0] > 0 ? static_cast<unsigned char>(thousands_separator(np.thousands_sep())) : 0U};
    const auto point {static_cast<unsigned char>(np.decimal_point())};
    separator = '\0';

    std::size_t i {};
    while (str[i] == ' ' || (str[i] >= '\t' && str[i] <= '\r'))
    {
        ++i;
    }
    const auto number_start {i};
    i += str[i] == '-' || str[i] == '+';

    const auto start {i};
    for (;; ++i)
    {
        const auto c {static_cast<std::make_unsigned_t<charT>>(str[i])};
        if (c >= '0' && c <= '9')
        {
            continue;
        }
        if (c != sep || c == 0U || i == start || str[i + 1] < '0' || str[i + 1] > '9')
        {
            break;
        }
        separator = static_cast<char>(sep);
    }
    integer_end = i;

    // The scan of the stream input ends at the end of the number, so a loop over a long text stays linear
    number_chars number {static_cast<char>(point), static_cast<char>(sep), false, true};
    for (i = number_start;; ++i)
    {
        const auto c {static_cast<std::make_unsigned_t<charT>>(str[i])};
        if (!number.take(c < 128U || c == point ? static_cast<char>(c) : '\0'))
        {
            return i;
        }
    }
}

inline int convert_pointer_pair_to_local_locale(char* first, const char* last, const std::locale& loc) noexcept
{
    const std::numpunct<char>& np = std::use_facet<std::numpunct<char>>(loc);

    const auto locale_decimal_point {np.decimal_point()};
    const auto locale_thousands_sep {thousands_separator(np.thousands_sep())};
    const auto grouping {np.grouping()};

    // Find the start of the number (skip sign if present)
    char* start = first;
    if (start < last && (*start == '-' || *start == '+'))
    {
        ++start;
    }

    // Find the actual end of the string
    auto string_end {start};
    while (string_end < last && *string_end != '\0')
    {
        ++string_end;
    }

    // A hex significand can have a letter before the point, so the point is found on its own
    bool hex {false};
    for (auto p {start}; p < string_end; ++p)
    {
        hex = hex || *p == 'p' || *p == 'P';
        if (*p == '.')
        {
            *p = locale_decimal_point;
        }
    }

    // The integer part is the digits after the sign, so an exponent, inf or nan has none after it.
    // A hex value has no grouping, as printf does
    auto int_end {start};
    while (!hex && int_end < string_end && *int_end >= '0' && *int_end <= '9')
    {
        ++int_end;
    }

    const auto digits {static_cast<std::size_t>(int_end - start)};
    std::size_t first_group {};
    const auto num_separators {locale_thousands_sep != '\0' ? count_separators(grouping, digits, first_group) : 0U};

    // If we need to add separators, shift content and insert them
    if (num_separators > 0U)
    {
        const auto new_length {static_cast<std::size_t>(string_end - first) + num_separators};

        // Check if we have enough space in the buffer
        if (new_length >= static_cast<std::size_t>(last - first))
        {
            // Not enough space, return error indicator
            return -1;
        }

        // Move the text after the integer part with its null terminator, and then the integer digits
        // from the right with a separator after each group. The offsets from first stay in the buffer
        auto old_pos {static_cast<std::size_t>(string_end - first) + 1U};
        auto new_pos {new_length + 1U};
        const auto int_end_pos {static_cast<std::size_t>(int_end - first)};
        while (old_pos > int_end_pos)
        {
            first[--new_pos] = first[--old_pos];
        }

        const auto int_begin_pos {static_cast<std::size_t>(start - first)};
        std::size_t separator {};
        std::size_t digits_in_group {};
        while (old_pos > int_begin_pos)
        {
            first[--new_pos] = first[--old_pos];

            // Insert a separator after each group from the right, but not after the leftmost digit
            if (++digits_in_group == group_size(grouping, separator, digits) && separator < num_separators)
            {
                first[--new_pos] = locale_thousands_sep;
                ++separator;
                digits_in_group = 0U;
            }
        }
    }

    return static_cast<int>(num_separators);
}

inline int convert_pointer_pair_to_local_locale(char* first, const char* last)
{
    const auto loc {std::locale()};
    return convert_pointer_pair_to_local_locale(first, last, loc);
}

} //namespace detail
} //namespace decimal
} //namespace boost

#endif //BOOST_DECIMAL_DETAIL_LOCALE_CONVERSION_HPP
