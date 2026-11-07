// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_FAST_FLOAT_BIG_UINT_HPP
#define BOOST_DECIMAL_DETAIL_FAST_FLOAT_BIG_UINT_HPP

#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/int128.hpp>
#include <boost/decimal/detail/power_tables.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstddef>
#include <cstdint>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace fast_float {

// A big integer for the exact compares of the conversions to binary, and for the exact digits of a double. It works on
// the words of the caller, one copy of the code for each capacity, and writes each word before it reads it.
class big_uint
{
    std::uint64_t* words_;
    std::size_t capacity_;
    std::size_t size_;

public:
    constexpr big_uint(std::uint64_t* words, const std::size_t capacity, const int128::uint128_t value) noexcept
        : words_ {words}, capacity_ {capacity}, size_ {value.high != 0U ? 2U : 1U}
    {
        BOOST_DECIMAL_ASSERT(capacity >= 2U);
        words_[0] = value.low;
        words_[1] = value.high;
    }

    // The number of words that the caller gave, which a debug build checks
    constexpr auto capacity() const noexcept -> std::size_t
    {
        return capacity_;
    }

    // A factor of zero would leave words of zero above the value, which compare() reads as a larger value
    BOOST_DECIMAL_CXX20_CONSTEXPR auto mul(const std::uint64_t factor) noexcept -> void
    {
        BOOST_DECIMAL_ASSERT(factor != 0U);
        std::uint64_t carry {};
        for (std::size_t i {}; i < size_; ++i)
        {
            const auto product {static_cast<int128::uint128_t>(words_[i]) * factor + int128::uint128_t {carry}};
            words_[i] = product.low;
            carry = product.high;
        }
        if (carry != 0U)
        {
            BOOST_DECIMAL_ASSERT(size_ < capacity());
            words_[size_++] = carry;
        }
    }

    BOOST_DECIMAL_CXX20_CONSTEXPR auto mul_pow5(std::int64_t n) noexcept -> void
    {
        constexpr std::uint64_t pow5_27 {UINT64_C(7450580596923828125)};
        for (; n >= 27; n -= 27)
        {
            mul(pow5_27);
        }
        const auto half {static_cast<std::uint64_t>(n / 2)};
        const auto other_half {static_cast<std::uint64_t>(n) - half};
        mul((pow10(half) >> half) * (pow10(other_half) >> other_half));
    }

    constexpr auto shl(const std::int64_t n) noexcept -> void
    {
        BOOST_DECIMAL_ASSERT(n >= 0);
        const auto word_shift {static_cast<std::size_t>(n / 64)};
        const auto bit_shift {static_cast<unsigned>(n % 64)};

        std::uint64_t top {};
        if (bit_shift != 0U)
        {
            top = words_[size_ - 1U] >> (64U - bit_shift);
            for (std::size_t i {size_ - 1U}; i > 0U; --i)
            {
                words_[i] = (words_[i] << bit_shift) | (words_[i - 1U] >> (64U - bit_shift));
            }
            words_[0] <<= bit_shift;
        }
        if (top != 0U)
        {
            BOOST_DECIMAL_ASSERT(size_ < capacity());
            words_[size_++] = top;
        }
        // Zero keeps one word, as compare() takes more words for a larger value
        if (word_shift != 0U && words_[size_ - 1U] != 0U)
        {
            BOOST_DECIMAL_ASSERT(size_ + word_shift <= capacity());
            for (std::size_t i {size_}; i > 0U; --i)
            {
                words_[i - 1U + word_shift] = words_[i - 1U];
            }
            for (std::size_t i {}; i < word_shift; ++i)
            {
                words_[i] = 0U;
            }
            size_ += word_shift;
        }
    }

    // Returns whether any of the bits shifted out was set
    constexpr auto shr(const std::int64_t n) noexcept -> bool
    {
        const auto word_shift {static_cast<std::size_t>(n / 64)};
        const auto bit_shift {static_cast<unsigned>(n % 64)};

        bool lost {};
        for (std::size_t i {}; i < word_shift && i < size_; ++i)
        {
            lost = lost || words_[i] != 0U;
            words_[i] = 0U;
        }
        if (word_shift >= size_)
        {
            size_ = 1U;
            return lost;
        }
        if (bit_shift != 0U)
        {
            lost = lost || (words_[word_shift] << (64U - bit_shift)) != 0U;
        }
        for (std::size_t i {word_shift}; i < size_; ++i)
        {
            const std::uint64_t next {i + 1U < size_ ? words_[i + 1U] : 0U};
            const auto shifted {bit_shift == 0U ? words_[i] : (words_[i] >> bit_shift) | (next << (64U - bit_shift))};
            words_[i] = 0U;
            words_[i - word_shift] = shifted;
        }
        size_ -= word_shift;
        trim();
        return lost;
    }

    // Returns the remainder
    BOOST_DECIMAL_CXX20_CONSTEXPR auto divide(const std::uint64_t divisor) noexcept -> std::uint64_t
    {
        std::uint64_t remainder {};
        for (std::size_t i {size_}; i > 0U; --i)
        {
            const int128::uint128_t current {remainder, words_[i - 1U]};
            const auto quotient {current / divisor};
            words_[i - 1U] = quotient.low;
            remainder = static_cast<std::uint64_t>(current - quotient * divisor);
        }
        trim();
        return remainder;
    }

    constexpr auto to_uint128() const noexcept -> int128::uint128_t
    {
        return {words_[1], words_[0]};
    }

    // The top word below size_ is not zero, unless the value is
    constexpr auto compare(const big_uint& other) const noexcept -> int
    {
        if (size_ != other.size_)
        {
            return size_ < other.size_ ? -1 : 1;
        }
        for (std::size_t i {size_}; i > 0U; --i)
        {
            if (words_[i - 1U] != other.words_[i - 1U])
            {
                return words_[i - 1U] < other.words_[i - 1U] ? -1 : 1;
            }
        }
        return 0;
    }

private:
    constexpr auto trim() noexcept -> void
    {
        while (size_ > 1U && words_[size_ - 1U] == 0U)
        {
            --size_;
        }
    }
};

// The sign of sig * 10^q - m * 2^e, exactly, in two big integers of capacity words each
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CXX20_CONSTEXPR auto compare_exact_words(std::uint64_t* value_words,
                                                                               std::uint64_t* other_words,
                                                                               const std::size_t capacity,
                                                                               const std::int64_t q,
                                                                               const int128::uint128_t sig,
                                                                               const int128::uint128_t m,
                                                                               const std::int64_t e) noexcept -> int
{
    big_uint value {value_words, capacity, sig};
    big_uint other {other_words, capacity, m};
    if (q >= 0)
    {
        value.mul_pow5(q);
    }
    else
    {
        other.mul_pow5(-q);
    }

    const auto two_shift {q - e};
    if (two_shift >= 0)
    {
        value.shl(two_shift);
    }
    else
    {
        other.shl(-two_shift);
    }

    return value.compare(other);
}

// The same with the words on the stack. It is not inline, because it runs only next to a boundary.
template <std::size_t capacity = 24>
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CXX20_CONSTEXPR auto compare_exact(const std::int64_t q,
                                                                         const int128::uint128_t sig,
                                                                         const int128::uint128_t m,
                                                                         const std::int64_t e) noexcept -> int
{
    std::uint64_t value_words[capacity];
    std::uint64_t other_words[capacity];
    return compare_exact_words(value_words, other_words, capacity, q, sig, m, e);
}

} // namespace fast_float
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_FAST_FLOAT_BIG_UINT_HPP
