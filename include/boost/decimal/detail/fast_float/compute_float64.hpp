// Copyright 2020-2023 Daniel Lemire
// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT64_HPP
#define BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT64_HPP

#include <boost/decimal/detail/fast_float/significand_tables.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/bit_cast.hpp>
#include <boost/decimal/detail/countl.hpp>
#include "../int128.hpp"
#include <boost/decimal/detail/integer_search_trees.hpp>
#include <boost/decimal/detail/power_tables.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <cstdint>
#include <cfloat>
#include <type_traits>
#include <cstring>
#include <cmath>
#endif

namespace boost {
namespace decimal {
namespace detail {
namespace fast_float {

constexpr auto umul128(const std::uint64_t x, const std::uint64_t y) noexcept
{
    return static_cast<int128::uint128_t>(x) * y;
}

// The binary32 and binary64 formats
template <typename T>
struct binary_format;

template <>
struct binary_format<double>
{
    using bits_type = std::uint64_t;
    static constexpr int mantissa_explicit_bits {52};
    static constexpr int minimum_exponent {-1023};
    static constexpr int infinite_power {0x7FF};
    static constexpr int min_exponent_round_to_even {-4};
    static constexpr int max_exponent_round_to_even {23};
    static constexpr int smallest_power_of_ten {-342};
    static constexpr int largest_power_of_ten {308};
};

template <>
struct binary_format<float>
{
    using bits_type = std::uint32_t;
    static constexpr int mantissa_explicit_bits {23};
    static constexpr int minimum_exponent {-127};
    static constexpr int infinite_power {0xFF};
    static constexpr int min_exponent_round_to_even {-17};
    static constexpr int max_exponent_round_to_even {10};
    static constexpr int smallest_power_of_ten {-65};
    static constexpr int largest_power_of_ten {38};
};

// The first row of significand_64 and significand_128
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE int smallest_power_of_five {-342};

// The explicit bits of the significand and the biased exponent of a binary value
struct adjusted_mantissa
{
    std::uint64_t mantissa;
    std::int32_t power2;

    constexpr auto operator==(const adjusted_mantissa& other) const noexcept -> bool
    {
        return mantissa == other.mantissa && power2 == other.power2;
    }
};

// The high 128 bits of w * 5^q, exact in the top bit_precision bits. Forced inline, as Clang calls it otherwise.
template <int bit_precision>
BOOST_DECIMAL_FORCE_INLINE constexpr auto compute_product_approximation(const std::int64_t q, const std::uint64_t w) noexcept -> int128::uint128_t
{
    const auto index {static_cast<std::size_t>(q - smallest_power_of_five)};
    auto first_product {umul128(w, significand_64[index])};

    constexpr std::uint64_t precision_mask {UINT64_MAX >> bit_precision};
    if ((first_product.high & precision_mask) == precision_mask)
    {
        const auto second_product {umul128(w, significand_128[index])};
        first_product += int128::uint128_t {second_product.high};
    }

    return first_product;
}

// floor(log2(5^q)) + q + 63 for q in (-400, 350)
constexpr auto binary_power(const std::int32_t q) noexcept -> std::int32_t
{
    return (((152170 + 65536) * q) >> 16) + 63;
}

// w * 2^shift * 10^q, round half to even, by the 128-bit product of Mushtak and Lemire. Forced inline for Clang. If
// wide, a non-null away tells if the bits below the rounding bit are at least 3 units from a boundary.
template <typename T, bool wide = false>
BOOST_DECIMAL_FORCE_INLINE constexpr auto compute_float(const std::int64_t q, std::uint64_t w, const std::int32_t shift = 0, bool* const away = nullptr) noexcept -> adjusted_mantissa
{
    using format = binary_format<T>;

    // Only the wide form reads shift and away, which GCC 9 would report as unused otherwise
    static_cast<void>(shift);
    static_cast<void>(away);

    // The wide form gives a float far below its range from the table too, as w * 2^shift can have 128 bits
    adjusted_mantissa answer {};
    if (w == 0U || q < (wide ? smallest_power_of_five : format::smallest_power_of_ten))
    {
        return answer;
    }
    if (q > format::largest_power_of_ten)
    {
        answer.power2 = format::infinite_power;
        return answer;
    }

    #ifdef BOOST_DECIMAL_HAS_STDBIT
    const int lz {std::countl_zero(w)};
    #else
    const int lz {countl_zero(w)};
    #endif
    w <<= lz;

    // The significand needs the implicit bit, a bit to round, and one bit that the shift by upperbit can drop
    const auto product {compute_product_approximation<format::mantissa_explicit_bits + 3>(q, w)};
    const auto upperbit {static_cast<int>(product.high >> 63U)};
    const auto drop {upperbit + 64 - format::mantissa_explicit_bits - 3};
    answer.mantissa = product.high >> drop;
    answer.power2 = binary_power(static_cast<std::int32_t>(q)) + upperbit - lz - format::minimum_exponent;
    BOOST_DECIMAL_IF_CONSTEXPR (wide)
    {
        answer.power2 += shift;
        const auto dropped {product.high & ((UINT64_C(1) << drop) - 1U)};
        if (away != nullptr)
        {
            *away = dropped >= 1U && dropped <= (UINT64_C(1) << drop) - 4U;
        }
    }

    if (answer.power2 <= 0)
    {
        // A subnormal value, or zero
        if (-answer.power2 + 1 >= 64)
        {
            answer.power2 = 0;
            answer.mantissa = 0;
            return answer;
        }

        answer.mantissa >>= -answer.power2 + 1;
        answer.mantissa += (answer.mantissa & 1U);
        answer.mantissa >>= 1U;
        answer.power2 = answer.mantissa < (UINT64_C(1) << format::mantissa_explicit_bits) ? 0 : 1;
        return answer;
    }

    // A tie is possible only when 5^q fits in a word. Then round down to the even value.
    if (product.low <= 1U && q >= format::min_exponent_round_to_even && q <= format::max_exponent_round_to_even &&
        (answer.mantissa & 3U) == 1U && (answer.mantissa << drop) == product.high)
    {
        answer.mantissa &= ~UINT64_C(1);
    }

    answer.mantissa += (answer.mantissa & 1U);
    answer.mantissa >>= 1U;
    if (answer.mantissa >= (UINT64_C(2) << format::mantissa_explicit_bits))
    {
        answer.mantissa = (UINT64_C(1) << format::mantissa_explicit_bits);
        ++answer.power2;
    }

    answer.mantissa &= ~(UINT64_C(1) << format::mantissa_explicit_bits);
    if (answer.power2 >= format::infinite_power)
    {
        answer.power2 = format::infinite_power;
        answer.mantissa = 0;
    }

    return answer;
}

// A big integer for the exact compares of the conversions to binary, and for the exact digits of a double
template <std::size_t capacity = 24>
class big_uint
{
    static_assert(capacity >= 2U, "the constructor writes two words");

    std::uint64_t words_[capacity] {};
    std::size_t size_ {};

public:
    constexpr explicit big_uint(const int128::uint128_t value) noexcept
    {
        words_[0] = value.low;
        words_[1] = value.high;
        size_ = value.high != 0U ? 2U : 1U;
    }

    BOOST_DECIMAL_CXX20_CONSTEXPR auto mul(const std::uint64_t factor) noexcept -> void
    {
        std::uint64_t carry {};
        for (std::size_t i {}; i < size_; ++i)
        {
            const auto product {umul128(words_[i], factor) + int128::uint128_t {carry}};
            words_[i] = product.low;
            carry = product.high;
        }
        if (carry != 0U)
        {
            BOOST_DECIMAL_ASSERT(size_ < capacity);
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
            BOOST_DECIMAL_ASSERT(size_ < capacity);
            words_[size_++] = top;
        }
        if (word_shift != 0U)
        {
            BOOST_DECIMAL_ASSERT(size_ + word_shift <= capacity);
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

// The sign of sig * 10^q - m * 2^e, exactly. It is not inline, because it runs only next to a boundary.
template <std::size_t capacity = 24>
BOOST_DECIMAL_NO_INLINE BOOST_DECIMAL_CXX20_CONSTEXPR auto compare_exact(const std::int64_t q,
                                                                         const int128::uint128_t sig,
                                                                         const int128::uint128_t m,
                                                                         const std::int64_t e) noexcept -> int
{
    big_uint<capacity> value {sig};
    big_uint<capacity> other {m};
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

template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compute_float_bounds(const std::int64_t q, const std::uint64_t w, const std::int32_t shift) noexcept -> adjusted_mantissa
{
    return compute_float<T, true>(q, w, shift);
}

// sig * 10^q for 65 to 128 bits. The leading 64 bits, or below the power table the 19 leading digits, give two bounds.
// If both round to the same value, it is the result. Else a halfway point is between, and an exact compare decides.
template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto compute_float_wide(const std::int64_t q, const int128::uint128_t sig) noexcept -> adjusted_mantissa
{
    using format = binary_format<T>;

    std::int64_t k {};
    std::int32_t s {};
    std::uint64_t power {1U};
    std::uint64_t w {};
    if (q >= smallest_power_of_five)
    {
        s = 64 - countl_zero(sig.high);
        w = static_cast<std::uint64_t>(sig >> static_cast<unsigned>(s));
    }
    else
    {
        k = num_digits(sig) - 19;
        power = pow10(static_cast<std::uint64_t>(k));
        w = static_cast<std::uint64_t>(sig / power);
    }

    // With the top bit of w set, w + 1 changes the product by less than 2 units, so it rounds the same if lower is away
    bool away {};
    const auto lower {compute_float<T, true>(q + k, w, s, &away)};
    if ((s != 0 && away) || sig == umul128(w, power) << static_cast<unsigned>(s))
    {
        return lower;
    }

    const auto upper {w == UINT64_MAX ? compute_float_bounds<T>(q + k, UINT64_C(1) << 63U, s + 1) : compute_float_bounds<T>(q + k, w + 1U, s)};
    if (lower == upper)
    {
        return upper;
    }

    // The halfway point (2m + 1) * 2^(e - 1) between lower and the next value
    const auto m {lower.power2 == 0 ? lower.mantissa : lower.mantissa | (UINT64_C(1) << format::mantissa_explicit_bits)};
    const std::int64_t e {(lower.power2 == 0 ? 1 : lower.power2) + format::minimum_exponent - format::mantissa_explicit_bits};

    const auto order {compare_exact(q, sig, 2U * m + 1U, e - 1)};
    return order < 0 || (order == 0 && (m & 1U) == 0U) ? lower : upper;
}

// 10^-22 to 10^22, each rounded to nearest
BOOST_DECIMAL_INLINE_CONSTEXPR_VARIABLE double double_powers_of_ten[] = {
    1e-22, 1e-21, 1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
    1e-10, 1e-9,  1e-8,  1e-7,  1e-6,  1e-5,  1e-4,  1e-3,  1e-2,  1e-1,  1e0,   1e1,
    1e2,   1e3,   1e4,   1e5,   1e6,   1e7,   1e8,   1e9,   1e10,  1e11,  1e12,  1e13,
    1e14,  1e15,  1e16,  1e17,  1e18,  1e19,  1e20,  1e21,  1e22
};

// Clinger, if FLT_EVAL_METHOD is 0: with w and 10^|q| exact in double, one multiply or divide rounds correctly. A float
// multiplies by the rounded 10^q, off by under 2 units: it rounds as the exact value unless a midpoint is that near.
template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto clinger_fast_path(const std::int64_t q, const std::uint64_t w, const bool negative, T& value) noexcept -> bool
{
    #if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
    if (q >= -22 && q <= 22 && w <= (UINT64_C(1) << 53U))
    {
        auto result {static_cast<double>(w)};
        BOOST_DECIMAL_IF_CONSTEXPR (std::is_same<T, float>::value)
        {
            result *= double_powers_of_ten[q + 22];
            if ((bit_cast<std::uint64_t>(result) & ((UINT64_C(1) << 29U) - 1U)) - ((UINT64_C(1) << 28U) - 2U) <= 4U)
            {
                return false;
            }
        }
        else
        {
            result = q < 0 ? result / double_powers_of_ten[22 - q] : result * double_powers_of_ten[q + 22];
        }
        value = static_cast<T>(negative ? -result : result);
        return true;
    }
    #else
    static_cast<void>(q);
    static_cast<void>(w);
    static_cast<void>(negative);
    static_cast<void>(value);
    #endif

    return false;
}

template <typename T>
BOOST_DECIMAL_CXX20_CONSTEXPR auto to_binary(const adjusted_mantissa am, const bool negative) noexcept -> T
{
    using format = binary_format<T>;
    using bits_type = typename format::bits_type;

    auto bits {static_cast<bits_type>(am.mantissa)};
    bits |= static_cast<bits_type>(static_cast<bits_type>(am.power2) << format::mantissa_explicit_bits);
    bits |= static_cast<bits_type>(static_cast<bits_type>(negative) << (sizeof(bits_type) * 8U - 1U));

    return bit_cast<T>(bits);
}

} // namespace fast_float
} // namespace detail
} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_FAST_FLOAT_COMPUTE_FLOAT64_HPP
