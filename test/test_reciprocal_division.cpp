// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The divisions by a reciprocal must give the quotient and the remainder of integer division

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>
#include <limits>
#include <random>

using namespace boost::decimal;
using boost::int128::uint128_t;
using boost::decimal::detail::u256;

static std::mt19937_64 rng(42);

constexpr rounding_mode modes[] {rounding_mode::fe_dec_to_nearest,
                                 rounding_mode::fe_dec_to_nearest_from_zero,
                                 rounding_mode::fe_dec_upward,
                                 rounding_mode::fe_dec_downward,
                                 rounding_mode::fe_dec_toward_zero};

auto ten(const std::uint64_t, const int n) -> std::uint64_t { return detail::pow10(static_cast<std::uint64_t>(n)); }
auto ten(const uint128_t&, const int n) -> uint128_t { return detail::pow10(uint128_t{static_cast<std::uint64_t>(n)}); }
auto ten(const u256&, const int n) -> u256 { return detail::pow10(u256{static_cast<std::uint64_t>(n)}); }

auto random_value(const std::uint64_t) -> std::uint64_t { return rng(); }
auto random_value(const uint128_t&) -> uint128_t { return uint128_t{rng(), rng()}; }
auto random_value(const u256&) -> u256 { return u256{rng(), rng(), rng(), rng()}; }

// A value in [0, bound)
template <typename Wide>
auto random_below(const Wide& bound) -> Wide
{
    const Wide value {random_value(Wide{})};
    return value - (value / bound) * bound;
}

// divmod_pow10_dispatch of one word against the operators / and %
template <typename T>
void check_narrow(const T coeff, const int shift)
{
    const T power {detail::pow10(static_cast<T>(shift))};
    const auto result {detail::impl::divmod_pow10_dispatch(coeff, shift, power)};
    BOOST_TEST_EQ(result.quotient, static_cast<T>(coeff / power));
    BOOST_TEST_EQ(result.remainder, static_cast<T>(coeff % power));
}

template <typename T>
void test_narrow()
{
    constexpr T max_value {(std::numeric_limits<T>::max)()};

    for (int shift {0}; shift <= std::numeric_limits<T>::digits10; ++shift)
    {
        const T power {detail::pow10(static_cast<T>(shift))};

        for (const T coeff : {T{0}, T{1}, static_cast<T>(max_value - 1U), max_value})
        {
            check_narrow(coeff, shift);
        }

        // An exact multiple and the values next to it
        const T max_multiplier {static_cast<T>(max_value / power)};
        for (const T multiplier : {T{1}, T{2}, T{5}, T{9}, T{10}, T{99}, static_cast<T>(max_multiplier / 2U), max_multiplier})
        {
            if (multiplier <= max_multiplier)
            {
                const T multiple {static_cast<T>(multiplier * power)};
                check_narrow(static_cast<T>(multiple - 1U), shift);
                check_narrow(multiple, shift);
                check_narrow(static_cast<T>(multiple + 1U), shift);
            }
        }

        // The values where the number of digits or of bits changes
        for (int digits {1}; digits <= std::numeric_limits<T>::digits10; ++digits)
        {
            const T threshold {detail::pow10(static_cast<T>(digits))};
            check_narrow(static_cast<T>(threshold - 1U), shift);
            check_narrow(threshold, shift);
        }
        for (int bits {1}; bits < std::numeric_limits<T>::digits; ++bits)
        {
            const T threshold {static_cast<T>(T{1} << bits)};
            check_narrow(static_cast<T>(threshold - 1U), shift);
            check_narrow(threshold, shift);
        }

        for (int i {0}; i < 2000; ++i)
        {
            const T coeff {static_cast<T>(rng())};
            check_narrow(coeff, shift);
            check_narrow(static_cast<T>(coeff >> (rng() % std::numeric_limits<T>::digits)), shift);
        }
    }
}

void check_uint128(const uint128_t& coeff, const int shift)
{
    const uint128_t power {ten(uint128_t{}, shift)};
    const auto result {detail::impl::divmod_pow10_dispatch(coeff, shift, power)};
    BOOST_TEST(result.quotient == coeff / power);
    BOOST_TEST(result.remainder == coeff % power);
}

void test_uint128()
{
    constexpr uint128_t max_value {UINT64_MAX, UINT64_MAX};

    for (int shift {0}; shift <= 38; ++shift)
    {
        const uint128_t power {ten(uint128_t{}, shift)};

        for (const uint128_t& coeff : {uint128_t{0U}, uint128_t{1U}, max_value - 1U, max_value})
        {
            check_uint128(coeff, shift);
        }

        const uint128_t max_multiplier {max_value / power};
        for (const uint128_t& multiplier : {uint128_t{1U}, uint128_t{2U}, uint128_t{5U}, uint128_t{9U}, uint128_t{10U}, max_multiplier / 2U, max_multiplier})
        {
            if (multiplier <= max_multiplier)
            {
                const uint128_t multiple {multiplier * power};
                check_uint128(multiple - 1U, shift);
                check_uint128(multiple, shift);
                check_uint128(multiple + 1U, shift);
            }
        }

        for (int digits {1}; digits <= 38; ++digits)
        {
            const uint128_t threshold {ten(uint128_t{}, digits)};
            check_uint128(threshold - 1U, shift);
            check_uint128(threshold, shift);
        }
        for (int bits {1}; bits < 128; ++bits)
        {
            const uint128_t threshold {uint128_t{1U} << bits};
            check_uint128(threshold - 1U, shift);
            check_uint128(threshold, shift);
        }

        for (int i {0}; i < 2000; ++i)
        {
            const uint128_t coeff {random_value(uint128_t{})};
            check_uint128(coeff, shift);
            check_uint128(coeff >> static_cast<int>(rng() % 128U), shift);
        }
    }
}

// The reference for the step of a quotient to the next value
auto steps_up(const rounding_mode mode, const bool negative, const bool odd, const bool inexact, const bool half, const bool above_half) -> bool
{
    switch (mode)
    {
        case rounding_mode::fe_dec_to_nearest:
            return above_half || (half && odd);
        case rounding_mode::fe_dec_to_nearest_from_zero:
            return above_half || half;
        case rounding_mode::fe_dec_upward:
            return inexact && !negative;
        case rounding_mode::fe_dec_downward:
            return inexact && negative;
        default:
            return false;
    }
}

// A product of 2p - 1 or 2p digits rounds to p digits. The reference divides with the operator /.
template <typename T, typename Sig, typename Wide, typename Finalize>
void check_product(const Wide& product, Finalize finalize)
{
    constexpr int p {std::numeric_limits<T>::digits10};
    int extra {product >= ten(Wide{}, 2 * p - 1) ? p : p - 1};
    const Wide power {ten(Wide{}, extra)};
    const Wide half {power >> 1};
    const Wide quotient {product / power};
    const Wide remainder {product - quotient * power};
    const bool odd {(static_cast<std::uint64_t>(quotient) & 1U) != 0U};

    for (const auto mode : modes)
    {
        for (const bool negative : {false, true})
        {
            Wide rounded {quotient};
            int exp {extra};
            if (steps_up(mode, negative, odd, remainder != Wide{}, remainder == half, remainder > half))
            {
                rounded = rounded + ten(Wide{}, 0);
            }
            // A quotient of all nines steps up to the next exponent
            if (rounded == ten(Wide{}, p))
            {
                rounded = ten(Wide{}, p - 1);
                ++exp;
            }

            const T expected {static_cast<Sig>(rounded), exp, negative};
            const T result {finalize(product, negative, mode)};
            BOOST_TEST(result == expected);
        }
    }
}

template <typename T, typename Sig, typename Wide, typename Finalize>
void test_products(Finalize finalize)
{
    constexpr int p {std::numeric_limits<T>::digits10};
    const Wide one {ten(Wide{}, 0)};
    const Wide min_sig {ten(Wide{}, p - 1)};
    const Wide sig_bound {ten(Wide{}, p)};
    const Wide max_product {(sig_bound - one) * (sig_bound - one)};

    for (const int extra : {p - 1, p})
    {
        const Wide power {ten(Wide{}, extra)};
        const Wide half {power >> 1};

        // An exact multiple needs the correction if the first quotient is one too low. The
        // quotients are the limits of the p digits, and the remainders are next to 0, 1/2 and 1.
        const Wide quotients[] {min_sig, min_sig + one, min_sig + one + one, sig_bound >> 1, sig_bound - one - one - one,
                                sig_bound - one - one, sig_bound - one};
        const Wide remainders[] {Wide{}, one, one + one, half - one, half, half + one, power - one - one, power - one};

        for (const auto& quotient : quotients)
        {
            for (const auto& remainder : remainders)
            {
                const Wide product {quotient * power + remainder};
                if (product <= max_product)
                {
                    check_product<T, Sig, Wide>(product, finalize);
                }
            }
        }

        for (int i {0}; i < 2000; ++i)
        {
            const Wide quotient {min_sig + random_below(sig_bound - min_sig)};
            const Wide product {quotient * power + remainders[rng() % 8U]};
            if (product <= max_product)
            {
                check_product<T, Sig, Wide>(product, finalize);
                check_product<T, Sig, Wide>(quotient * power + random_below(power), finalize);
            }
        }
    }

    // The values next to the threshold between 2p - 1 and 2p digits
    const Wide threshold {ten(Wide{}, 2 * p - 1)};
    check_product<T, Sig, Wide>(threshold, finalize);
    Wide distance {one};
    for (int bits {0}; bits < 3 * p; ++bits)
    {
        check_product<T, Sig, Wide>(threshold - distance, finalize);
        check_product<T, Sig, Wide>(threshold + distance, finalize);
        distance = distance + distance;
    }

    // The products of two significands of p digits
    for (int i {0}; i < 4000; ++i)
    {
        const Wide lhs {min_sig + random_below(sig_bound - min_sig)};
        const Wide rhs {min_sig + random_below(sig_bound - min_sig)};
        check_product<T, Sig, Wide>(lhs * rhs, finalize);
    }
    check_product<T, Sig, Wide>(min_sig * min_sig, finalize);
    check_product<T, Sig, Wide>(max_product, finalize);
}

template <typename T>
void test_products_32()
{
    test_products<T, std::uint32_t, std::uint64_t>([](const std::uint64_t product, const bool negative, const rounding_mode mode)
        { return detail::impl::mul_finalize_u64<T>(product, 0, negative, mode); });
}

template <typename T>
void test_products_64()
{
    test_products<T, std::uint64_t, uint128_t>([](const uint128_t& product, const bool negative, const rounding_mode mode)
        { return detail::impl::mul_finalize_u128<T>(product, 0, negative, mode); });
}

template <typename T>
void test_products_128()
{
    test_products<T, uint128_t, u256>([](const u256& product, const bool negative, const rounding_mode mode)
        { return detail::impl::mul_finalize_u256<T>(product, 0, negative, mode); });
}

// from_bits keeps a coefficient above the largest one, which IEEE 754 3.5.2 reads as zero.
// This checks only that the product finishes, not its value.
void test_coefficient_above_maximum()
{
    const decimal32_t a {from_bits(UINT32_C(0x6CBFFFFF))};
    const decimal64_t b {from_bits(UINT64_C(0x6C77FFFFFFFFFFFF))};
    const decimal128_t c {from_bits(uint128_t{UINT64_C(0x6C107FFFFFFFFFFF), UINT64_MAX})};

    for (const auto mode : modes)
    {
        fesetround(mode);
        BOOST_TEST(isfinite(a * a));
        BOOST_TEST(isfinite(b * b));
        BOOST_TEST(isfinite(c * c));
    }
    fesetround(rounding_mode::fe_dec_to_nearest);
}

int main()
{
    test_narrow<std::uint32_t>();
    test_narrow<std::uint64_t>();
    test_uint128();

    test_products_32<decimal32_t>();
    test_products_32<decimal_fast32_t>();
    test_products_64<decimal64_t>();
    test_products_64<decimal_fast64_t>();
    test_products_128<decimal128_t>();
    test_products_128<decimal_fast128_t>();

    test_coefficient_above_maximum();

    return boost::report_errors();
}
