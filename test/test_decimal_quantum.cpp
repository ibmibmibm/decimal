// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "mini_to_chars.hpp"
#include "testing_config.hpp"
#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <bitset>
#include <limits>
#include <random>
#include <climits>
#include <cmath>
#include <cerrno>
#include <iostream>
#include <iomanip>

using namespace boost::decimal;

#if !defined(BOOST_DECIMAL_REDUCE_TEST_DEPTH)
static constexpr auto N = static_cast<std::size_t>(1024); // Number of trials
#else
static constexpr auto N = static_cast<std::size_t>(1024 >> 4U); // Number of trials
#endif

// NOLINTNEXTLINE : Seed with a constant for repeatability
static std::mt19937_64 rng(42); // NOSONAR : Global rng is not const

template <typename Dec>
void test_nonfinite_samequantum()
{
    constexpr Dec one(1);
    constexpr Dec two(2);
    BOOST_TEST(samequantum(std::numeric_limits<Dec>::infinity(), std::numeric_limits<Dec>::infinity()));
    BOOST_TEST(samequantum(std::numeric_limits<Dec>::quiet_NaN(), std::numeric_limits<Dec>::quiet_NaN()));
    BOOST_TEST(samequantum(std::numeric_limits<Dec>::signaling_NaN(), std::numeric_limits<Dec>::signaling_NaN()));
    BOOST_TEST(!samequantum(std::numeric_limits<Dec>::infinity(), std::numeric_limits<Dec>::quiet_NaN()));
    BOOST_TEST(!samequantum(one, std::numeric_limits<Dec>::infinity()));
    BOOST_TEST(!samequantum(one, std::numeric_limits<Dec>::quiet_NaN()));
    BOOST_TEST(!samequantum(one, std::numeric_limits<Dec>::signaling_NaN()));
    BOOST_TEST(!samequantum(std::numeric_limits<Dec>::infinity(), one));
    BOOST_TEST(!samequantum(std::numeric_limits<Dec>::quiet_NaN(), one));
    BOOST_TEST(!samequantum(std::numeric_limits<Dec>::signaling_NaN(), one));
    BOOST_TEST(samequantum(one, two));
}

template <typename Dec>
void test_same_quantum()
{
    std::uniform_int_distribution<std::int64_t> sig(1'000'000, 9'999'999);
    std::uniform_int_distribution<std::int32_t> exp(std::numeric_limits<Dec>::min_exponent10 + 19,
                                                    std::numeric_limits<Dec>::max_exponent10 - 19);

    constexpr auto max_iter {std::is_same<Dec, decimal128_t>::value ? N / 4 : N};
    for (std::size_t i {}; i < max_iter; ++i)
    {
        auto exp1 {exp(rng)};
        auto exp2 {exp(rng)};

        const Dec val1 {sig(rng), exp1};
        const Dec val2 {sig(rng), exp2};


        if (exp1 == exp2)
        {
            if (!BOOST_TEST(samequantum(val1, val2)))
            {
                // LCOV_EXCL_START
                std::cerr << "Val 1: " << val1
                          << "\nVal 2: " << val2 << std::endl;
                // LCOV_EXCL_STOP
            }
        }
        else
        {
            if (!BOOST_TEST(!samequantum(val1, val2)))
            {
                // LCOV_EXCL_START
                std::cerr << "Val 1: " << val1
                          << "\nVal 2: " << val2 << std::endl;
                // LCOV_EXCL_STOP
            }
        }
    }
}

// quantexp returns the quantum exponent, not the biased encoded exponent field,
// so the IEEE types report back exactly the exponent they were constructed with.
template <typename Dec>
void test_quantexp()
{
    // Loop through every possible exponent
    for (auto i {std::numeric_limits<Dec>::min_exponent10}; i < std::numeric_limits<Dec>::max_exponent10 - std::numeric_limits<Dec>::digits10; ++i)
    {
        const Dec val1 {1, i};

        if (isinf(val1))
        {
            continue; // LCOV_EXCL_LINE
        }
        else
        {
            if (!BOOST_TEST_EQ(quantexp(val1), i))
            // LCOV_EXCL_START
            {
                std::cerr << "Val: " << val1
                          << "\nExp 1: " << i << std::endl;
            }
            // LCOV_EXCL_STOP
        }
    }
}

// Spot values from the bug report that the biased return value got wrong.
// The IEEE types preserve the quantum they were constructed with.
template <typename Dec>
void test_quantexp_spot()
{
    BOOST_TEST_EQ(quantexp(Dec {3, 2}), 2);
    BOOST_TEST_EQ(quantexp(Dec {3000000, -4}), -4);
    BOOST_TEST_EQ(quantexp(Dec {5, 0}), 0);

    // The quantum exponent belongs to the cohort member, not to the value, so
    // two equal values with different quanta report different exponents.
    BOOST_TEST((Dec {3, 2} == Dec {3000000, -4}));
    BOOST_TEST_NE(quantexp(Dec {3, 2}), quantexp(Dec {3000000, -4}));
}

// samequantum compares two exponents, so it is immune to a uniform bias offset.
// Holding it against quantexp pins the two functions to the same definition.
template <typename Dec>
void test_quantexp_samequantum_agreement()
{
    std::uniform_int_distribution<std::int64_t> sig(1, 9'999'999);
    std::uniform_int_distribution<std::int32_t> exp(std::numeric_limits<Dec>::min_exponent10 + 19,
                                                    std::numeric_limits<Dec>::max_exponent10 - 19);

    constexpr auto max_iter {std::is_same<Dec, decimal128_t>::value ? N / 4 : N};
    for (std::size_t i {}; i < max_iter; ++i)
    {
        const Dec val1 {sig(rng), exp(rng)};
        const Dec val2 {sig(rng), exp(rng)};

        if (!isfinite(val1) || !isfinite(val2))
        {
            continue; // LCOV_EXCL_LINE
        }

        if (!BOOST_TEST_EQ(samequantum(val1, val2), quantexp(val1) == quantexp(val2)))
        // LCOV_EXCL_START
        {
            std::cerr << "Val 1: " << val1
                      << "\nVal 2: " << val2 << std::endl;
        }
        // LCOV_EXCL_STOP
    }
}

template <typename Dec>
void test_nonfinite_quantexp()
{
    BOOST_TEST_EQ(quantexp(std::numeric_limits<Dec>::infinity()), INT_MIN);
    BOOST_TEST_EQ(quantexp(-std::numeric_limits<Dec>::infinity()), INT_MIN);
    BOOST_TEST_EQ(quantexp(std::numeric_limits<Dec>::quiet_NaN()), INT_MIN);
    BOOST_TEST_EQ(quantexp(-std::numeric_limits<Dec>::quiet_NaN()), INT_MIN);
    BOOST_TEST_EQ(quantexp(std::numeric_limits<Dec>::signaling_NaN()), INT_MIN);
    BOOST_TEST_EQ(quantexp(-std::numeric_limits<Dec>::signaling_NaN()), INT_MIN);
}

template <typename Dec>
void test_quantize()
{
    using sig_type = typename Dec::significand_type;

    std::uniform_int_distribution<std::uint64_t> sig(1'000'000, 9'999'999);
    std::uniform_int_distribution<std::int32_t> exp(std::numeric_limits<Dec>::min_exponent10 + std::numeric_limits<Dec>::digits10 + 1,
                                                    std::numeric_limits<Dec>::max_exponent10 - std::numeric_limits<Dec>::digits10 - 1);

    constexpr auto max_iter {std::is_same<Dec, decimal128_t>::value ? N / 4 : N};
    for (std::size_t i {}; i < max_iter; ++i)
    {
        const auto sig1 {static_cast<sig_type>(sig(rng))};
        const auto sig2 {static_cast<sig_type>(sig(rng))};
        const auto exp1 {exp(rng)};

        const Dec val1 {sig1, exp1};
        // Both values share the quantum exponent so the IEEE 754-2008 quantize
        // operation reduces to a no-op that returns lhs unchanged. Both
        // significands have the same digit count, so the fast types normalize
        // them identically and end up with the same stored exponent.
        const Dec val2 {sig2, exp1};

        if (!BOOST_TEST_EQ(quantize(val1, val2), val1))
        {
            // LCOV_EXCL_START
            std::cerr << std::setprecision(std::numeric_limits<Dec>::digits10)
                      << "Val 1: " << val1
                      << "\nVal 2: " << val2
                      << "\n Func: " << quantize(val1, val2) << std::endl;
            // LCOV_EXCL_STOP
        }
    }
}

template <typename Dec>
void test_quantize_spot()
{
    // IEEE 754-2008 5.3.2: quantize(x, y) returns a value equal to x (modulo
    // rounding) but with the quantum exponent of y.

    const auto check = [](Dec x, Dec y, Dec expected)
    {
        if (!BOOST_TEST_EQ(quantize(x, y), expected))
        {
            // LCOV_EXCL_START
            std::cerr << std::setprecision(std::numeric_limits<Dec>::digits10)
                      << "x:        " << x
                      << "\ny:        " << y
                      << "\nexpected: " << expected
                      << "\nactual:   " << quantize(x, y) << std::endl;
            // LCOV_EXCL_STOP
        }
    };

    // Same exponent: identity
    check(Dec{1, -1}, Dec{1, -1}, Dec{1, -1});

    // Smaller rhs exponent: significand grows, value is exact
    check(Dec{217}, Dec{1, -1}, Dec{2170, -1});

    // Larger rhs exponent: significand shrinks with rounding
    check(Dec{217}, Dec{1, 1}, Dec{22, 1});
    check(Dec{217}, Dec{1, 2}, Dec{2, 2});

    // Rounding all digits away rounds to zero or one based on the leading digit
    check(Dec{1, -1}, Dec{1, 1}, Dec{0, 1});
    check(Dec{1, -1}, Dec{1, 2}, Dec{0, 2});

    // Sign of x is preserved on non-zero results
    check(Dec{-217}, Dec{1, -1}, Dec{-2170, -1});
}

template <typename Dec>
void test_nonfinite_quantize()
{
    const Dec one{1};
    std::uniform_int_distribution<std::int64_t> dist(1'000'000, 9'999'999);

    BOOST_TEST(isnan(quantize(std::numeric_limits<Dec>::quiet_NaN(), one)));
    BOOST_TEST(isnan(quantize(one, std::numeric_limits<Dec>::quiet_NaN())));
    BOOST_TEST(isnan(quantize(std::numeric_limits<Dec>::signaling_NaN(), one)));
    BOOST_TEST(isnan(quantize(one, std::numeric_limits<Dec>::signaling_NaN())));
    BOOST_TEST(isnan(quantize(std::numeric_limits<Dec>::infinity(), one)));
    BOOST_TEST(isnan(quantize(one, std::numeric_limits<Dec>::infinity())));
    BOOST_TEST(isinf(quantize(std::numeric_limits<Dec>::infinity() * dist(rng),
                              std::numeric_limits<Dec>::infinity() * dist(rng))));
}

// IEEE 754 3.5.2 gives a coefficient above the precision the value zero, as from_bid does. from_bits
// keeps such a coefficient, and quantize gives the result of the constructor: this is library behavior
template <typename T, typename Bits, typename Sig>
void test_quantize_non_canonical(const Bits bits, const Sig sig)
{
    for (const bool negative : {false, true})
    {
        const auto x {negative ? -from_bits(bits) : from_bits(bits)};
        BOOST_TEST(to_bid(quantize(x, T{1, 0})) == to_bid(T{sig, 0, negative}));
        BOOST_TEST(to_bid(quantize(x, T{1, 1})) == to_bid(T{static_cast<Sig>((sig + 5U) / 10U), 1, negative}));

        const auto standard {negative ? -from_bid<T>(bits) : from_bid<T>(bits)};
        BOOST_TEST(to_bid(quantize(standard, T{1, 1})) == to_bid(T{0U, 1, negative}));
    }
}

// A coefficient above 64 bits, q * 100 + tail, quantized two digits up. The remainder tells an exact tie
// from the values next to it, in each rounding mode and for both signs. up says that the result is q + 1
void test_quantize_wide_rounding()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    using boost::int128::uint128_t;
    struct test_case { uint128_t q; unsigned tail; rounding_mode mode; bool up_positive; bool up_negative; };
    const uint128_t q {(uint128_t{1U} << 64U) + 2U};
    const test_case cases[] {
        {q, 49U, rounding_mode::fe_dec_to_nearest, false, false},
        {q, 50U, rounding_mode::fe_dec_to_nearest, false, false},
        {q + 1U, 50U, rounding_mode::fe_dec_to_nearest, true, true},
        {q, 51U, rounding_mode::fe_dec_to_nearest, true, true},
        {q, 49U, rounding_mode::fe_dec_to_nearest_from_zero, false, false},
        {q, 50U, rounding_mode::fe_dec_to_nearest_from_zero, true, true},
        {q, 1U, rounding_mode::fe_dec_toward_zero, false, false},
        {q, 99U, rounding_mode::fe_dec_toward_zero, false, false},
        {q, 1U, rounding_mode::fe_dec_upward, true, false},
        {q, 0U, rounding_mode::fe_dec_upward, false, false},
        {q, 1U, rounding_mode::fe_dec_downward, false, true},
        {q, 0U, rounding_mode::fe_dec_downward, false, false},
    };

    for (const auto& c : cases)
    {
        fesetround(c.mode);
        for (const bool negative : {false, true})
        {
            const decimal128_t x {c.q * 100U + c.tail, 0, negative};
            const auto up {negative ? c.up_negative : c.up_positive};
            const auto result {quantize(x, decimal128_t{1, 2})};
            BOOST_TEST(to_bid(result) == to_bid(decimal128_t{up ? c.q + 1U : c.q, 2, negative}));
            BOOST_TEST_EQ(quantexp(result), 2);
        }
    }
    fesetround(rounding_mode::fe_dec_default);
    #endif
}

// The fast types are absent from every test below on purpose. They normalize on
// construction so they have no cohorts, and samequantum, quantexp, and quantize are all
// rejected for them at compile time. See test_illegal_fast_samequantum.cpp,
// test_illegal_fast_quantexp.cpp, and test_illegal_fast_quantize.cpp, and use rescale to
// round a fast type to a fixed number of fractional digits.
int main()
{
    test_same_quantum<decimal32_t>();
    test_nonfinite_samequantum<decimal32_t>();
    test_quantexp<decimal32_t>();
    test_quantexp_spot<decimal32_t>();
    test_quantexp_samequantum_agreement<decimal32_t>();
    test_nonfinite_quantexp<decimal32_t>();
    test_quantize<decimal32_t>();
    test_quantize_spot<decimal32_t>();
    test_nonfinite_quantize<decimal32_t>();

    test_same_quantum<decimal64_t>();
    test_nonfinite_samequantum<decimal64_t>();
    test_quantexp<decimal64_t>();
    test_quantexp_spot<decimal64_t>();
    test_quantexp_samequantum_agreement<decimal64_t>();
    test_nonfinite_quantexp<decimal64_t>();
    test_quantize<decimal64_t>();
    test_quantize_spot<decimal64_t>();
    test_nonfinite_quantize<decimal64_t>();

    test_same_quantum<decimal128_t>();
    test_nonfinite_samequantum<decimal128_t>();
    test_quantexp<decimal128_t>();
    test_quantexp_spot<decimal128_t>();
    test_quantexp_samequantum_agreement<decimal128_t>();
    test_nonfinite_quantexp<decimal128_t>();
    test_quantize<decimal128_t>();
    test_quantize_spot<decimal128_t>();
    test_nonfinite_quantize<decimal128_t>();

    test_quantize_wide_rounding();

    using boost::int128::uint128_t;
    test_quantize_non_canonical<decimal32_t>(std::uint32_t{UINT32_C(0x60000000) | (UINT32_C(101) << 21) | UINT32_C(0x1FFFFF)}, UINT32_C(10485759));
    test_quantize_non_canonical<decimal64_t>(std::uint64_t{(UINT64_C(3) << 61) | (UINT64_C(398) << 51) | ((UINT64_C(1) << 51) - 1U)},
                                             (UINT64_C(1) << 53) + (UINT64_C(1) << 51) - 1U);
    test_quantize_non_canonical<decimal128_t>((uint128_t{6176U} << 113) | ((uint128_t{1U} << 113) - 1U), (uint128_t{1U} << 113) - 1U);

    return boost::report_errors();
}
