// Copyright 2023 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/decimal.hpp>
#include <boost/decimal/detail/fast_float/compute_float80_128.hpp>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <climits>
#include <cstring>

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#include <boost/core/lightweight_test.hpp>

void test_compute_float32()
{
    using namespace boost::decimal::detail::fast_float;
    const auto compute_float32 = [](const std::int64_t q, const std::uint64_t w, const bool negative)
    { return to_binary<float>(compute_float<float>(q, w), negative); };

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<std::int64_t> dist (1, 1);

    // Trivial verifcation
    BOOST_TEST_EQ(compute_float32(1 * dist(gen), 1, false), 1e1F);
    BOOST_TEST_EQ(compute_float32(0 * dist(gen), 1, true), -1e0F);
    BOOST_TEST_EQ(compute_float32(38 * dist(gen), 1, false), 1e38F);

    // out of range
    BOOST_TEST_EQ(compute_float32(310 * dist(gen), 5, false), HUGE_VALF);
    BOOST_TEST_EQ(compute_float32(310 * dist(gen), 5, true), -HUGE_VALF);
    BOOST_TEST_EQ(compute_float32(1000 * dist(gen), 5, false), HUGE_VALF);
    BOOST_TEST_EQ(compute_float32(1000 * dist(gen), 5, true), -HUGE_VALF);
    BOOST_TEST_EQ(compute_float32(-325 * dist(gen), 5, false), 0.0F);

    // Composite
    BOOST_TEST_EQ(compute_float32(10 * dist(gen), 123456789, false), 123456789e10F);
    BOOST_TEST_EQ(compute_float32(20 * dist(gen), UINT64_C(444444444), false), 444444444e20F);
}

void test_compute_float64()
{
    using namespace boost::decimal::detail::fast_float;
    const auto compute_float64 = [](const std::int64_t q, const std::uint64_t w, const bool negative)
    { return to_binary<double>(compute_float<double>(q, w), negative); };

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<std::uint64_t> dist (1, 1);

    // Trivial verification
    BOOST_TEST_EQ(compute_float64(1, dist(gen), false), 1e1);
    BOOST_TEST_EQ(compute_float64(0, dist(gen), true), -1e0);
    BOOST_TEST_EQ(compute_float64(308, dist(gen), false), 1e308);

    // out of range
    BOOST_TEST_EQ(compute_float64(310, 5U * dist(gen), false), HUGE_VAL);
    BOOST_TEST_EQ(compute_float64(310, 5U * dist(gen), true), -HUGE_VAL);
    BOOST_TEST_EQ(compute_float64(1000, 5U * dist(gen), false), HUGE_VAL);
    BOOST_TEST_EQ(compute_float64(1000, 5U * dist(gen), true), -HUGE_VAL);
    BOOST_TEST_EQ(compute_float64(-325, 5U * dist(gen), false), 0.0);
    BOOST_TEST_EQ(compute_float64(static_cast<std::int64_t>(dist(gen)) * 50U, 0, false), 0.0);
    BOOST_TEST_EQ(compute_float64(static_cast<std::int64_t>(dist(gen)) * 50U, 0, true), 0.0);
    BOOST_TEST_EQ(compute_float64(300, UINT64_MAX, false), HUGE_VAL * static_cast<double>(dist(gen)));

    // Composite
    BOOST_TEST_EQ(compute_float64(10 * static_cast<std::int64_t>(dist(gen)), 123456789, false), 123456789e10);
    BOOST_TEST_EQ(compute_float64(100 * static_cast<std::int64_t>(dist(gen)), UINT64_C(4444444444444444444), false), 4444444444444444444e100);
    BOOST_TEST_EQ(compute_float64(100 * static_cast<std::int64_t>(dist(gen)), std::numeric_limits<std::uint64_t>::max(), false), 18446744073709551615e100);
    BOOST_TEST_EQ(compute_float64(100 * static_cast<std::int64_t>(dist(gen)), UINT64_C(10000000000000000000), false), 10000000000000000000e100);
}

void test_compute_float80_128()
{
    using boost::decimal::detail::fast_float::compute_float80_128;

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<std::uint64_t> dist (1, 1);

    BOOST_TEST_EQ(compute_float80_128(1, dist(gen), false), 1e1L);
    BOOST_TEST_EQ(compute_float80_128(0, dist(gen), true), -1e0L);
    BOOST_TEST_EQ(compute_float80_128(-1, dist(gen), false), 1e-1L);
    BOOST_TEST_EQ(compute_float80_128(100, std::uint64_t {0}, false), 0.0L);
}

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

// Ryu converts only the types wider than double
template <typename T>
void test_ryu()
{
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<T> dist (1, 2);

    using namespace boost::decimal::detail::ryu;

    BOOST_TEST(floating_point_to_fd128(T(0) * dist(gen)).mantissa == 0U);
    BOOST_TEST(floating_point_to_fd128(std::numeric_limits<T>::infinity() * dist(gen)).exponent == fd128_exceptional_exponent);
    BOOST_TEST(floating_point_to_fd128(std::numeric_limits<T>::quiet_NaN() * dist(gen)).exponent == fd128_exceptional_exponent);
    BOOST_TEST(floating_point_to_fd128(-std::numeric_limits<T>::infinity() * dist(gen)).exponent == fd128_exceptional_exponent);
    BOOST_TEST(floating_point_to_fd128(-std::numeric_limits<T>::quiet_NaN() * dist(gen)).exponent == fd128_exceptional_exponent);

    // The lower end of the interval has trailing zeros in its last digits, so the shortest digits are 646
    #if BOOST_DECIMAL_LDBL_BITS == 80
    BOOST_TEST(floating_point_to_fd128(6.46e26L).mantissa == 646U);
    BOOST_TEST(floating_point_to_fd128(1.0435e26L).mantissa == 10435U);
    BOOST_TEST(floating_point_to_fd128(656177904623787416832.0L).mantissa == UINT64_C(6561779046237874168));
    BOOST_TEST(floating_point_to_fd128(static_cast<long double>(UINT64_C(12500000000000004688)) * 32768).mantissa ==
               UINT64_C(4096000000000001536));
    #endif

    // The removed digits are 560 and not a tie, so the last digit rounds up
    #if BOOST_DECIMAL_LDBL_BITS == 80
    const auto fd {floating_point_to_fd128(74787340373893776834560.0L)};
    BOOST_TEST(fd.mantissa == (unsigned_128_type {UINT64_C(4), UINT64_C(1000364079055570371)}));
    #else
    const auto fd {floating_point_to_fd128(34623143365181110578251157253715394560.0L)};
    BOOST_TEST(fd.mantissa == (unsigned_128_type {UINT64_C(1876924362740321), UINT64_C(17418691559299806659)}));
    #endif
    BOOST_TEST_EQ(fd.exponent, 3);

    // A binary128 integer below 2^120 keeps every digit, so the digits are exact
    const auto exact {generic_binary_to_decimal(unsigned_128_type {UINT64_C(0x4072000000000000), UINT64_C(1)}, 112, 15, false)};
    BOOST_TEST(exact.mantissa == (unsigned_128_type {UINT64_C(2251799813685248), UINT64_C(8)}));
    BOOST_TEST_EQ(exact.exponent, 0);
    BOOST_TEST_EQ(exact.direction, 0);
}

// Shortest digits one past the precision that end in 5 round toward the exact value of the binary number
void test_shortest_tie()
{
    using namespace boost::decimal;
    using namespace boost::decimal::literals;

    BOOST_TEST_EQ(decimal32_t {1.0000025F}, 1.000003_DF);
    BOOST_TEST_EQ(decimal32_t {1.0000035F}, 1.000003_DF);
    BOOST_TEST_EQ(decimal_fast32_t {1.0000025F}, 1.000003_DFF);
    BOOST_TEST_EQ(decimal64_t {2.9576661961133885e-8}, 2.957666196113389e-8_DD);
    BOOST_TEST_EQ(decimal64_t {1.1430767782883795e-192}, 1.143076778288379e-192_DD);
    BOOST_TEST_EQ(decimal_fast64_t {2.9576661961133885e-8}, 2.957666196113389e-8_DDF);
    BOOST_TEST_EQ(decimal32_t {5.9139005e-26}, 5.913901e-26_DF);
    BOOST_TEST_EQ(decimal32_t {0.87616956F}, 0.8761696_DF);
    BOOST_TEST_EQ(decimal64_t {-2.9362768111642555e224}, "-2.936276811164256e224"_DD);
    BOOST_TEST_EQ(decimal64_t {2750498739004538.5}, 2750498739004538_DD);
    BOOST_TEST_EQ(decimal_fast32_t {2.4941545e31}, 2.494155e31_DFF);

    // A float on the tie goes to the even digit, and one with two digits past the precision rounds both
    BOOST_TEST_EQ(decimal32_t {574335.25F}, 574335.2_DF);
    BOOST_TEST_EQ(decimal32_t {108435896.0F}, 1.084359e8_DF);

    // A subnormal double has its exact value below the tie, and a normal decimal result
    BOOST_TEST_EQ(decimal64_t {1.9103049429316615e-308}, 1.910304942931661e-308_DD);
    BOOST_TEST_EQ(decimal_fast64_t {1.4615434887007215e-308}, 1.461543488700721e-308_DDF);

    #if BOOST_DECIMAL_LDBL_BITS == 80
    BOOST_TEST_EQ(decimal64_t {268.41196881091365L}, 268.4119688109137_DD);
    BOOST_TEST_EQ(decimal64_t {6.1440403750342675L}, 6.144040375034267_DD);
    #elif BOOST_DECIMAL_LDBL_BITS == 128
    BOOST_TEST_EQ(decimal128_t {0.039266412274453850025464945980491025L}, "0.03926641227445385002546494598049103"_DL);
    BOOST_TEST_EQ(decimal128_t {0.19676412089341422684781472488646195L}, "0.1967641208934142268478147248864619"_DL);
    #endif

    // A wide number on a midpoint is a tie to even, and one next to 10^-42 carries to it without trailing zeros
    #if BOOST_DECIMAL_LDBL_BITS == 80 || BOOST_DECIMAL_LDBL_BITS == 128
    BOOST_TEST_EQ(decimal32_t {12345675.0L}, 1.234568e7_DF);
    BOOST_TEST_EQ(decimal_fast32_t {12345665.0L}, 1.234566e7_DFF);
    BOOST_TEST_EQ(decimal64_t {12345678901234565.0L}, 1.234567890123456e16_DD);
    BOOST_TEST_EQ(quantexp(decimal32_t {1e-42L}), -42);
    BOOST_TEST_EQ(quantexp(decimal64_t {1e-42L}), -42);
    #endif
}

// A value beside a midpoint of a type narrower than float rounds toward it, also when float rounds it onto the midpoint
void test_narrow_midpoint()
{
    using namespace boost::decimal;
    using namespace boost::decimal::literals;

    #ifdef BOOST_DECIMAL_HAS_FLOAT16
    BOOST_TEST(static_cast<std::float16_t>(2049.000000000001_DD) == static_cast<std::float16_t>(2050));
    BOOST_TEST(static_cast<std::float16_t>(2050.999999999999_DD) == static_cast<std::float16_t>(2050));
    BOOST_TEST(static_cast<std::float16_t>("2049.000000000000000000000000000001"_DL) == static_cast<std::float16_t>(2050));
    BOOST_TEST(static_cast<std::float16_t>("2.980232238769531250000000000000001e-8"_DL) > static_cast<std::float16_t>(0));
    BOOST_TEST(static_cast<std::float16_t>("2.98023223876953125e-8"_DL) == static_cast<std::float16_t>(0));
    BOOST_TEST(static_cast<std::float16_t>("9.853363e-3"_DF) == static_cast<std::float16_t>(0x1.42cp-7));
    BOOST_TEST(static_cast<std::float16_t>("9.853363e-3"_DFF) == static_cast<std::float16_t>(0x1.42cp-7));
    BOOST_TEST(static_cast<std::float16_t>("1.00048828126"_DDF) == static_cast<std::float16_t>(0x1.004p+0));
    BOOST_TEST(static_cast<std::float16_t>("1.00048828126"_DLF) == static_cast<std::float16_t>(0x1.004p+0));
    BOOST_TEST(static_cast<std::float16_t>("-1.762509346008301e-4"_DD) == static_cast<std::float16_t>(-0xb.8ep-16));
    BOOST_TEST(static_cast<std::float16_t>("-2.974271774291992e-4"_DD) == static_cast<std::float16_t>(-0x9.bep-15));
    BOOST_TEST(static_cast<std::float16_t>("-1.8921529734809417935e-2"_DL) == static_cast<std::float16_t>(-0x9.bp-9));
    #endif

    #ifdef BOOST_DECIMAL_HAS_BRAINFLOAT16
    BOOST_TEST(static_cast<std::bfloat16_t>(257.0000000001_DD) == static_cast<std::bfloat16_t>(258));
    BOOST_TEST(static_cast<std::bfloat16_t>(258.9999999999_DD) == static_cast<std::bfloat16_t>(258));
    BOOST_TEST(static_cast<std::bfloat16_t>("8.240022e-10"_DF) == static_cast<std::bfloat16_t>(0x1.c6p-31));
    BOOST_TEST(static_cast<std::bfloat16_t>("8.240022e-10"_DFF) == static_cast<std::bfloat16_t>(0x1.c6p-31));
    BOOST_TEST(static_cast<std::bfloat16_t>("1.00390625001"_DDF) == static_cast<std::bfloat16_t>(0x1.02p+0));
    BOOST_TEST(static_cast<std::bfloat16_t>("1.00390625001"_DL) == static_cast<std::bfloat16_t>(0x1.02p+0));
    BOOST_TEST(static_cast<std::bfloat16_t>("1.00390625001"_DLF) == static_cast<std::bfloat16_t>(0x1.02p+0));
    #endif
}

// A 16-bit type gives its own shortest digits, not those of the same value as a float
void test_narrow_shortest()
{
    using namespace boost::decimal;
    using namespace boost::decimal::literals;

    #ifdef BOOST_DECIMAL_HAS_FLOAT16
    BOOST_TEST_EQ(decimal32_t {static_cast<std::float16_t>(0.1F)}, 0.1_DF);
    BOOST_TEST_EQ(decimal64_t {static_cast<std::float16_t>(0.1F)}, 0.1_DD);
    BOOST_TEST_EQ(decimal32_t {static_cast<std::float16_t>(16384)}, 16380_DF);
    BOOST_TEST_EQ(decimal32_t {static_cast<std::float16_t>(0.0078125F)}, 0.007812_DF);
    BOOST_TEST_EQ(decimal32_t {std::numeric_limits<std::float16_t>::min()}, 6.104e-5_DF);
    BOOST_TEST_EQ(decimal32_t {std::numeric_limits<std::float16_t>::denorm_min()}, 6e-8_DF);
    BOOST_TEST_EQ(decimal128_t {static_cast<std::float16_t>(-0.015625F)}, -0.01563_DL);
    BOOST_TEST_EQ(decimal32_t {static_cast<std::float16_t>(-0.0F)}, 0_DF);
    #endif

    #ifdef BOOST_DECIMAL_HAS_BRAINFLOAT16
    BOOST_TEST_EQ(decimal32_t {static_cast<std::bfloat16_t>(0.1F)}, 0.1_DF);
    BOOST_TEST_EQ(decimal32_t {static_cast<std::bfloat16_t>(0.03125F)}, 0.0312_DF);
    BOOST_TEST_EQ(decimal32_t {std::numeric_limits<std::bfloat16_t>::max()}, 3.39e38_DF);
    BOOST_TEST_EQ(decimal32_t {std::numeric_limits<std::bfloat16_t>::denorm_min()}, 9e-41_DF);
    BOOST_TEST_EQ(decimal64_t {-std::numeric_limits<std::bfloat16_t>::denorm_min()}, -9e-41_DD);
    #endif
}

// Each finite 16-bit value, also -0, comes back with the same bits through decimal32_t
template <typename T>
void test_narrow_round_trip()
{
    for (std::uint32_t i {}; i < 65536U; ++i)
    {
        const auto bits {static_cast<std::uint16_t>(i)};
        T x {};
        std::memcpy(&x, &bits, sizeof(bits));
        if (!std::isfinite(static_cast<float>(x)))
        {
            continue;
        }
        const auto back {static_cast<T>(boost::decimal::decimal32_t {x})};
        std::uint16_t back_bits {};
        std::memcpy(&back_bits, &back, sizeof(back_bits));
        BOOST_TEST_EQ(back_bits, bits);
    }
}

// Each positive finite 16-bit value x = m * 2^e gives the nearest of the shortest decimals in its rounding interval.
// compare_exact gives the sign of d * 10^k - n * 2^b, so the test does not use the conversion it checks.
template <typename T>
void test_narrow_nearest_shortest(const int fraction_bits, const int bias)
{
    using boost::decimal::detail::fast_float::compare_exact;
    const auto max_biased {(1 << (15 - fraction_bits)) - 1};
    for (std::uint32_t bits {1U}; bits < (static_cast<std::uint32_t>(max_biased) << fraction_bits); ++bits)
    {
        const auto biased {static_cast<int>(bits >> fraction_bits)};
        const auto fraction {bits & ((1U << fraction_bits) - 1U)};
        const std::uint64_t m {biased == 0 ? fraction : fraction | (1U << fraction_bits)};
        const std::int64_t e {(biased == 0 ? 1 : biased) - bias - fraction_bits};
        const auto narrow_below {fraction == 0U && biased > 1};
        const auto inclusive {m % 2U == 0U};

        const auto narrow_bits {static_cast<std::uint16_t>(bits)};
        T x {};
        std::memcpy(&x, &narrow_bits, sizeof(narrow_bits));
        int k {};
        auto d {static_cast<std::uint64_t>(boost::decimal::frexp10(boost::decimal::decimal32_t {x}, &k))};
        while (d % 10U == 0U)
        {
            d /= 10U;
            ++k;
        }

        // c * 10^j is in the interval between the midpoints with the neighbors of x
        const auto inside {[&](const std::uint64_t c, const std::int64_t j)
        {
            const auto above_low {narrow_below ? compare_exact(j, c, 4U * m - 1U, e - 2) : compare_exact(j, c, 2U * m - 1U, e - 1)};
            const auto below_high {compare_exact(j, c, 2U * m + 1U, e - 1)};
            return (above_low > 0 || (inclusive && above_low == 0)) && (below_high < 0 || (inclusive && below_high == 0));
        }};
        BOOST_TEST(inside(d, k));
        if (d >= 10U)
        {
            BOOST_TEST(!inside(d / 10U, k + 1) && !inside(d / 10U + 1U, k + 1));
        }
        // A neighbor with as many digits is not nearer to x than d
        BOOST_TEST(!inside(d + 1U, k) || compare_exact(k - 1, (2U * d + 1U) * 5U, m, e) >= 0);
        BOOST_TEST(!inside(d - 1U, k) || compare_exact(k - 1, (2U * d - 1U) * 5U, m, e) <= 0);
        // Below 1 * 10^k, the decimal with one digit is 9 * 10^(k - 1)
        BOOST_TEST(d != 1U || !inside(9U, k - 1) || compare_exact(k - 2, 95U, m, e) <= 0);
    }
}

template <typename T>
void test_generic_binary_to_decimal()
{
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<T> dist (1, 2);

    using namespace boost::decimal::detail::ryu;

    std::uint64_t result[5] {};
    generic_computePow5(56, result);
    BOOST_TEST_EQ(result[0], 0);
    BOOST_TEST_EQ(result[1], 5206161169240293376);
    BOOST_TEST_EQ(result[2], 4575641699882439235);

    // https://github.com/boostorg/decimal/issues/242
    BOOST_DECIMAL_IF_CONSTEXPR (std::is_same<T, double>::value)
    {
        // Consider the following integral patterns potentially stored
        // in a given decimal64_t.

        // 0x6A3A25E507BB83D9
        // 0x6A3800D0288A63E9
        // 0x6A380BF1150D7F35
        // 0x6A3B71802C99CB39

        // https://bugs.llvm.org/show_bug.cgi?id=21629
        std::array<std::uint64_t, 4> test_values = {{UINT64_C(0x6A3A25E507BB83D9),
                                                     UINT64_C(0x6A3800D0288A63E9),
                                                     UINT64_C(0x6A380BF1150D7F35),
                                                     UINT64_C(0x6A3B71802C99CB39)}};

        for (const auto bit_values : test_values)
        {
            const std::uint64_t ui = bit_values;

            // And we convert it via direct-memory access to a decimal64_t.

            // In this particular example we assume that this number could
            // be the result of some calculation, initialization or any other
            // valid operations.

            const boost::decimal::decimal64_t dec = *reinterpret_cast<const boost::decimal::decimal64_t*>(&ui);

            // Now let's convert this decimal64_t value to built-in double.
            // ISSUE: This cast seems to result in 0.
            const auto dbl = static_cast<double>(dec);
            const boost::decimal::decimal64_t return_dec(dbl);

            if (!BOOST_TEST(abs(dec - return_dec) <= 1 ))
            {
                // LCOV_EXCL_START
                std::cerr << std::scientific << std::setprecision(std::numeric_limits<boost::decimal::decimal64_t>::digits10)
                          << "       Dec: " << dec
                          << "\n       Dbl: " << dbl
                          << "\nReturn Dec: " << return_dec
                          << "\nDist: " << abs(dec - return_dec) << std::endl;
                // LCOV_EXCL_STOP
            }
        }
    }
}

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

void test_parser()
{
    using namespace boost::decimal::detail;

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<std::uint64_t> dist (1, 2);

    const char* pos_str = "+12345";
    std::uint64_t sig {dist(gen)};
    std::int32_t exp {};
    bool sign {};

    auto res = parser(pos_str, pos_str, sign, sig, exp);
    BOOST_TEST(res.ec == std::errc::invalid_argument);
    res = parser(pos_str, pos_str + std::strlen(pos_str), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc::invalid_argument);

    const char* nan_str = "nan";
    res = parser(nan_str, nan_str + std::strlen(nan_str), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc::not_supported);

    const char* no_trailing = "12345";
    res = parser(no_trailing, no_trailing + std::strlen(no_trailing), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* all_zeros = "0.00000001";
    res = parser(all_zeros, all_zeros + std::strlen(all_zeros), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* big_sig = "123456789012345678901234567890";
    res = parser(big_sig, big_sig + std::strlen(big_sig), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* big_sig_with_frac = "123456789012345678901234567890.123";
    res = parser(big_sig_with_frac, big_sig_with_frac + std::strlen(big_sig_with_frac), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* big_exp = "12345.6789e+1000000";
    res = parser(big_exp, big_exp + std::strlen(big_exp), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST_EQ(res.ptr, big_exp + std::strlen(big_exp));
    BOOST_TEST_EQ(exp, 999996); // 1000000 less the 4 fraction digits

    const char* zeros = "0.00000000";
    res = parser(zeros, zeros + std::strlen(zeros), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* fives = "5.555555555555555555555555555555e+05";
    res = parser(fives, fives + std::strlen(fives), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* sixes = "6.6666666666666666666666666666666e+06";
    res = parser(sixes, sixes + std::strlen(sixes), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* sevens = "7.777777777777777777777777777777e+07";
    res = parser(sevens, sevens + std::strlen(sevens), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* eights = "8.888888888888888888888888888888e+08";
    res = parser(eights, eights + std::strlen(eights), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    const char* nines = "9.99999999999999999999999999999999e+09";
    res = parser(nines, nines + std::strlen(nines), sign, sig, exp);
    BOOST_TEST(res.ec == std::errc());

    // A type wider than 64 bits takes the digits in groups of 19 decimal or 15 hex digits
    using boost::int128::uint128_t;
    uint128_t wide_sig {};
    const char* group_and_four = "12345678901234567890123";
    res = parser(group_and_four, group_and_four + std::strlen(group_and_four), sign, wide_sig, exp);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST(wide_sig == uint128_t{UINT64_C(1234567890123456789)} * 10000U + 123U);

    const uint128_t thirty_eight_digits {uint128_t{UINT64_C(1234567890123456789)} * UINT64_C(10000000000000000000) + UINT64_C(123456789012345678)};
    const char* two_groups = "12345678901234567890123456789012345678";
    res = parser(two_groups, two_groups + std::strlen(two_groups), sign, wide_sig, exp);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST(wide_sig == thirty_eight_digits);

    const char* hex_group = "fffffffffffffffff";
    res = parser(hex_group, hex_group + std::strlen(hex_group), sign, wide_sig, exp, boost::decimal::chars_format::hex);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST(wide_sig == (uint128_t{1U} << 68U) - 1U);

    // Only hex digits can be too many for the type
    const char* hex_33 = "100000000000000000000000000000000";
    res = parser(hex_33, hex_33 + std::strlen(hex_33), sign, wide_sig, exp, boost::decimal::chars_format::hex);
    BOOST_TEST(res.ec == std::errc::result_out_of_range);
    res = parser(hex_33, hex_33 + 17, sign, sig, exp, boost::decimal::chars_format::hex);
    BOOST_TEST(res.ec == std::errc::result_out_of_range);
    res = parser(hex_33, hex_33 + 16, sign, sig, exp, boost::decimal::chars_format::hex);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST_EQ(sig, UINT64_C(1) << 60U);

    // The error comes at the end of the number, after its exponent
    const char* hex_exp = "1.0000000000000000p12x";
    res = parser(hex_exp, hex_exp + std::strlen(hex_exp), sign, sig, exp, boost::decimal::chars_format::hex);
    BOOST_TEST(res.ec == std::errc::result_out_of_range);
    BOOST_TEST_EQ(res.ptr, hex_exp + 21);

    // Eight digits at a time stop at a character just below '0', just above '9' or outside ASCII
    for (const char stop : {'/', ':', 'e', static_cast<char>(0xB9)})
    {
        for (std::size_t k {1}; k < 16U; ++k)
        {
            char digits[] = "1234567890123456";
            digits[k] = stop;
            std::uint64_t expected {};
            for (std::size_t j {}; j < k; ++j)
            {
                expected = expected * 10U + static_cast<std::uint64_t>(digits[j] - '0');
            }

            res = parser(digits, digits + 16, sign, sig, exp, boost::decimal::chars_format::fixed);
            BOOST_TEST(res.ec == std::errc());
            BOOST_TEST_EQ(res.ptr, digits + k);
            BOOST_TEST_EQ(sig, expected);

            res = parser(digits, digits + 16, sign, wide_sig, exp, boost::decimal::chars_format::fixed);
            BOOST_TEST(res.ec == std::errc());
            BOOST_TEST_EQ(res.ptr, digits + k);
            BOOST_TEST(wide_sig == uint128_t{expected});
        }
    }

    const char* fraction_groups = "1.2345678901234567890123456789012345678";
    res = parser(fraction_groups, fraction_groups + std::strlen(fraction_groups), sign, wide_sig, exp);
    BOOST_TEST(res.ec == std::errc());
    BOOST_TEST(wide_sig == thirty_eight_digits);
    BOOST_TEST_EQ(exp, -37);
}

// Each byte value at each place of eight digits gives the same result as one digit at a time.
// The accumulator starts empty, or with digits just below or at the end of a 64-bit group
template <typename Unsigned_Integer>
void test_push_eight(const std::size_t start)
{
    using namespace boost::decimal::detail;

    for (std::size_t k {}; k < 8U; ++k)
    {
        for (int b {}; b < 256; ++b)
        {
            char digits[] = "12345678";
            digits[k] = static_cast<char>(b);

            significand_accumulator<Unsigned_Integer> eight {};
            significand_accumulator<Unsigned_Integer> one {};
            for (std::size_t j {}; j < start; ++j)
            {
                eight.push('9', 10);
                one.push('9', 10);
            }

            const bool accepted {eight.push_eight(digits)};
            BOOST_TEST_EQ(accepted, is_integer_char(digits[k]));
            for (std::size_t j {}; j < 8U && accepted; ++j)
            {
                one.push(digits[j], 10);
            }

            Unsigned_Integer eight_value {};
            Unsigned_Integer one_value {};
            BOOST_TEST(eight.get(eight_value, start + 8U, 10) == std::errc());
            BOOST_TEST(one.get(one_value, start + 8U, 10) == std::errc());
            BOOST_TEST(eight_value == one_value);
        }
    }
}

void test_hex_integer()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "2a";
    auto r1 = boost::decimal::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent, boost::decimal::chars_format::hex);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(significand, UINT64_C(42));
    BOOST_TEST_EQ(exponent, 0);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "-1a3b5c7d9";
    auto r2 = boost::decimal::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::decimal::chars_format::hex);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, 0);
    BOOST_TEST_EQ(significand, UINT64_C(7041566681));

    auto r3 = boost::decimal::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::decimal::chars_format::scientific);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);
}

void test_hex_scientific()
{
    std::uint64_t significand {};
    std::int64_t  exponent {};
    bool sign {};

    const char* val1 = "2ap+5";
    auto r1 = boost::decimal::detail::parser(val1, val1 + std::strlen(val1), sign, significand, exponent, boost::decimal::chars_format::hex);
    BOOST_TEST(r1.ec == std::errc());
    BOOST_TEST_EQ(sign, false);
    BOOST_TEST_EQ(significand, UINT64_C(42));
    BOOST_TEST_EQ(exponent, 5);

    significand = 0;
    exponent = 0;
    sign = false;

    const char* val2 = "-1.3a2bp-10";
    auto r2 = boost::decimal::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::decimal::chars_format::hex);
    BOOST_TEST(r2.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, -14);
    BOOST_TEST_EQ(significand, UINT64_C(80427));

    auto r3 = boost::decimal::detail::parser(val2, val2 + std::strlen(val2), sign, significand, exponent, boost::decimal::chars_format::scientific);
    BOOST_TEST(r3.ec == std::errc::invalid_argument);

    const char* val4 = "-1.3A2BP-10";
    auto r4 = boost::decimal::detail::parser(val4, val4 + std::strlen(val4), sign, significand, exponent, boost::decimal::chars_format::hex);
    BOOST_TEST(r4.ec == std::errc());
    BOOST_TEST_EQ(sign, true);
    BOOST_TEST_EQ(exponent, -14);
    BOOST_TEST_EQ(significand, UINT64_C(80427));
}

void test_from_chars()
{
    using namespace boost::decimal::detail;

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<std::uint64_t> dist (1, 2);

    const char* pos_str = "+12345";
    std::uint64_t sig {dist(gen)};

    auto res = from_chars_integer_impl<std::uint64_t, std::uint64_t>(pos_str, pos_str, sig, 10);
    BOOST_TEST(res.ec == std::errc::invalid_argument);
    res = from_chars_integer_impl<std::uint64_t, std::uint64_t>(pos_str, pos_str + std::strlen(pos_str), sig, 10);
    BOOST_TEST(res.ec == std::errc::invalid_argument);

    const char* neg_str = "-12345";
    res = from_chars_integer_impl<std::uint64_t, std::uint64_t>(neg_str, neg_str + std::strlen(neg_str), sig, 10);
    BOOST_TEST(res.ec == std::errc::invalid_argument);

    auto signed_sig {static_cast<std::int64_t>(dist(gen))};
    res = from_chars_integer_impl<std::int64_t, std::uint64_t>(neg_str, neg_str + 1, signed_sig, 10);
    BOOST_TEST(res.ec == std::errc::invalid_argument);

    res = from_chars_integer_impl<std::int64_t, std::uint64_t>(neg_str, neg_str + std::strlen(neg_str), signed_sig, 10);
    BOOST_TEST(res.ec == std::errc());

    const char* big_sig = "-123456789123456789123456789123456789";
    res = from_chars_integer_impl<std::int64_t, std::uint64_t>(big_sig, big_sig + std::strlen(big_sig), signed_sig, 10);
    BOOST_TEST(res.ec == std::errc::result_out_of_range);
}

// A negative zero or infinity stays negative in the binary type
template <typename Dec, typename Float>
void test_signed_zero_and_infinity()
{
    const Dec zero {0};
    const Dec inf {std::numeric_limits<Dec>::infinity()};

    BOOST_TEST(!std::signbit(static_cast<Float>(zero)));
    BOOST_TEST(std::signbit(static_cast<Float>(-zero)));
    BOOST_TEST_EQ(static_cast<Float>(inf), std::numeric_limits<Float>::infinity());
    BOOST_TEST_EQ(static_cast<Float>(-inf), -std::numeric_limits<Float>::infinity());
}

template <typename Float>
void test_signed_zero_and_infinity_to()
{
    using namespace boost::decimal;

    test_signed_zero_and_infinity<decimal32_t, Float>();
    test_signed_zero_and_infinity<decimal64_t, Float>();
    test_signed_zero_and_infinity<decimal128_t, Float>();
    test_signed_zero_and_infinity<decimal_fast32_t, Float>();
    test_signed_zero_and_infinity<decimal_fast64_t, Float>();
    test_signed_zero_and_infinity<decimal_fast128_t, Float>();
}

int main()
{
    test_compute_float32();
    test_compute_float64();
    test_compute_float80_128();
    test_generic_binary_to_decimal<float>();
    test_generic_binary_to_decimal<double>();

    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    test_generic_binary_to_decimal<long double>();
    #endif

    #if BOOST_DECIMAL_LDBL_BITS == 80 || BOOST_DECIMAL_LDBL_BITS == 128
    test_ryu<long double>();
    #endif

    test_shortest_tie();
    test_narrow_midpoint();
    test_narrow_shortest();
    #ifdef BOOST_DECIMAL_HAS_FLOAT16
    test_narrow_round_trip<std::float16_t>();
    test_narrow_nearest_shortest<std::float16_t>(10, 15);
    #endif
    #ifdef BOOST_DECIMAL_HAS_BRAINFLOAT16
    test_narrow_round_trip<std::bfloat16_t>();
    test_narrow_nearest_shortest<std::bfloat16_t>(7, 127);
    #endif

    test_parser();
    test_push_eight<std::uint64_t>(0U);
    test_push_eight<std::uint64_t>(11U);
    test_push_eight<boost::int128::uint128_t>(0U);
    test_push_eight<boost::int128::uint128_t>(11U);
    test_push_eight<boost::int128::uint128_t>(12U);
    test_hex_integer();
    test_hex_scientific();

    test_from_chars();

    test_signed_zero_and_infinity_to<float>();
    test_signed_zero_and_infinity_to<double>();
    #ifndef BOOST_DECIMAL_UNSUPPORTED_LONG_DOUBLE
    test_signed_zero_and_infinity_to<long double>();
    #endif

    return boost::report_errors();
}
