// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Checks the error bounds of the first tries of exp, exp2, erf and erfc against the same code in more words.
// A bound that is too small gives a wrong result only next to a rounding boundary, so the result alone rarely shows it.

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <cstdint>

using namespace boost::decimal;
using namespace boost::decimal::detail;
using boost::int128::uint128_t;

namespace {

constexpr int count {20000};

std::uint64_t g_state {UINT64_C(0x9E3779B97F4A7C15)};

auto next_word() -> std::uint64_t
{
    g_state ^= g_state << 13U;
    g_state ^= g_state >> 7U;
    g_state ^= g_state << 17U;
    return g_state;
}

// A significand of p digits
auto next_sig(const int p) -> uint128_t
{
    uint128_t sig {next_word() % 9U + 1U};
    for (int i {1}; i < p; ++i)
    {
        sig = sig * 10U + next_word() % 10U;
    }
    return sig;
}

auto distance(const uint128_t& a, const uint128_t& b) -> uint128_t
{
    return a > b ? a - b : b - a;
}

// The first try in the given words has a bound of d * 2^-bits + base units of d, and three words add 8 units of theirs
template <int words, bool Precise, bool Short, bool Base2>
void test_exp(const int p, const unsigned bits, const unsigned base)
{
    const auto& c {Base2 ? fx_table::log10_2 : fx_table::log10_e};
    constexpr int reference {3};
    const auto scale {pow10(static_cast<uint128_t>(fx_scale10(reference) - fx_scale10(words)))};

    for (int k_shift {-2}; k_shift < 22; ++k_shift)
    {
        for (int i {0}; i < count / 24; ++i)
        {
            const auto sig {next_sig(p)};
            const int k {p + k_shift};
            const auto whole {static_cast<std::uint64_t>(exp_detail::exp_whole(sig, k))};
            if (whole >= 880U)
            {
                continue;
            }
            const bool neg {(next_word() & 1U) != 0U};
            const auto y1 {exp_detail::exp_log10_signed<words>(sig, k, whole, c, neg)};
            const auto y3 {exp_detail::exp_log10_signed<reference>(sig, k, whole, c, neg)};
            if (y1.exp10 != y3.exp10)
            {
                continue;
            }
            const auto d1 {exp_detail::exp10_digits(exp_detail::exp10_frac<words, Precise, Short>(y1.f), fx_scale10(words))};
            const auto d3 {exp_detail::exp10_digits(exp_detail::exp10_frac(y3.f), fx_scale10(reference))};
            BOOST_TEST(distance(d1 * scale, d3) <= ((d1 >> bits) + base) * scale + 8U);
        }
    }
}

void test_erf_above_one()
{
    using T = decimal128_t;
    for (int i {0}; i < count; ++i)
    {
        const int p {34};
        const int k {p - 1 - static_cast<int>(next_word() % 2U)};
        const auto sig {next_sig(p)};
        const auto whole {static_cast<std::uint64_t>(exp_detail::exp_whole(sig, k))};
        if (whole < 1U || whole >= 120U)
        {
            continue;
        }

        for (int kind {0}; kind < 3; ++kind)
        {
            const bool complement {kind != 0};
            const bool neg {kind == 1 || (kind == 0 && (next_word() & 1U) != 0U)};
            const auto r3 {erf_detail::erf_above_one<T, 3>(sig, k, whole, complement, neg)};
            if (r3.kind != exp_detail::exp_kind::digits)
            {
                continue;
            }
            uint128_t d2 {};
            int e2 {};
            std::uint64_t limit {};
            if (kind == 2)
            {
                const auto r2 {erf_detail::erfc_above_one_digits2(sig, k, whole)};
                d2 = r2.digits;
                e2 = r2.exp10;
                limit = static_cast<std::uint64_t>(((r2.digits >> 64U) * (whole < 12U ? 32U : 67U)) >> 62U) + 2U;
            }
            else
            {
                const auto r2 {erf_detail::erf_above_one<T, 2>(sig, k, whole, complement, neg)};
                if (r2.kind != exp_detail::exp_kind::digits)
                {
                    continue;
                }
                d2 = r2.digits;
                e2 = r2.exp10;
                limit = complement ? 3U : 10U;
            }

            // The two values can lie on the two sides of a power of ten
            auto d3 {static_cast<uint128_t>(r3.digits)};
            if (e2 == r3.exp10 + 1)
            {
                d3 /= 10U;
            }
            else if (e2 == r3.exp10 - 1)
            {
                d2 /= 10U;
            }
            BOOST_TEST(e2 >= r3.exp10 - 1 && e2 <= r3.exp10 + 1);
            BOOST_TEST(distance(d2, d3) <= limit);
        }
    }
}

void test_erf_below_one()
{
    using T = decimal128_t;
    for (int i {0}; i < count; ++i)
    {
        const int p {34};
        const int e {-p - static_cast<int>(next_word() % 4U)};
        const auto sig {next_sig(p)};
        const bool complement {(next_word() & 1U) != 0U};
        const bool neg {(next_word() & 1U) != 0U};
        const auto r2 {erf_detail::erf_small_x<T, 2>(sig, e, complement, neg)};
        const auto r3 {erf_detail::erf_small_x<T, 3>(sig, e, complement, neg)};
        auto d2 {static_cast<uint128_t>(r2.digits)};
        auto d3 {static_cast<uint128_t>(r3.digits)};
        const std::uint64_t limit {complement ? UINT64_C(250) : static_cast<std::uint64_t>(d2 >> 118U) + 2U};
        if (r2.exp10 == r3.exp10 + 1)
        {
            d3 /= 10U;
        }
        else if (r2.exp10 == r3.exp10 - 1)
        {
            d2 /= 10U;
        }
        BOOST_TEST(distance(d2, d3) <= limit);
    }
}

} // namespace

int main()
{
    test_exp<1, false, true, false>(7, 37U, 2U);
    test_exp<1, false, true, true>(7, 37U, 2U);
    test_exp<1, false, false, false>(16, 58U, 3U);
    test_exp<1, false, false, true>(16, 58U, 3U);
    test_exp<2, true, false, false>(34, 119U, 2U);
    test_exp<2, true, false, true>(34, 119U, 2U);
    test_erf_above_one();
    test_erf_below_one();

    return boost::report_errors();
}
