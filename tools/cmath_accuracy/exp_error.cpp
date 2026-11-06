// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// exp_error <n>: the one-word digits d1 of exp, exp2, sinh and cosh for decimal64_t against three words, in units
// of d1.
// See README.md in this directory.

#include <boost/decimal.hpp>
#include <cstdio>
#include <cmath>
#include <random>
using namespace boost::decimal;
using namespace boost::decimal::detail;
using namespace boost::decimal::detail::exp_detail;
using boost::int128::uint128_t;

template <int W, bool Base2> static uint128_t digits_of(const uint128_t &sig, int k, bool neg, int &ex)
{
    const auto whole{static_cast<std::uint64_t>(exp_whole(sig, k))};
    auto y{exp_log10(whole, exp_frac<W>(sig, k, whole), Base2 ? fx_table::log10_2 : fx_table::log10_e)};
    if (neg)
    {
        y.exp10 = -y.exp10 - 1;
        y.f = fx_sub(fx_one<W>(), y.f);
    }
    const auto t{exp10_frac(y.f)};
    ex = y.exp10;
    return fx_scale(pow10(static_cast<uint128_t>(W == 1 ? 18 : 37)) << (t.a + 1U), fx_shr(t.m, 1U));
}

int main(int argc, char **argv)
{
    const long n{argc > 1 ? std::atol(argv[1]) : 1000000};
    std::mt19937_64 rng{7};
    for (int b2 = 0; b2 < 2; ++b2)
        for (double hi : {1e-15, 1e-3, 1.0, 50.0, 880.0})
        {
            double maxe{};
            for (long i{}; i < n; ++i)
            {
                std::uniform_real_distribution<double> d{std::log10(hi) - 3, std::log10(hi)};
                auto x{decimal64_t{std::pow(10.0, d(rng))}};
                int e{};
                auto s{frexp10(x, &e)};
                s = s - s % 100000U + rng() % 100000U;
                const bool neg{(rng() & 1U) != 0U};
                int e1{}, e3{};
                // The first try runs only for these exponents
                if (-e < 12 || -e > 32)
                    continue;
                const uint128_t d1{b2 ? exp_word_digits<true, false>(s, -e, neg, e1) : exp_word_digits<false, false>(s, -e, neg, e1)};
                const auto d3{b2 ? digits_of<3, true>(s, -e, neg, e3) : digits_of<3, false>(s, -e, neg, e3)};
                const auto q{d3 / pow10(static_cast<uint128_t>(19))};
                const auto rr{d3 - q * pow10(static_cast<uint128_t>(19))};
                const double frac{static_cast<double>(rr) / 1e19};
                const double err{std::fabs(static_cast<double>(static_cast<std::int64_t>(
                                               static_cast<std::uint64_t>(d1) - static_cast<std::uint64_t>(q))) -
                                           frac)};
                // when the two differ on exp10 skip
                if (e1 != e3)
                    continue;
                if (err > maxe)
                    maxe = err;
            }
            std::printf("%s |x| in [%g, %g): max %.2f units of d1\n", b2 ? "exp2" : "exp", hi / 1000, hi, maxe);
        }

    // sinh for |x| >= 1 and cosh for |x| >= 0.1: the last range is next above ln(20), where sinh is next below 10
    for (int ch = 0; ch < 2; ++ch)
        for (double hi : {1.0, 3.0, 23.0, 880.0, 0.0})
        {
            if (!ch && hi == 1.0)
                continue;
            double maxe{};
            for (long i{}; i < n; ++i)
            {
                int e{-15};
                std::uint64_t s{UINT64_C(2995732273553991) + rng() % UINT64_C(2500000000000)};
                if (hi > 0.0)
                {
                    std::uniform_real_distribution<double> d{hi == 1.0 ? 0.1 : hi == 3.0 ? 1.0 : hi == 23.0 ? 3.0 : 23.0, hi};
                    s = frexp10(decimal64_t{d(rng)}, &e);
                    s = s - s % 100000U + rng() % 100000U;
                }
                const auto whole{static_cast<std::uint64_t>(exp_whole(uint128_t{s}, -e))};
                if (!ch && whole == 0U)
                    continue;
                int e1{};
                const uint128_t d1{sinh_cosh_word_digits<false>(s, -e, ch != 0, e1)};
                const auto r{sinh_cosh_parts<3>(whole, exp_frac<3>(uint128_t{s}, -e, whole), ch != 0)};
                const auto d3{fx_scale(pow10(static_cast<uint128_t>(37)) * 4U, r.f)};
                const auto scale{pow10(static_cast<uint128_t>(e1 - (r.exp10 - 37)))};
                const auto a{d1 * scale};
                const double err{static_cast<double>(a > d3 ? a - d3 : d3 - a) / static_cast<double>(scale)};
                if (err > maxe)
                    maxe = err;
            }
            if (hi > 0.0)
                std::printf("%s |x| below %g: max %.2f units of d1\n", ch ? "cosh" : "sinh", hi, maxe);
            else
                std::printf("%s |x| next above ln(20): max %.2f units of d1\n", ch ? "cosh" : "sinh", maxe);
        }
}
