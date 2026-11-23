// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// erf_error <count> <seed> <lo> <hi>: largest difference of the two-word and three-word digits of decimal128_t erf and
// erfc for x in [lo, hi] and x >= 1, in units of the two-word digits. See README.md in this directory.

#include <boost/decimal.hpp>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
using namespace boost::decimal;
using boost::int128::uint128_t;
// |a 10^ea - b 10^eb| in units of the smaller exponent for |ea - eb| <= 1, as test_first_try_bounds.cpp compares them:
// c of the larger one and f of the smaller one give 10 c - f = 10 (c - f / 10) - f % 10
static uint128_t units_apart(const uint128_t a, const int ea, const uint128_t b, const int eb)
{
    if (ea == eb)
        return a > b ? a - b : b - a;
    const uint128_t c = ea > eb ? a : b, f = ea > eb ? b : a;
    const uint128_t q = f / 10U, r = f % 10U, g = c > q ? c - q : q - c;
    if (g > (std::numeric_limits<uint128_t>::max)() / 16U)
        return (std::numeric_limits<uint128_t>::max)();
    return c >= q ? (g * 10U > r ? g * 10U - r : r - g * 10U) : g * 10U + r;
}
int main(int argc, char **v)
{
    if (argc != 5)
    {
        std::puts("usage: erf_error count seed lo hi");
        return 2;
    }
    long n = std::atol(v[1]);
    std::mt19937_64 g(std::atol(v[2]));
    double lo = std::atof(v[3]), hi = std::atof(v[4]);
    std::uniform_real_distribution<double> u(lo, hi);
    double mx[4]{};
    long cnt[4]{};
    for (long i = 0; i < n; ++i)
    {
        char b[64];
        std::snprintf(b, sizeof b, "%.6f", u(g));
        std::string s = b;
        for (int k = 0; k < 28; ++k)
            s += char('0' + g() % 10);
        decimal128_t x{};
        from_chars(s.data(), s.data() + s.size(), x);
        int e{};
        const auto sig{static_cast<uint128_t>(frexp10(x, &e))};
        const int k{-e};
        const auto whole{static_cast<std::uint64_t>(detail::exp_detail::exp_whole(sig, k))};
        if (whole >= 120U || whole < 1U)
            continue;
        for (int c = 0; c < 4; ++c)
        {
            const bool comp = c & 1, neg = c & 2;
            auto r2{detail::erf_detail::erf_above_one<decimal128_t, 2>(sig, k, whole, comp && neg, neg)};
            if (comp && !neg)
            {
                // Positive erfc takes its own path in the first try
                const auto d2{detail::erf_detail::erfc_above_one_digits2(sig, k, whole)};
                r2.digits = d2.digits;
                r2.exp10 = d2.exp10;
            }
            const auto r3{detail::erf_detail::erf_above_one<decimal128_t, 3>(sig, k, whole, comp, neg)};
            if ((r2.kind != detail::exp_detail::exp_kind::digits && (!comp || neg)) ||
                r3.kind != detail::exp_detail::exp_kind::digits)
                continue;
            const int de = r2.exp10 - r3.exp10;
            if (de < -1 || de > 1)
            {
                std::printf("exp10 mismatch %d\n", de);
                continue;
            }
            const uint128_t diff = units_apart(r2.digits, r2.exp10, r3.digits, r3.exp10);
            const double dd = static_cast<double>(diff);
            if (dd > mx[c])
            {
                mx[c] = dd;
                if (dd > 100)
                    std::printf("c%d %s diff %.0f\n", c, s.c_str(), dd);
            }
            ++cnt[c];
        }
    }
    for (int c = 0; c < 4; ++c)
        std::printf("comp %d neg %d n %ld max %.0f\n", c & 1, (c >> 1) & 1, cnt[c], mx[c]);
}
