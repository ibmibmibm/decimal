// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// zeta_error <n>: error of the first try of riemann_zeta for x > 0 against MPFR, in units of its digits. The library
// allows 24 units of 2^-62 in one word, and 44 units of the digits in two words, or 116 within 1/2 of the pole.
// x stops where the library returns a value just above 1. See README.md in this directory.

#include <boost/decimal.hpp>
#include <mpfr.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
using namespace boost::decimal;
namespace rz = boost::decimal::detail::riemann_zeta_detail;
using u128 = boost::int128::uint128_t;
static std::string u128s(u128 v)
{
    std::string s;
    do
    {
        s.insert(s.begin(), char('0' + static_cast<int>((v % 10U).low)));
        v /= 10U;
    } while (v != 0U);
    return s;
}
template <typename T, bool Precise> void run(const char *tn, long n)
{
    constexpr int p = std::numeric_limits<T>::digits10;
    struct R
    {
        const char *name;
        int kind;
        double a, b;
    } rs[] = {{"(1e-9,0.01)", 3, 1, 9},   {"(0.01,0.5)", 0, 0.01, 0.5},
              {"(0.5,1.5)", 0, 0.5, 1.5}, {"1+-10^-u", 1, 1, double(p - 1)},
              {"(1.5,4)", 0, 1.5, 4},     {"(4,16)", 0, 4, 16},
              {"(16,64)", 0, 16, 64},     {"(64,128)", 0, 64, 128},
              {"near ints", 2, 2, 128}};
    std::mt19937_64 g(11);
    std::uniform_real_distribution<double> d(0, 1);
    mpfr_t xm, zm, dm, t, ul;
    mpfr_inits2(300, xm, zm, dm, t, ul, (mpfr_ptr)0);
    for (auto &r : rs)
    {
        double worst = 0;
        std::string wx;
        for (long i = 0; i < n; ++i)
        {
            double v =
                r.kind == 0   ? r.a + (r.b - r.a) * d(g)
                : r.kind == 3 ? std::pow(10.0, -(r.a + (r.b - r.a) * d(g)))
                : r.kind == 1
                    ? 1 + ((g() & 1) ? 1 : -1) * std::pow(10.0, -(r.a + (r.b - r.a) * d(g))) * (1 + d(g))
                    : std::floor(r.a + (r.b - r.a) * d(g)) + ((g() & 1) ? 1 : -1) * std::pow(10.0, -(1 + 14 * d(g)));
            char b[64];
            std::snprintf(b, sizeof b, "%.17e", v);
            T x;
            from_chars(b, b + std::strlen(b), x);
            x += T{static_cast<long long>(g() % 1000000000), -30};
            if (!(x > T{0}) || x == T{1} || x > T{p <= 7 ? 21 : p <= 16 ? 51 : 111})
                continue;
            int e;
            const auto sig = static_cast<u128>(frexp10(x, &e));
            const bool small{x > T{5, -1} && x < T{15, -1}};
            const bool neg{x < T{1}};
            const auto d2 =
                rz::zeta_positive_digits<std::numeric_limits<T>::digits10 == 34 ? 2 : 1, Precise>(sig, e, small, neg);
            mpfr_set_str(xm, (u128s(sig) + "e" + std::to_string(e)).c_str(), 10, MPFR_RNDN);
            mpfr_zeta(zm, xm, MPFR_RNDN);
            mpfr_set_str(dm, (std::string(d2.neg ? "-" : "") + u128s(d2.digits)).c_str(), 10, MPFR_RNDN);
            mpfr_set_si(ul, 10, MPFR_RNDN);
            mpfr_t pw;
            mpfr_init2(pw, 300);
            mpfr_set_si(pw, d2.exp10, MPFR_RNDN);
            mpfr_pow(ul, ul, pw, MPFR_RNDN);
            mpfr_clear(pw);
            mpfr_div(t, zm, ul, MPFR_RNDN);
            mpfr_sub(t, dm, t, MPFR_RNDN);
            mpfr_abs(t, t, MPFR_RNDN);
            const double err = mpfr_get_d(t, MPFR_RNDN);
            if (err > worst)
            {
                worst = err;
                wx = b;
            }
        }
        std::printf("%s P=%d %-12s digits max %.2f at %s\n", tn, Precise ? 1 : 0, r.name, worst, wx.c_str());
        fflush(stdout);
    }
    mpfr_clears(xm, zm, dm, t, ul, (mpfr_ptr)0);
}
int main(int argc, char **argv)
{
    const long n = argc > 1 ? std::atol(argv[1]) : 8000;
    run<decimal32_t, false>("d32", n);
    run<decimal64_t, true>("d64", n);
    run<decimal128_t, true>("d128", n);
}
