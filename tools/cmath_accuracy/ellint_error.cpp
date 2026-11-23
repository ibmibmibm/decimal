// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// ellint_error <n>: error of ellint_fit, the one-word first try of K and E, against the AGM in MPFR. It prints units of
// 2^-62 and the largest ratio to the bound that the library takes. See README.md in this directory.

#include <boost/decimal.hpp>
#include <mpfr.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
using namespace boost::decimal;
namespace ed = boost::decimal::detail::ellint_detail;
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
// K = pi / (2 a_inf), E = K (1 - sum 2^(n-1) c_n^2) from a = 1, b = k', c_0 = k
static void ke(mpfr_t K, mpfr_t E, mpfr_t k)
{
    mpfr_t a, b, c, t, s, p;
    mpfr_inits2(500, a, b, c, t, s, p, (mpfr_ptr)0);
    mpfr_set_ui(a, 1, MPFR_RNDN);
    mpfr_sqr(t, k, MPFR_RNDN);
    mpfr_ui_sub(b, 1, t, MPFR_RNDN);
    mpfr_sqrt(b, b, MPFR_RNDN);
    mpfr_div_2ui(s, t, 1, MPFR_RNDN);
    mpfr_set_ui(p, 1, MPFR_RNDN);
    for (int i = 0; i < 40; ++i)
    {
        mpfr_sub(c, a, b, MPFR_RNDN);
        mpfr_div_2ui(c, c, 1, MPFR_RNDN);
        mpfr_add(t, a, b, MPFR_RNDN);
        mpfr_mul(b, a, b, MPFR_RNDN);
        mpfr_sqrt(b, b, MPFR_RNDN);
        mpfr_div_2ui(a, t, 1, MPFR_RNDN);
        mpfr_sqr(t, c, MPFR_RNDN);
        mpfr_mul(t, t, p, MPFR_RNDN);
        mpfr_add(s, s, t, MPFR_RNDN);
        mpfr_mul_2ui(p, p, 1, MPFR_RNDN);
    }
    mpfr_const_pi(K, MPFR_RNDN);
    mpfr_div(K, K, a, MPFR_RNDN);
    mpfr_div_2ui(K, K, 1, MPFR_RNDN);
    mpfr_ui_sub(t, 1, s, MPFR_RNDN);
    mpfr_mul(E, K, t, MPFR_RNDN);
    mpfr_clears(a, b, c, t, s, p, (mpfr_ptr)0);
}
template <typename T> void run(const char *tn, long n)
{
    constexpr int p = std::numeric_limits<T>::digits10;
    std::mt19937_64 g(3);
    std::uniform_real_distribution<double> d(0, 1);
    mpfr_t km, K, E, dm, t;
    mpfr_inits2(500, km, K, E, dm, t, (mpfr_ptr)0);
    const char *names[] = {"(0,0.866)", "(0.866,0.99)", "1-10^-u"};
    for (int r = 0; r < 3; ++r)
        for (int kind = 0; kind < 2; ++kind)
        {
            double worst = 0, ratio = 0;
            std::string wx;
            for (long i = 0; i < n; ++i)
            {
                const double v = r == 0   ? 0.866 * d(g)
                                 : r == 1 ? 0.866 + 0.124 * d(g)
                                          : 1 - std::pow(10.0, -(2 + (p - 2) * d(g)));
                char b[64];
                std::snprintf(b, sizeof b, "%.*e", p - 1, v);
                T x;
                from_chars(b, b + std::strlen(b), x);
                if (!(x < T{1}))
                    continue;
                int e;
                const auto sig = static_cast<u128>(frexp10(x, &e));
                std::uint64_t fd_err{};
                int fd_exp10{};
                const auto fd_digits = ed::ellint_fit(sig, -e, kind ? ed::ellint_kind::second : ed::ellint_kind::first, fd_err, fd_exp10);
                mpfr_set_str(km, (u128s(sig) + "e" + std::to_string(e)).c_str(), 10, MPFR_RNDN);
                ke(K, E, km);
                mpfr_set_str(dm, (u128s(fd_digits) + "e" + std::to_string(fd_exp10)).c_str(), 10, MPFR_RNDN);
                mpfr_sub(t, dm, kind ? E : K, MPFR_RNDN);
                mpfr_abs(t, t, MPFR_RNDN);
                mpfr_mul_2si(t, t, 62, MPFR_RNDN);
                const double err = mpfr_get_d(t, MPFR_RNDN);
                if (err > worst)
                {
                    worst = err;
                    wx = b;
                }
                ratio = std::max(ratio, err / static_cast<double>(fd_err));
            }
            std::printf("%s %s %-13s max %.2f units at %s, max err/bound %.3f\n", tn, kind ? "E" : "K", names[r], worst,
                        wx.c_str(), ratio);
        }
    mpfr_clears(km, K, E, dm, t, (mpfr_ptr)0);
}
int main(int argc, char **argv)
{
    const long n = argc > 1 ? std::atol(argv[1]) : 20000;
    run<decimal32_t>("d32", n);
    run<decimal64_t>("d64", n);
}
