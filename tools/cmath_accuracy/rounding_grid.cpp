// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Compares the <cmath> functions against MPFR at fixed grids of arguments, for one type and in all five rounding modes.
// See README.md in this directory.

#include <boost/decimal.hpp>
#include <mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <array>
using namespace boost::decimal;
static int g_precision = 500;
struct big_float
{
    mpfr_t v;
    big_float(int p = g_precision)
    {
        mpfr_init2(v, p);
    }
    ~big_float()
    {
        mpfr_clear(v);
    }
    big_float(const big_float &) = delete;
};
template <class U> std::string ustr(U s)
{
    if (s == 0)
        return "0";
    std::string r;
    while (s != 0)
    {
        r.insert(r.begin(), char('0' + static_cast<int>(s % 10)));
        s /= 10;
    }
    return r;
}
template <class T> void to_mpfr(T x, mpfr_t o)
{
    int e;
    auto s = frexp10(x, &e);
    std::string st = ustr(s) + "e" + std::to_string(e);
    mpfr_set_str(o, st.c_str(), 10, MPFR_RNDN);
    if (signbit(x))
        mpfr_neg(o, o, MPFR_RNDN);
}
template <class T> std::string show(T x)
{
    if (isnan(x))
        return "nan";
    if (isinf(x))
        return signbit(x) ? "-inf" : "inf";
    int e;
    auto s = frexp10(x, &e);
    return std::string(signbit(x) ? "-" : "") + ustr(s) + "e" + std::to_string(e);
}
template <class T> T list_argument(const std::string &s)
{
    T x{};
    const bool neg = s[0] == '-';
    from_chars(s.data() + (neg ? 1 : 0), s.data() + s.size(), x);
    return neg ? -x : x;
}
#define MPFR1(f) [](mpfr_t o, mpfr_t a) { f(o, a, MPFR_RNDN); }
#define DEC1(f) [](T x) { return f(x); }
#ifdef BOUND_MARGIN
#include "bound_margin.hpp"
#else
static void margin_events(mpfr_t, const char *, const std::string &)
{
}
#endif
// With GRID_REC set to a file, the file gets a line for each check: name|arguments|mode|result
static FILE *g_rec = std::getenv("GRID_REC") ? std::fopen(std::getenv("GRID_REC"), "w") : nullptr;
static long g_bad_total;
const rounding_mode g_modes[] = {rounding_mode::fe_dec_to_nearest, rounding_mode::fe_dec_upward,
                                 rounding_mode::fe_dec_downward, rounding_mode::fe_dec_to_nearest_from_zero,
                                 rounding_mode::fe_dec_toward_zero};
static bool g_strict = true, g_fixed_tolerance = false;
template <class T> struct is_fast
{
    static constexpr bool v = false;
};
template <> struct is_fast<decimal_fast32_t>
{
    static constexpr bool v = true;
};
template <> struct is_fast<decimal_fast64_t>
{
    static constexpr bool v = true;
};
template <> struct is_fast<decimal_fast128_t>
{
    static constexpr bool v = true;
};
static constexpr int check_precision = 60000;
// Checks r against the value t in mode index mode, with zero, subnormal and overflow results.
// The direction dir is 0 for nearest, 1 for up and 2 for down.
static bool g_overflow, g_underflow;
template <class T> bool check(T r, mpfr_t t, int mode, const char *name, const std::string &arg, long &bad)
{
    if (g_rec)
        std::fprintf(g_rec, "%s|%s|%d|%s\n", name, arg.c_str(), mode, show(r).c_str());
    margin_events(t, name, arg);
    static big_float d(check_precision), u(check_precision), mx(check_precision), mn(check_precision),
        a0(check_precision);
    auto say = [&](const char *w, double du)
    {
        if (bad++ < 8)
            printf("  %s(%s) mode %d %s err %.6f got %s\n", name, arg.c_str(), mode, w, du, show(r).c_str());
        return false;
    };
    if (mpfr_nan_p(t))
        return isnan(r) ? true : say("NaN expected", 0);
    if (isnan(r))
        return say("unexpected NaN", 0);
    const bool neg = mpfr_signbit(t);
    // effective direction: 0 nearest, 1 up, 2 down
    const int dir = mode == 1 ? 1 : mode == 2 ? 2 : mode == 4 ? (neg ? 1 : 2) : 0;
    if (mpfr_zero_p(t) && g_underflow)
    {
        if (is_fast<T>::v)
            return true;
        const T dm = std::numeric_limits<T>::denorm_min();
        const bool ok = dir == 0 || (dir == 2 && !neg) || (dir == 1 && neg) ? (r == T{0}) : abs(r) == dm;
        return ok && signbit(r) == neg ? true : say("underflow", 0);
    }
    if (mpfr_zero_p(t))
    {
        if (r != T{0})
            return say("zero expected", 0);
        return signbit(r) == neg ? true : say("zero sign", 0);
    }
    to_mpfr(std::numeric_limits<T>::max(), mx.v);
    mpfr_abs(u.v, t, MPFR_RNDN);
    if (mpfr_inf_p(t) && !g_overflow)
    {
        return (isinf(r) && signbit(r) == neg) ? true : say("pole expected", 0);
    }
    if (mpfr_inf_p(t) || mpfr_cmp(u.v, mx.v) > 0)
    {
        if (is_fast<T>::v)
            return true;
        const T mxv = std::numeric_limits<T>::max();
        const bool toinf = dir == 0 || (dir == 1 && !neg) || (dir == 2 && neg);
        if (toinf)
        {
            if (isinf(r) && signbit(r) == neg)
                return true;
            if (!mpfr_inf_p(t) && dir == 0 && abs(r) == mxv)
                return true;
            return say("inf expected", 0);
        }
        return r == (neg ? -mxv : mxv) ? true : say("max expected", 0);
    }
    if (isinf(r))
        return say("inf unexpected", 0);
    int e;
    auto s = frexp10(r, &e);
    int etiny;
    {
        auto q = frexp10(std::numeric_limits<T>::denorm_min(), &etiny);
        int dg = 0;
        for (auto w = q; w != 0; w /= 10)
            ++dg;
        etiny += dg - 1;
    }
    to_mpfr(std::numeric_limits<T>::min(), mn.v);
    const bool sub = mpfr_cmp(u.v, mn.v) < 0 || r == T{0};
    if (sub)
    {
        if (is_fast<T>::v)
            return true;
        e = etiny;
    }
    to_mpfr(r, d.v);
    mpfr_sub(d.v, d.v, t, MPFR_RNDN);
    mpfr_ui_pow_ui(u.v, 10, e < 0 ? -e : e, MPFR_RNDN);
    if (e < 0)
        mpfr_mul(d.v, d.v, u.v, MPFR_RNDN);
    else
        mpfr_div(d.v, d.v, u.v, MPFR_RNDN);
    double du = mpfr_get_d(d.v, MPFR_RNDN);
    static big_float tl(check_precision), w1(check_precision), w2(check_precision);
    bool ok;
    if (!g_strict || g_fixed_tolerance)
    {
        const double tol = 1e-9;
        if (dir == 1)
            ok = du >= -tol && du < 1 + tol;
        else if (dir == 2)
            ok = du <= tol && du > -1 - tol;
        else
            ok = du >= -0.5 - tol && du <= 0.5 + tol;
        return ok ? true : say("misround", du);
    }
    mpfr_set_ui_2exp(tl.v, 1, -(g_precision - 300), MPFR_RNDN);
    if (sub && abs(r) == std::numeric_limits<T>::denorm_min() && std::fabs(std::fabs(du) - 1) < 1e-9)
        return true;
    if (dir == 1)
    {
        mpfr_neg(w1.v, tl.v, MPFR_RNDN);
        mpfr_ui_sub(w2.v, 1, tl.v, MPFR_RNDN);
        ok = mpfr_cmp(d.v, w1.v) >= 0 && mpfr_cmp(d.v, w2.v) < 0;
    }
    else if (dir == 2)
    {
        mpfr_sub_ui(w2.v, tl.v, 1, MPFR_RNDN);
        ok = mpfr_cmp(d.v, tl.v) <= 0 && mpfr_cmp(d.v, w2.v) > 0;
    }
    else
    {
        mpfr_abs(w1.v, d.v, MPFR_RNDN);
        mpfr_set_d(w2.v, 0.5, MPFR_RNDN);
        mpfr_add(w2.v, w2.v, tl.v, MPFR_RNDN);
        ok = mpfr_cmp(w1.v, w2.v) <= 0;
    }
    if (!ok)
        return say("misround", du);
    if (dir == 0)
    {
        mpfr_abs(w1.v, d.v, MPFR_RNDN);
        mpfr_sub_d(w1.v, w1.v, 0.5, MPFR_RNDN);
        mpfr_abs(w1.v, w1.v, MPFR_RNDN);
        if (mpfr_cmp(w1.v, tl.v) < 0)
        { // an exact tie
            if (mode == 0)
            {
                if (!sub && static_cast<int>(s % 2U) != 0)
                    return say("tie not even", du);
            }
            else
            {
                const bool away = neg ? du < 0 : du > 0;
                if (!away)
                    return say("tie not away", du);
            }
        }
    }
    return true;
}
template <class T> std::vector<T> small_grid(int lo, int hi, bool thin)
{
    std::vector<T> v;
    std::vector<long long> ms{1,  2,   3,   5,   7,   9,   11,  12,  15,   25,   50,    75,
                              99, 101, 125, 333, 500, 707, 866, 999, 1001, 9999, 99999, 999999};
    std::vector<long long> m2;
    for (size_t i = 0; i < ms.size(); i += thin ? 3 : 1)
        m2.push_back(ms[i]);
    for (auto m : m2)
    {
        int dig = 0;
        for (long long t = m; t; t /= 10)
            ++dig;
        for (int e = lo; e <= hi; ++e)
            v.push_back(T{m, e - dig + 1});
    }
    return v;
}
// F or E by descending AGM with a continuous phase, |k| < 1
static void agm_ellint(mpfr_t k, mpfr_t phi, mpfr_t oF, mpfr_t oE)
{
    big_float a, b, c, p, t, t2, pi, sum, sn, two;
    mpfr_const_pi(pi.v, MPFR_RNDN);
    mpfr_set_ui(a.v, 1, MPFR_RNDN);
    mpfr_sqr(t.v, k, MPFR_RNDN);
    mpfr_ui_sub(t.v, 1, t.v, MPFR_RNDN);
    mpfr_sqrt(b.v, t.v, MPFR_RNDN);
    mpfr_abs(c.v, k, MPFR_RNDN);
    mpfr_set(p.v, phi, MPFR_RNDN);
    mpfr_sqr(sum.v, c.v, MPFR_RNDN);
    mpfr_div_ui(sum.v, sum.v, 2, MPFR_RNDN);
    mpfr_set_ui(sn.v, 0, MPFR_RNDN);
    mpfr_set_ui(two.v, 1, MPFR_RNDN);
    for (int n = 1; n < 40; ++n)
    {
        if (mpfr_zero_p(c.v) || mpfr_get_exp(c.v) < -g_precision - 10)
            break;
        mpfr_div(t.v, b.v, a.v, MPFR_RNDN);
        mpfr_tan(t2.v, p.v, MPFR_RNDN);
        mpfr_mul(t.v, t.v, t2.v, MPFR_RNDN);
        mpfr_atan(t.v, t.v, MPFR_RNDN);
        mpfr_div(t2.v, p.v, pi.v, MPFR_RNDN);
        mpfr_rint(t2.v, t2.v, MPFR_RNDN);
        mpfr_mul(t2.v, t2.v, pi.v, MPFR_RNDN);
        mpfr_add(t.v, t.v, t2.v, MPFR_RNDN);
        mpfr_add(p.v, p.v, t.v, MPFR_RNDN);
        mpfr_add(t.v, a.v, b.v, MPFR_RNDN);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        mpfr_sub(c.v, a.v, b.v, MPFR_RNDN);
        mpfr_div_ui(c.v, c.v, 2, MPFR_RNDN);
        mpfr_mul(t2.v, a.v, b.v, MPFR_RNDN);
        mpfr_sqrt(b.v, t2.v, MPFR_RNDN);
        mpfr_set(a.v, t.v, MPFR_RNDN);
        mpfr_mul_ui(two.v, two.v, 2, MPFR_RNDN);
        mpfr_sqr(t.v, c.v, MPFR_RNDN);
        mpfr_mul(t.v, t.v, two.v, MPFR_RNDN);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        mpfr_add(sum.v, sum.v, t.v, MPFR_RNDN);
        mpfr_sin(t.v, p.v, MPFR_RNDN);
        mpfr_mul(t.v, t.v, c.v, MPFR_RNDN);
        mpfr_add(sn.v, sn.v, t.v, MPFR_RNDN);
    }
    mpfr_mul(t.v, two.v, a.v, MPFR_RNDN);
    mpfr_div(oF, p.v, t.v, MPFR_RNDN);
    mpfr_ui_sub(t.v, 1, sum.v, MPFR_RNDN);
    mpfr_mul(t.v, t.v, oF, MPFR_RNDN);
    mpfr_add(oE, t.v, sn.v, MPFR_RNDN);
}
static void ellint_values(mpfr_t k, mpfr_t phi, mpfr_t oF, mpfr_t oE)
{
    big_float ak, s, t, pi;
    mpfr_abs(ak.v, k, MPFR_RNDN);
    int c = mpfr_cmp_ui(ak.v, 1);
    if (c < 0)
    {
        agm_ellint(ak.v, phi, oF, oE);
        return;
    }
    mpfr_const_pi(pi.v, MPFR_RNDN);
    if (c == 0)
    {
        mpfr_div(t.v, phi, pi.v, MPFR_RNDN);
        mpfr_rint(t.v, t.v, MPFR_RNDN);
        big_float pr, sn;
        mpfr_mul(pr.v, t.v, pi.v, MPFR_RNDN);
        mpfr_sub(pr.v, phi, pr.v, MPFR_RNDN);
        mpfr_sin(sn.v, pr.v, MPFR_RNDN);
        mpfr_mul_ui(t.v, t.v, 2, MPFR_RNDN);
        mpfr_add(oE, t.v, sn.v, MPFR_RNDN);
        big_float hp;
        mpfr_div_ui(hp.v, pi.v, 2, MPFR_RNDN);
        mpfr_abs(s.v, phi, MPFR_RNDN);
        if (mpfr_cmp(s.v, hp.v) >= 0)
            mpfr_set_inf(oF, mpfr_sgn(phi));
        else
        {
            mpfr_sin(s.v, phi, MPFR_RNDN);
            mpfr_atanh(oF, s.v, MPFR_RNDN);
        }
        return;
    }
    mpfr_set_nan(oF);
    mpfr_set_nan(oE);
}
static void comp_ellint_values(mpfr_t k, mpfr_t oK, mpfr_t oE)
{
    big_float a, b, c, t, sum, two, pi, ak;
    mpfr_abs(ak.v, k, MPFR_RNDN);
    int r = mpfr_cmp_ui(ak.v, 1);
    if (r > 0)
    {
        mpfr_set_nan(oK);
        mpfr_set_nan(oE);
        return;
    }
    if (r == 0)
    {
        mpfr_set_inf(oK, 1);
        mpfr_set_ui(oE, 1, MPFR_RNDN);
        return;
    }
    mpfr_set_ui(a.v, 1, MPFR_RNDN);
    mpfr_sqr(t.v, ak.v, MPFR_RNDN);
    mpfr_ui_sub(t.v, 1, t.v, MPFR_RNDN);
    mpfr_sqrt(b.v, t.v, MPFR_RNDN);
    mpfr_set(c.v, ak.v, MPFR_RNDN);
    mpfr_sqr(sum.v, c.v, MPFR_RNDN);
    mpfr_div_ui(sum.v, sum.v, 2, MPFR_RNDN);
    mpfr_set_ui(two.v, 1, MPFR_RNDN);
    for (int n = 1; n < 40 && !mpfr_zero_p(c.v) && mpfr_get_exp(c.v) >= -g_precision - 10; ++n)
    {
        mpfr_add(t.v, a.v, b.v, MPFR_RNDN);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        mpfr_sub(c.v, a.v, b.v, MPFR_RNDN);
        mpfr_div_ui(c.v, c.v, 2, MPFR_RNDN);
        big_float t2;
        mpfr_mul(t2.v, a.v, b.v, MPFR_RNDN);
        mpfr_sqrt(b.v, t2.v, MPFR_RNDN);
        mpfr_set(a.v, t.v, MPFR_RNDN);
        mpfr_mul_ui(two.v, two.v, 2, MPFR_RNDN);
        mpfr_sqr(t.v, c.v, MPFR_RNDN);
        mpfr_mul(t.v, t.v, two.v, MPFR_RNDN);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        mpfr_add(sum.v, sum.v, t.v, MPFR_RNDN);
    }
    mpfr_const_pi(pi.v, MPFR_RNDN);
    mpfr_mul_ui(t.v, a.v, 2, MPFR_RNDN);
    mpfr_div(oK, pi.v, t.v, MPFR_RNDN);
    mpfr_ui_sub(t.v, 1, sum.v, MPFR_RNDN);
    mpfr_mul(oE, t.v, oK, MPFR_RNDN);
}
static void comp_ellint_e(mpfr_t k, mpfr_t o)
{
    big_float x;
    comp_ellint_values(k, x.v, o);
}
template <class T> void check_ellint(const char *tn)
{
    auto K = small_grid<T>(-3, 0, false);
    auto P = small_grid<T>(-3, 1, true);
    for (int j = 1; j < std::numeric_limits<T>::digits10; j += 2)
    {
        K.push_back(T{1} - T{1, -j});
        K.push_back(T{1} + T{1, -j});
    }
    K.push_back(T{1});
    P.push_back(T{15707963267948966LL, -16});
    for (T c : {T{1}, T{2}, T{3}, T{10}, T{100}, T{1, 4}, T{1, 10}})
        P.push_back(c);
    long bF = 0, bE = 0, n = 0;
    big_float a, b, oF, oE;
    for (auto k0 : K)
        for (auto p0 : P)
            for (int sg = 0; sg < 4; ++sg)
            {
                T k = (sg & 1) ? -k0 : k0, p = (sg & 2) ? -p0 : p0;
                to_mpfr(k, a.v);
                to_mpfr(p, b.v);
                mpfr_clear_flags();
                ellint_values(a.v, b.v, oF.v, oE.v);
                g_overflow = g_underflow = false;
                std::string as = show(k) + "," + show(p);
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T f = ellint_1(k, p);
                    margin_events(oF.v, "ellint_1", as);
                    T e = ellint_2(k, p);
                    margin_events(oE.v, "ellint_2", as);
                    fesetround(g_modes[0]);
                    ++n;
                    check(f, oF.v, m, "ellint_1", as, bF);
                    check(e, oE.v, m, "ellint_2", as, bE);
                }
            }
    printf("%s ellint n=%ld badF=%ld badE=%ld\n", tn, n, bF, bE);
    fflush(stdout);
    g_bad_total += bF + bE;
    long bC = 0;
    n = 0;
    for (auto k0 : K)
        for (int sg = 0; sg < 2; ++sg)
        {
            T k = sg ? -k0 : k0;
            to_mpfr(k, a.v);
            comp_ellint_e(a.v, oE.v);
            g_overflow = g_underflow = false;
            std::string as = show(k);
            for (int m = 0; m < 5; ++m)
            {
                fesetround(g_modes[m]);
                T e = comp_ellint_2(k);
                fesetround(g_modes[0]);
                ++n;
                check(e, oE.v, m, "comp_ellint_2", as, bC);
            }
        }
    printf("%s comp_ellint_2 n=%ld bad=%ld\n", tn, n, bC);
    fflush(stdout);
    g_bad_total += bC;
}
static void poly_value(int kind, unsigned n, unsigned mm, mpfr_t x, mpfr_t o)
{
    big_float p0, p1, p2, t1, t2;
    if (kind == 1)
    {
        mpfr_abs(t1.v, x, MPFR_RNDN);
        if (mpfr_cmp_ui(t1.v, 1) > 0)
        {
            mpfr_set_nan(o);
            return;
        }
    }
    if (kind == 4)
    {
        big_float c;
        mpfr_abs(t1.v, x, MPFR_RNDN);
        if (mpfr_cmp_ui(t1.v, 1) > 0)
        {
            mpfr_set_nan(o);
            return;
        }
        if (mm > n)
        {
            mpfr_set_ui(o, 0, MPFR_RNDN);
            return;
        }
        mpfr_sqr(c.v, x, MPFR_RNDN);
        mpfr_ui_sub(c.v, 1, c.v, MPFR_RNDN);
        mpfr_sqrt(c.v, c.v, MPFR_RNDN);
        mpfr_set_ui(p0.v, 1, MPFR_RNDN);
        for (unsigned i = 0; i < mm; ++i)
        {
            mpfr_mul_ui(p0.v, p0.v, 2 * i + 1, MPFR_RNDN);
            mpfr_mul(p0.v, p0.v, c.v, MPFR_RNDN);
        }
        if (n == mm)
        {
            mpfr_set(o, p0.v, MPFR_RNDN);
            return;
        }
        mpfr_mul(p1.v, x, p0.v, MPFR_RNDN);
        mpfr_mul_ui(p1.v, p1.v, 2 * mm + 1, MPFR_RNDN);
        for (unsigned l = mm + 1; l < n; ++l)
        {
            mpfr_mul(t2.v, x, p1.v, MPFR_RNDN);
            mpfr_mul_ui(t2.v, t2.v, 2 * l + 1, MPFR_RNDN);
            mpfr_mul_ui(t1.v, p0.v, l + mm, MPFR_RNDN);
            mpfr_sub(p2.v, t2.v, t1.v, MPFR_RNDN);
            mpfr_div_ui(p2.v, p2.v, l - mm + 1, MPFR_RNDN);
            mpfr_swap(p0.v, p1.v);
            mpfr_swap(p1.v, p2.v);
        }
        // A root of the polynomial has no sign, where x < 0 makes this -0
        mpfr_set(o, p1.v, MPFR_RNDN);
        if (mpfr_zero_p(o))
            mpfr_set_ui(o, 0, MPFR_RNDN);
        return;
    }
    mpfr_set_ui(p0.v, 1, MPFR_RNDN);
    if (n == 0)
    {
        mpfr_set_ui(o, 1, MPFR_RNDN);
        return;
    }
    if (kind == 0)
        mpfr_mul_ui(p1.v, x, 2, MPFR_RNDN);
    else if (kind == 1)
        mpfr_set(p1.v, x, MPFR_RNDN);
    else
        mpfr_ui_sub(p1.v, 1 + mm, x, MPFR_RNDN);
    for (unsigned k = 1; k < n; ++k)
    {
        if (kind == 0)
        {
            mpfr_mul(t1.v, x, p1.v, MPFR_RNDN);
            mpfr_mul_ui(t1.v, t1.v, 2, MPFR_RNDN);
            mpfr_mul_ui(t2.v, p0.v, 2 * k, MPFR_RNDN);
            mpfr_sub(p2.v, t1.v, t2.v, MPFR_RNDN);
        }
        else if (kind == 1)
        {
            mpfr_mul(t1.v, x, p1.v, MPFR_RNDN);
            mpfr_mul_ui(t1.v, t1.v, 2 * k + 1, MPFR_RNDN);
            mpfr_mul_ui(t2.v, p0.v, k, MPFR_RNDN);
            mpfr_sub(p2.v, t1.v, t2.v, MPFR_RNDN);
            mpfr_div_ui(p2.v, p2.v, k + 1, MPFR_RNDN);
        }
        else
        {
            mpfr_ui_sub(t1.v, 2 * k + 1 + mm, x, MPFR_RNDN);
            mpfr_mul(t1.v, t1.v, p1.v, MPFR_RNDN);
            mpfr_mul_ui(t2.v, p0.v, k + mm, MPFR_RNDN);
            mpfr_sub(p2.v, t1.v, t2.v, MPFR_RNDN);
            mpfr_div_ui(p2.v, p2.v, k + 1, MPFR_RNDN);
        }
        mpfr_swap(p0.v, p1.v);
        mpfr_swap(p1.v, p2.v);
    }
    mpfr_set(o, p1.v, MPFR_RNDN);
}
template <class T> void check_assoc(const char *tn)
{
    auto X = small_grid<T>(-5, 3, false);
    auto Xs = small_grid<T>(-3, 0, false);
    Xs.push_back(T{1});
    long bL = 0, bP = 0, n = 0;
    big_float a, res;
    g_overflow = g_underflow = false;
    for (unsigned nn : {0u, 1u, 2u, 3u, 5u, 8u, 16u, 40u, 100u})
        for (unsigned mm : {0u, 1u, 2u, 3u, 5u, 20u, 100u})
            for (auto x0 : X)
                for (int sg = 0; sg < 2; ++sg)
                {
                    T x = sg ? -x0 : x0;
                    to_mpfr(x, a.v);
                    poly_value(3, nn, mm, a.v, res.v);
                    std::string as = std::to_string(nn) + "," + std::to_string(mm) + "," + show(x);
                    for (int m = 0; m < 5; ++m)
                    {
                        fesetround(g_modes[m]);
                        T r = assoc_laguerre(nn, mm, x);
                        fesetround(g_modes[0]);
                        ++n;
                        check(r, res.v, m, "assoc_laguerre", as, bL);
                    }
                }
    printf("%s assoc_laguerre n=%ld bad=%ld\n", tn, n, bL);
    fflush(stdout);
    n = 0;
    for (unsigned nn : {0u, 1u, 2u, 3u, 5u, 8u, 16u, 40u, 100u})
        for (unsigned mm : {0u, 1u, 2u, 3u, 4u, 5u, 8u, 20u, 100u})
            for (auto x0 : Xs)
                for (int sg = 0; sg < 2; ++sg)
                {
                    T x = sg ? -x0 : x0;
                    to_mpfr(x, a.v);
                    poly_value(4, nn, mm, a.v, res.v);
                    std::string as = std::to_string(nn) + "," + std::to_string(mm) + "," + show(x);
                    for (int m = 0; m < 5; ++m)
                    {
                        fesetround(g_modes[m]);
                        T r = assoc_legendre(nn, mm, x);
                        fesetround(g_modes[0]);
                        ++n;
                        check(r, res.v, m, "assoc_legendre", as, bP);
                    }
                }
    printf("%s assoc_legendre n=%ld bad=%ld\n", tn, n, bP);
    fflush(stdout);
    g_bad_total += bL + bP;
}
// Tiny arguments of these functions give results nearer to a rounding boundary than the strict window of check, so
// they take a fixed tolerance. The tolerance also applies to all checks of check_extremes.
template <class T, class FD, class FM>
void check_unary(const char *name, const std::vector<T> &G, FD fd, FM fm, int prec, double lim = 1e308)
{
    g_precision = prec;
    std::string nm = name;
    g_fixed_tolerance = nm == "log1p_tiny" || nm == "expm1" || nm == "tanh" || nm == "erf" || nm == "erfc" ||
                        nm == "zeta" || nm == "atan" || nm == "exp2" || nm == "exp" || nm == "cosh";
    long bad = 0, n = 0;
    big_float a, t;
    for (auto x0 : G)
        for (int sg = 0; sg < 2; ++sg)
        {
            T x = sg ? -x0 : x0;
            to_mpfr(x, a.v);
            if (mpfr_cmp_d(a.v, lim) > 0 || mpfr_cmp_d(a.v, -lim) < 0)
                continue;
            mpfr_clear_flags();
            fm(t.v, a.v);
            g_overflow = mpfr_overflow_p();
            g_underflow = mpfr_underflow_p();
            std::string as = show(x);
            for (int m = 0; m < 5; ++m)
            {
                fesetround(g_modes[m]);
                T r = fd(x);
                fesetround(g_modes[0]);
                ++n;
                check(r, t.v, m, name, as, bad);
            }
        }
    g_fixed_tolerance = false;
    printf("%s n=%ld bad=%ld\n", name, n, bad);
    fflush(stdout);
    g_bad_total += bad;
}
// inputs near min, denorm_min and max, and the bands that give results there
template <class T> void check_extremes(const char *tn)
{
    int edn;
    {
        auto q = frexp10(std::numeric_limits<T>::denorm_min(), &edn);
        int dg = 0;
        for (auto w = q; w != 0; w /= 10)
            ++dg;
        edn += dg - 1;
    }
    g_strict = false;
    int emax = std::numeric_limits<T>::max_exponent10;
    std::vector<T> X;
    const int Pd = std::numeric_limits<T>::digits10;
    int hiprec = Pd > 30 ? 45000 : (Pd > 10 ? 3000 : 1000);
    for (long long m : {1LL, 2LL, 3LL, 7LL, 99LL, 101LL, 999LL, 1001LL, 123457LL})
    {
        int dig = 0;
        for (long long t = m; t; t /= 10)
            ++dig;
        for (int e : {edn, edn + 1, edn + 2, edn + 5, edn + Pd - 1, edn + Pd, edn + Pd + 1, edn + 2 * Pd, -emax / 2,
                      emax - Pd, emax - 5, emax - 2, emax - 1})
        {
            T v{m, e - dig + 1};
            if (!isinf(v) && v != T{0})
                X.push_back(v);
        }
    }
    int lim = emax;
    for (int e : {lim / 10, lim / 5, lim / 3, lim / 2, lim, lim * 2, lim * 3})
        for (long long m : {1LL, 2LL, 5LL, 23LL, 1000LL})
            X.push_back(T{m * e, 0});
    printf("== ext %s n=%zu\n", tn, X.size());
    // check_unary takes each value with both signs, so this also checks the sign of the result for -0
    X.push_back(T{0});
    check_unary<T>("exp", X, DEC1(exp), MPFR1(mpfr_exp), hiprec);
    check_unary<T>("expm1", X, DEC1(expm1), MPFR1(mpfr_expm1), hiprec);
    check_unary<T>("log1p", X, DEC1(log1p), MPFR1(mpfr_log1p), hiprec);
    check_unary<T>("log", X, DEC1(log), MPFR1(mpfr_log), hiprec);
    check_unary<T>("log2", X, DEC1(log2), MPFR1(mpfr_log2), hiprec);
    check_unary<T>("log10", X, DEC1(log10), MPFR1(mpfr_log10), hiprec);
    check_unary<T>("sin", X, DEC1(sin), MPFR1(mpfr_sin), hiprec);
    check_unary<T>("cos", X, DEC1(cos), MPFR1(mpfr_cos), hiprec);
    check_unary<T>("tan", X, DEC1(tan), MPFR1(mpfr_tan), hiprec);
    check_unary<T>("asin", X, DEC1(asin), MPFR1(mpfr_asin), hiprec);
    check_unary<T>("acos", X, DEC1(acos), MPFR1(mpfr_acos), hiprec);
    check_unary<T>("atan", X, DEC1(atan), MPFR1(mpfr_atan), hiprec);
    check_unary<T>("sinh", X, DEC1(sinh), MPFR1(mpfr_sinh), hiprec);
    check_unary<T>("cosh", X, DEC1(cosh), MPFR1(mpfr_cosh), hiprec);
    check_unary<T>("tanh", X, DEC1(tanh), MPFR1(mpfr_tanh), hiprec);
    check_unary<T>("asinh", X, DEC1(asinh), MPFR1(mpfr_asinh), hiprec);
    check_unary<T>("acosh", X, DEC1(acosh), MPFR1(mpfr_acosh), hiprec);
    check_unary<T>("atanh", X, DEC1(atanh), MPFR1(mpfr_atanh), hiprec);
    check_unary<T>("cbrt", X, DEC1(cbrt), MPFR1(mpfr_cbrt), hiprec);
    check_unary<T>("sqrt", X, DEC1(sqrt), MPFR1(mpfr_sqrt), hiprec);
    check_unary<T>("erf", X, DEC1(erf), MPFR1(mpfr_erf), 3000);
    check_unary<T>("erfc", X, DEC1(erfc), MPFR1(mpfr_erfc), 3000);
    check_unary<T>("tgamma", X, DEC1(tgamma), MPFR1(mpfr_gamma), 3000, 1e6);
    check_unary<T>(
        "lgamma", X, DEC1(lgamma),
        [](mpfr_t o, mpfr_t a)
        {
            int s;
            mpfr_lgamma(o, &s, a, MPFR_RNDN);
        },
        3000, 1e6);
    check_unary<T>("zeta", X, DEC1(riemann_zeta), MPFR1(mpfr_zeta), 3000, 1e4);
    check_unary<T>(
        "comp_ellint_1", X, DEC1(comp_ellint_1),
        [](mpfr_t o, mpfr_t k)
        {
            big_float x;
            comp_ellint_values(k, o, x.v);
        },
        3000);
    check_unary<T>("comp_ellint_2", X, DEC1(comp_ellint_2), [](mpfr_t o, mpfr_t k) { comp_ellint_e(k, o); }, 3000);
    g_precision = 500;
    g_strict = true;
}

template <class T> std::vector<T> wide_grid(bool small)
{
    std::vector<T> v;
    std::vector<long long> ms{1,    2,     3,     4,     5,      6,      7,       8,      9,    11,   12,
                              15,   16,    25,    31,    32,     50,     64,      75,     99,   101,  125,
                              256,  333,   500,   512,   625,    999,    1001,    1024,   3142, 5000, 7071,
                              9999, 10001, 31416, 99999, 100001, 999999, 1000001, 9999999};
    const int P = std::numeric_limits<T>::digits10;
    int PP = P > 18 ? 18 : P;
    long long big = 1;
    for (int i = 1; i < PP; ++i)
        big *= 10;
    for (long long d : {1LL, 2LL, 3LL})
    {
        ms.push_back(big + d);
        ms.push_back(big * 10 - d);
        ms.push_back(big * 5 + d);
        ms.push_back(big * 5 - d);
    }
    if (small)
    {
        std::vector<long long> m2;
        for (size_t i = 0; i < ms.size(); i += 3)
            m2.push_back(ms[i]);
        ms = m2;
    }
    int lo = small ? -5 : -20, hi = small ? 3 : 6;
    for (auto m : ms)
    {
        int dig = 0;
        for (long long t = m; t; t /= 10)
            ++dig;
        for (int e = lo; e <= hi; ++e)
            v.push_back(T{m, e - dig + 1});
    }
    return v;
}
template <class T, class FD, class FM>
void check_binary(const char *name, const std::vector<T> &H, const std::vector<T> &K, FD fd, FM fm, int sgs = 4)
{
    long bad = 0, n = 0;
    big_float a, b, t;
    for (auto x0 : H)
        for (auto y0 : K)
            for (int sg = 0; sg < sgs; ++sg)
            {
                T x = (sg & 1) ? -x0 : x0, y = (sg & 2) ? -y0 : y0;
                to_mpfr(x, a.v);
                to_mpfr(y, b.v);
                mpfr_clear_flags();
                fm(t.v, a.v, b.v);
                g_overflow = mpfr_overflow_p();
                g_underflow = mpfr_underflow_p();
                std::string as = show(x) + "," + show(y);
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T r = fd(x, y);
                    fesetround(g_modes[0]);
                    ++n;
                    check(r, t.v, m, name, as, bad);
                }
            }
    printf("%s n=%ld bad=%ld\n", name, n, bad);
    fflush(stdout);
    g_bad_total += bad;
}
template <class T, class FD> void check_polynomial(const char *name, const std::vector<T> &H, FD fd, int kind)
{
    g_overflow = g_underflow = false;
    long bad = 0, n = 0;
    big_float a, res;
    for (unsigned l : {2u, 3u, 4u, 5u, 8u, 16u, 40u, 100u})
        for (auto x0 : H)
            for (int sg = 0; sg < 2; ++sg)
            {
                T x = sg ? -x0 : x0;
                to_mpfr(x, a.v);
                poly_value(kind, l, 0, a.v, res.v);
                std::string as = std::to_string(l) + "," + show(x);
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T r = fd(l, x);
                    fesetround(g_modes[0]);
                    ++n;
                    check(r, res.v, m, name, as, bad);
                }
            }
    printf("%s n=%ld bad=%ld\n", name, n, bad);
    fflush(stdout);
    g_bad_total += bad;
}
template <class T> void check_main(const char *tn)
{
    g_precision = 700;
    auto G = wide_grid<T>(false);
    auto H = wide_grid<T>(true);
    printf("== %s G=%zu H=%zu\n", tn, G.size(), H.size());
    std::vector<T> H2;
    for (size_t i = 0; i < H.size(); i += 4)
        H2.push_back(H[i]);
    for (int m : {1, 2, 3, 5, 7})
        for (int e = -3; e <= 2; ++e)
            H2.push_back(T{m, e});
    check_unary<T>("exp", G, DEC1(exp), MPFR1(mpfr_exp), 700, 1e5);
    check_unary<T>("exp2", G, DEC1(exp2), MPFR1(mpfr_exp2), 700, 1e5);
    check_unary<T>("expm1", G, DEC1(expm1), MPFR1(mpfr_expm1), 700, 1e5);
    check_unary<T>("log", G, DEC1(log), MPFR1(mpfr_log), 700);
    check_unary<T>("log2", G, DEC1(log2), MPFR1(mpfr_log2), 700);
    check_unary<T>("log10", G, DEC1(log10), MPFR1(mpfr_log10), 700);
    check_unary<T>("log1p", G, DEC1(log1p), MPFR1(mpfr_log1p), 700);
    check_unary<T>("sin", G, DEC1(sin), MPFR1(mpfr_sin), 700);
    check_unary<T>("cos", G, DEC1(cos), MPFR1(mpfr_cos), 700);
    check_unary<T>("tan", G, DEC1(tan), MPFR1(mpfr_tan), 700);
    check_unary<T>("asin", G, DEC1(asin), MPFR1(mpfr_asin), 700);
    check_unary<T>("acos", G, DEC1(acos), MPFR1(mpfr_acos), 700);
    check_unary<T>("atan", G, DEC1(atan), MPFR1(mpfr_atan), 700);
    check_unary<T>("sinh", G, DEC1(sinh), MPFR1(mpfr_sinh), 700, 1e5);
    check_unary<T>("cosh", G, DEC1(cosh), MPFR1(mpfr_cosh), 700, 1e5);
    check_unary<T>("tanh", G, DEC1(tanh), MPFR1(mpfr_tanh), 700, 1e5);
    check_unary<T>("asinh", G, DEC1(asinh), MPFR1(mpfr_asinh), 700);
    check_unary<T>("acosh", G, DEC1(acosh), MPFR1(mpfr_acosh), 700);
    check_unary<T>("atanh", G, DEC1(atanh), MPFR1(mpfr_atanh), 700);
    check_unary<T>("cbrt", G, DEC1(cbrt), MPFR1(mpfr_cbrt), 700);
    check_unary<T>("sqrt", G, DEC1(sqrt), MPFR1(mpfr_sqrt), 700);
    check_unary<T>("erf", G, DEC1(erf), MPFR1(mpfr_erf), 700, 1e4);
    check_unary<T>("erfc", G, DEC1(erfc), MPFR1(mpfr_erfc), 700, 1e4);
    check_unary<T>("tgamma", G, DEC1(tgamma), MPFR1(mpfr_gamma), 700, 800);
    check_unary<T>(
        "lgamma", G, DEC1(lgamma),
        [](mpfr_t o, mpfr_t a)
        {
            int s;
            mpfr_lgamma(o, &s, a, MPFR_RNDN);
        },
        700, 1e6);
    check_unary<T>("zeta", G, DEC1(riemann_zeta), MPFR1(mpfr_zeta), 700, 300);
    check_binary<T>(
        "atan2", H2, H2, [](T a, T b) { return atan2(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_atan2(o, a, b, MPFR_RNDN); });
    check_binary<T>(
        "pow", H2, H2, [](T a, T b) { return pow(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_pow(o, a, b, MPFR_RNDN); });
    check_binary<T>(
        "hypot", H2, H2, [](T a, T b) { return hypot(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_hypot(o, a, b, MPFR_RNDN); });
    check_binary<T>(
        "beta", H2, H2, [](T a, T b) { return beta(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_beta(o, a, b, MPFR_RNDN); }, 1);
    g_precision = 700;
    check_polynomial<T>("hermite", H, [](unsigned l, T x) { return hermite(l, x); }, 0);
    check_polynomial<T>("legendre", H, [](unsigned l, T x) { return legendre(l, x); }, 1);
    check_polynomial<T>("laguerre", H, [](unsigned l, T x) { return laguerre(l, x); }, 2);
}
// exact results: ties and exact values in every mode
template <class T> void check_exact()
{
    g_precision = 1500;
    std::vector<T> B;
    auto addm = [&](long long m)
    {
        for (int e : {-6, -3, -1, 0, 1, 3, 6})
            B.push_back(T{m, e});
    };
    for (long long m :
         {2LL,      3LL,       4LL,       5LL,       6LL,        7LL,         8LL,          9LL,          11LL,
          12LL,     13LL,      15LL,      16LL,      25LL,       32LL,        35LL,         64LL,         75LL,
          99LL,     101LL,     125LL,     128LL,     175LL,      225LL,       256LL,        625LL,        1001LL,
          1024LL,   1025LL,    3125LL,    4096LL,    10001LL,    15625LL,     12345LL,      65536LL,      78125LL,
          123457LL, 1048576LL, 1953125LL, 3162278LL, 31622777LL, 123456789LL, 3162277660LL, 31622776602LL})
        addm(m);
    for (long long m :
         {3162277660168379LL, 99999999LL, 9999999999LL, 100000001LL, 10000000001LL, 99999995LL, 316227766016838LL})
        if (m < 1000000000000000000LL)
        {
            B.push_back(T{m, 0});
            B.push_back(T{m, -8});
        }
    std::vector<T> Y;
    for (int i = -12; i <= 16; ++i)
        Y.push_back(T{i});
    for (int m : {5, 15, 25, 35, -5, -15, -25})
        Y.push_back(T{m, -1});
    for (int m : {25, 75, 125, -25, -75})
        Y.push_back(T{m, -2});
    Y.push_back(T{100});
    Y.push_back(T{1000});
    check_binary<T>(
        "pow_x", B, Y, [](T a, T b) { return pow(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_pow(o, a, b, MPFR_RNDN); }, 2);
    {
        long bad = 0, n = 0;
        big_float a, t;
        for (auto x0 : B)
            for (int sg = 0; sg < 2; ++sg)
                for (int k : {-1000, -100, -16, -12, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0,  1,   2,   3,
                              4,     5,    6,   7,   8,  9,  10, 11, 12, 13, 16, 17, 25, 50, 100, 1000})
                {
                    T x = sg ? -x0 : x0;
                    to_mpfr(x, a.v);
                    mpfr_clear_flags();
                    mpfr_pow_si(t.v, a.v, k, MPFR_RNDN);
                    g_overflow = mpfr_overflow_p();
                    g_underflow = mpfr_underflow_p();
                    std::string as = show(x) + ",int " + std::to_string(k);
                    for (int m = 0; m < 5; ++m)
                    {
                        fesetround(g_modes[m]);
                        T r = pow(x, k);
                        fesetround(g_modes[0]);
                        ++n;
                        check(r, t.v, m, "pow_int", as, bad);
                    }
                }
        printf("pow_int n=%ld bad=%ld\n", n, bad);
        g_bad_total += bad;
    }
    check_binary<T>(
        "hypot_x", B, B, [](T a, T b) { return hypot(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_hypot(o, a, b, MPFR_RNDN); }, 1);
    {
        long bad = 0, n = 0;
        big_float a, b, c, t, t2;
        g_overflow = g_underflow = false;
        std::vector<T> J;
        for (long long m : {1LL, 3LL, 5LL, 7LL, 11LL, 1001LL, 3333LL, 12345LL, 333333LL, 3333333LL, 33333333LL,
                            3333333333LL, 5555555555LL})
            for (int e : {-3, 0})
                J.push_back(T{m, e});
        for (auto j : J)
            for (auto tr : {std::array<int, 3>{1, 2, 2},
                            {2, 3, 6},
                            {1, 4, 8},
                            {4, 4, 7},
                            {2, 6, 9},
                            {6, 6, 7},
                            {0, 3, 4},
                            {5, 0, 12}})
            {
                T x = j * T{tr[0]}, y = j * T{tr[1]}, z = j * T{tr[2]};
                to_mpfr(x, a.v);
                to_mpfr(y, b.v);
                to_mpfr(z, c.v);
                mpfr_sqr(t.v, a.v, MPFR_RNDN);
                mpfr_sqr(t2.v, b.v, MPFR_RNDN);
                mpfr_add(t.v, t.v, t2.v, MPFR_RNDN);
                mpfr_sqr(t2.v, c.v, MPFR_RNDN);
                mpfr_add(t.v, t.v, t2.v, MPFR_RNDN);
                mpfr_sqrt(t.v, t.v, MPFR_RNDN);
                std::string as = show(x) + "," + show(y) + "," + show(z);
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T r = hypot(x, y, z);
                    fesetround(g_modes[0]);
                    ++n;
                    check(r, t.v, m, "hypot3_x", as, bad);
                }
            }
        printf("hypot3_x n=%ld bad=%ld\n", n, bad);
        g_bad_total += bad;
    }
    std::vector<T> I;
    for (int i = -400; i <= 400; ++i)
        I.push_back(T{i});
    for (int i = -20; i <= 20; ++i)
        I.push_back(T{2 * i + 1, 0} / T{2});
    check_unary<T>("exp2_x", I, DEC1(exp2), MPFR1(mpfr_exp2), 1500);
    std::vector<T> Pw2, Pw10;
    {
        T a{1};
        for (int i = 0; i < 200 && !isinf(a); ++i)
        {
            Pw2.push_back(a);
            a = a * T{2};
        }
    }
    for (int e = -90; e <= 90; ++e)
        Pw10.push_back(T{1, e});
    check_unary<T>("log2_x", Pw2, DEC1(log2), MPFR1(mpfr_log2), 1500);
    check_unary<T>("log10_x", Pw10, DEC1(log10), MPFR1(mpfr_log10), 1500);
    std::vector<T> F;
    for (int i = 1; i <= 170; ++i)
        F.push_back(T{i});
    for (int i : {171, 200, 300, 500, 1000, 1500, 1754, 2000})
        F.push_back(T{i});
    for (int i = 1; i <= 30; ++i)
        F.push_back(T{2 * i - 1} / T{2});
    check_unary<T>("tgamma_x", F, DEC1(tgamma), MPFR1(mpfr_gamma), 1500, 5000);
    {
        std::vector<T> A;
        for (int i : {1, 2, 3, 4, 5, 6, 8, 10, 16, 20, 25, 32, 40, 50, 64, 100, 125, 200, 1000})
            A.push_back(T{i});
        check_binary<T>(
            "beta_x", A, A, [](T a, T b) { return beta(a, b); },
            [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_beta(o, a, b, MPFR_RNDN); }, 1);
    }
    check_binary<T>(
        "atan2_x", B, Y, [](T a, T b) { return atan2(a, b); },
        [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_atan2(o, a, b, MPFR_RNDN); }, 4);
}

template <class T> T from_mpfr(mpfr_t v)
{
    const int P = std::numeric_limits<T>::digits10;
    mpfr_exp_t ex;
    char *st = mpfr_get_str(nullptr, &ex, 10, P, v, MPFR_RNDN);
    std::string o(st);
    mpfr_free_str(st);
    bool n = false;
    if (o[0] == '-')
    {
        n = true;
        o = o.substr(1);
    }
    std::string z = o.substr(0, 1) + "." + o.substr(1) + "e" + std::to_string(static_cast<long>(ex) - 1);
    T r{};
    boost::decimal::from_chars(z.data(), z.data() + z.size(), r);
    return n ? -r : r;
}
template <class T> T ulp_of(T x)
{
    int e;
    frexp10(x, &e);
    return T{1, e};
}
template <class T> void add_around(std::vector<T> &v, T b)
{
    if (isnan(b) || isinf(b))
        return;
    const int P = std::numeric_limits<T>::digits10;
    v.push_back(b);
    T u = ulp_of(b == T{0} ? T{1} : b);
    if (b == T{0})
        u = T{1, -P};
    std::vector<long long> ks{1, 2, 3, 5, 9, 10, 11, 50, 99, 100};
    long long p10 = 1000;
    for (int i = 3; i < P - 1 && i < 18; ++i)
    {
        ks.push_back(p10);
        p10 *= 10;
    }
    for (auto k : ks)
    {
        T d = u * T{k};
        T a = b + d, c = b - d;
        if (!isinf(a))
            v.push_back(a);
        if (!isinf(c) && c != T{0})
            v.push_back(c);
    }
}
template <class T> std::vector<T> edges()
{
    std::vector<T> v;
    const char *E[] = {"0.00390625", "0.015625", "0.0625", "0.1", "0.11", "0.125", "0.2",   "0.25",  "0.3",
                       "0.375",      "0.5",      "0.625",  "0.7", "0.71", "0.75",  "0.866", "0.875", "0.9",
                       "1",          "1.25",     "1.5",    "1.6", "1.75", "2",     "2.5",   "3",     "4",
                       "5",          "6",        "7",      "8",   "9",    "10",    "12",    "16",    "20",
                       "32",         "64",       "100",    "256", "1000", "10000", "100000"};
    for (auto e : E)
    {
        std::string st(e);
        T x{};
        boost::decimal::from_chars(st.data(), st.data() + st.size(), x);
        add_around(v, x);
    }
    return v;
}
// x in [lo, hi] with f(x) = target, by bisection
template <class F> bool bisect(F f, double lo, double hi, mpfr_t out)
{
    big_float a, b, m, fa, fm;
    mpfr_set_d(a.v, lo, MPFR_RNDN);
    mpfr_set_d(b.v, hi, MPFR_RNDN);
    f(fa.v, a.v);
    f(fm.v, b.v);
    if (mpfr_nan_p(fa.v) || mpfr_nan_p(fm.v) || mpfr_sgn(fa.v) * mpfr_sgn(fm.v) > 0)
        return false;
    int sa = mpfr_sgn(fa.v);
    for (int i = 0; i < 260; ++i)
    {
        mpfr_add(m.v, a.v, b.v, MPFR_RNDN);
        mpfr_div_ui(m.v, m.v, 2, MPFR_RNDN);
        f(fm.v, m.v);
        if (mpfr_sgn(fm.v) == 0)
            break;
        if (mpfr_sgn(fm.v) == sa)
            mpfr_set(a.v, m.v, MPFR_RNDN);
        else
            mpfr_set(b.v, m.v, MPFR_RNDN);
    }
    mpfr_set(out, m.v, MPFR_RNDN);
    return true;
}
template <class T, class FD>
void roots(const char *name, int kind, unsigned n, unsigned mm, double lo, double hi, FD fd, long &bad, long &cnt)
{
    big_float x, y, r;
    std::vector<T> v;
    const int N = 400;
    for (int i = 0; i < N; ++i)
    {
        double a = lo + (hi - lo) * i / N, b = lo + (hi - lo) * (i + 1) / N;
        if (bisect([&](mpfr_t o, mpfr_t q) { poly_value(kind, n, mm, q, o); }, a, b, r.v))
            add_around(v, from_mpfr<T>(r.v));
    }
    g_overflow = g_underflow = false;
    for (T q : v)
    {
        to_mpfr(q, x.v);
        poly_value(kind, n, mm, x.v, y.v);
        std::string as = std::to_string(n) + "," + std::to_string(mm) + "," + show(q);
        for (int m = 0; m < 5; ++m)
        {
            fesetround(g_modes[m]);
            T res = fd(q);
            fesetround(g_modes[0]);
            ++cnt;
            check(res, y.v, m, name, as, bad);
        }
    }
}
// beta(x, y) with one argument at or above 10^15, from the Stirling series of ln(gamma(M + m)) - ln(gamma(M))
static void large_beta_value(mpfr_t o, mpfr_t x, mpfr_t y)
{
    static const long num[] = {1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611, 854513, -236364091};
    static const long den[] = {6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138, 2730};
    const bool x_large = mpfr_cmp(x, y) > 0;
    const int prec = g_precision + 2 * static_cast<int>(mpfr_get_exp(x_large ? x : y));
    mpfr_t m, big, sum, a, b, t;
    mpfr_inits2(prec, m, big, sum, a, b, t, static_cast<mpfr_ptr>(nullptr));
    mpfr_set(m, x_large ? y : x, MPFR_RNDN);
    mpfr_set(big, x_large ? x : y, MPFR_RNDN);
    mpfr_add(a, big, m, MPFR_RNDN);
    if (mpfr_integer_p(m) && mpfr_cmp_ui(m, 1000) <= 0)
    {
        // beta(m, M) = (m - 1)! / (M (M + 1) ... (M + m - 1)) for an integer m
        mpfr_set_ui(t, 1, MPFR_RNDN);
        mpfr_set(b, big, MPFR_RNDN);
        for (unsigned long j = 1; mpfr_cmp_ui(m, j) > 0; ++j)
        {
            mpfr_mul_ui(t, t, j, MPFR_RNDN);
            mpfr_add_ui(b, b, 1, MPFR_RNDN);
        }
        mpfr_set(sum, big, MPFR_RNDN);
        for (mpfr_set(a, big, MPFR_RNDN); mpfr_cmp(a, b) < 0;)
        {
            mpfr_add_ui(a, a, 1, MPFR_RNDN);
            mpfr_mul(sum, sum, a, MPFR_RNDN);
        }
        mpfr_div(t, t, sum, MPFR_RNDN);
        mpfr_set_prec(o, prec);
        mpfr_set(o, t, MPFR_RNDN);
        mpfr_clears(m, big, sum, a, b, t, static_cast<mpfr_ptr>(nullptr));
        return;
    }
    // (M + m - 1/2) ln(M + m) - (M - 1/2) ln(M) - m
    mpfr_sub_d(t, a, 0.5, MPFR_RNDN);
    mpfr_log(sum, a, MPFR_RNDN);
    mpfr_mul(sum, sum, t, MPFR_RNDN);
    mpfr_sub_d(t, big, 0.5, MPFR_RNDN);
    mpfr_log(b, big, MPFR_RNDN);
    mpfr_mul(t, t, b, MPFR_RNDN);
    mpfr_sub(sum, sum, t, MPFR_RNDN);
    mpfr_sub(sum, sum, m, MPFR_RNDN);
    for (int j = 1; j <= 12; ++j)
    {
        // B(2j) / (2j (2j - 1)) * ((M + m)^(1 - 2j) - M^(1 - 2j))
        mpfr_ui_div(t, 1, a, MPFR_RNDN);
        mpfr_pow_ui(t, t, static_cast<unsigned long>(2 * j - 1), MPFR_RNDN);
        mpfr_ui_div(b, 1, big, MPFR_RNDN);
        mpfr_pow_ui(b, b, static_cast<unsigned long>(2 * j - 1), MPFR_RNDN);
        mpfr_sub(t, t, b, MPFR_RNDN);
        mpfr_mul_si(t, t, num[j - 1], MPFR_RNDN);
        mpfr_div_ui(t, t, static_cast<unsigned long>(den[j - 1] * 2 * j * (2 * j - 1)), MPFR_RNDN);
        mpfr_add(sum, sum, t, MPFR_RNDN);
    }
    mpfr_gamma(t, m, MPFR_RNDN);
    mpfr_neg(sum, sum, MPFR_RNDN);
    mpfr_exp(sum, sum, MPFR_RNDN);
    mpfr_mul(t, t, sum, MPFR_RNDN);
    mpfr_set_prec(o, prec);
    mpfr_set(o, t, MPFR_RNDN);
    mpfr_clears(m, big, sum, a, b, t, static_cast<mpfr_ptr>(nullptr));
}
template <class T> void check_near()
{
    g_precision = 1500;
    auto E = edges<T>();
    big_float pi, t, mx, mn, dn;
    mpfr_const_pi(pi.v, MPFR_RNDN);
    to_mpfr(std::numeric_limits<T>::max(), mx.v);
    to_mpfr(std::numeric_limits<T>::min(), mn.v);
    to_mpfr(std::numeric_limits<T>::denorm_min(), dn.v);
    auto addm = [&](std::vector<T> &v, mpfr_t q) { add_around(v, from_mpfr<T>(q)); };
    // multiples of pi/2 and pi/4
    std::vector<T> Tr = E;
    for (long n :
         {1L,  2L,  3L,  4L,   5L,   6L,    7L,    8L,      9L,       10L,         11L,         12L,           13L,
          16L, 31L, 32L, 100L, 101L, 1000L, 1001L, 100000L, 1000001L, 1000000000L, 1000000001L, 1000000000000L})
    {
        mpfr_mul_si(t.v, pi.v, n, MPFR_RNDN);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        addm(Tr, t.v);
        mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
        addm(Tr, t.v);
    }
    check_unary<T>("sin", Tr, DEC1(sin), MPFR1(mpfr_sin), 3000);
    check_unary<T>("cos", Tr, DEC1(cos), MPFR1(mpfr_cos), 3000);
    check_unary<T>("tan", Tr, DEC1(tan), MPFR1(mpfr_tan), 3000);
    g_precision = 1500;
    std::vector<T> Inv = E;
    for (double d : {0.5, 0.7071067811865476, 0.8660254037844386, 0.9, 0.99, 0.999, 0.9999})
    {
        mpfr_set_d(t.v, d, MPFR_RNDN);
        addm(Inv, t.v);
    }
    mpfr_const_pi(t.v, MPFR_RNDN);
    mpfr_div_ui(t.v, t.v, 4, MPFR_RNDN);
    mpfr_sin(t.v, t.v, MPFR_RNDN);
    addm(Inv, t.v);
    check_unary<T>("asin", Inv, DEC1(asin), MPFR1(mpfr_asin), 1500);
    check_unary<T>("acos", Inv, DEC1(acos), MPFR1(mpfr_acos), 1500);
    check_unary<T>("atan", Inv, DEC1(atan), MPFR1(mpfr_atan), 1500);
    check_unary<T>("asinh", Inv, DEC1(asinh), MPFR1(mpfr_asinh), 1500);
    check_unary<T>("atanh", Inv, DEC1(atanh), MPFR1(mpfr_atanh), 1500);
    check_unary<T>("acosh", Inv, DEC1(acosh), MPFR1(mpfr_acosh), 1500);
    g_precision = 1500;
    std::vector<T> Lg = E;
    mpfr_set_ui(t.v, 1, MPFR_RNDN);
    mpfr_exp(t.v, t.v, MPFR_RNDN);
    addm(Lg, t.v);
    mpfr_sub_ui(t.v, t.v, 1, MPFR_RNDN);
    addm(Lg, t.v);
    add_around(Lg, std::numeric_limits<T>::max() / T{2});
    add_around(Lg, std::numeric_limits<T>::min() * T{1000});
    check_unary<T>("log", Lg, DEC1(log), MPFR1(mpfr_log), 1500);
    check_unary<T>("log2", Lg, DEC1(log2), MPFR1(mpfr_log2), 1500);
    check_unary<T>("log10", Lg, DEC1(log10), MPFR1(mpfr_log10), 1500);
    {
        std::vector<T> L1;
        for (T q : Lg)
            if (abs(q) > T{1, -100})
                L1.push_back(q);
        check_unary<T>("log1p", L1, DEC1(log1p), MPFR1(mpfr_log1p), 1500);
    }
    g_precision = 1500;
    std::vector<T> Ex = E;
    {
        big_float a;
        mpfr_log(a.v, mx.v, MPFR_RNDN);
        addm(Ex, a.v);
        mpfr_log(a.v, mn.v, MPFR_RNDN);
        addm(Ex, a.v);
        mpfr_log(a.v, dn.v, MPFR_RNDN);
        addm(Ex, a.v);
        mpfr_mul_ui(t.v, mx.v, 2, MPFR_RNDN);
        mpfr_log(a.v, t.v, MPFR_RNDN);
        addm(Ex, a.v);
        mpfr_set_ui(t.v, 2, MPFR_RNDN);
        mpfr_log(a.v, t.v, MPFR_RNDN);
        addm(Ex, a.v);
        mpfr_set_ui(t.v, 10, MPFR_RNDN);
        mpfr_log(a.v, t.v, MPFR_RNDN);
        addm(Ex, a.v);
    }
    check_unary<T>("exp", Ex, DEC1(exp), MPFR1(mpfr_exp), 1500);
    check_unary<T>("expm1", Ex, DEC1(expm1), MPFR1(mpfr_expm1), 1500);
    check_unary<T>("sinh", Ex, DEC1(sinh), MPFR1(mpfr_sinh), 1500);
    check_unary<T>("cosh", Ex, DEC1(cosh), MPFR1(mpfr_cosh), 1500);
    check_unary<T>("tanh", Ex, DEC1(tanh), MPFR1(mpfr_tanh), 1500);
    g_precision = 1500;
    std::vector<T> E2 = E;
    {
        big_float a;
        mpfr_log2(a.v, mx.v, MPFR_RNDN);
        addm(E2, a.v);
        mpfr_log2(a.v, mn.v, MPFR_RNDN);
        addm(E2, a.v);
        mpfr_log2(a.v, dn.v, MPFR_RNDN);
        addm(E2, a.v);
        mpfr_set_ui(t.v, 10, MPFR_RNDN);
        mpfr_log2(a.v, t.v, MPFR_RNDN);
        addm(E2, a.v);
    }
    check_unary<T>("exp2", E2, DEC1(exp2), MPFR1(mpfr_exp2), 1500);
    g_precision = 1500;
    std::vector<T> Er = E;
    {
        big_float a;
        auto f = [&](mpfr_t o, mpfr_t q)
        {
            mpfr_erfc(o, q, MPFR_RNDN);
            mpfr_sub(o, o, mn.v, MPFR_RNDN);
        };
        if (bisect(f, 1, 1000, a.v))
            addm(Er, a.v);
        auto f2 = [&](mpfr_t o, mpfr_t q)
        {
            mpfr_erfc(o, q, MPFR_RNDN);
            mpfr_sub(o, o, dn.v, MPFR_RNDN);
        };
        if (bisect(f2, 1, 1000, a.v))
            addm(Er, a.v);
        for (int i = 1; i <= 16; ++i)
            add_around(Er, T{i} / T{2});
    }
    check_unary<T>("erf", Er, DEC1(erf), MPFR1(mpfr_erf), 3000, 1e4);
    check_unary<T>("erfc", Er, DEC1(erfc), MPFR1(mpfr_erfc), 3000, 1e4);
    g_precision = 1500;
    std::vector<T> Gm = E;
    for (int n = -60; n <= 20; ++n)
        if (n != 0)
            add_around(Gm, T{n});
    for (int n = -20; n <= 20; ++n)
        add_around(Gm, T{2 * n + 1} / T{2});
    {
        big_float a;
        auto f = [&](mpfr_t o, mpfr_t q)
        {
            mpfr_gamma(o, q, MPFR_RNDN);
            mpfr_sub(o, o, mx.v, MPFR_RNDN);
        };
        if (bisect(f, 3, 3000, a.v))
            addm(Gm, a.v);
    }
    check_unary<T>("tgamma", Gm, DEC1(tgamma), MPFR1(mpfr_gamma), 1500, 5000);
    g_precision = 1500;
    std::vector<T> Lz = E;
    for (int n = -60; n <= 4; ++n)
        if (n != 0)
            add_around(Lz, T{n});
    {
        big_float a;
        auto f = [&](mpfr_t o, mpfr_t q)
        {
            int s;
            mpfr_lgamma(o, &s, q, MPFR_RNDN);
        };
        for (double l = -10; l < 0; l += 0.25)
            if (bisect(f, l + 1e-9, l + 0.25 - 1e-9, a.v))
                addm(Lz, a.v);
    }
    check_unary<T>(
        "lgamma", Lz, DEC1(lgamma),
        [](mpfr_t o, mpfr_t a)
        {
            int s;
            mpfr_lgamma(o, &s, a, MPFR_RNDN);
        },
        1500, 1e6);
    g_precision = 1000;
    std::vector<T> Zt = E;
    for (int n = -60; n <= 4; n += 2)
        add_around(Zt, T{n});
    add_around(Zt, T{1});
    add_around(Zt, T{3});
    {
        for (int j = 1; j < std::numeric_limits<T>::digits10; j += 3)
        {
            add_around(Zt, T{1} + T{1, -j});
            add_around(Zt, T{1} - T{1, -j});
        }
    }
    check_unary<T>("zeta", Zt, DEC1(riemann_zeta), MPFR1(mpfr_zeta), 1000, 200);
    g_precision = 1500;
    long bad = 0, cnt = 0;
    for (unsigned n : {2u, 3u, 5u, 8u, 16u, 40u})
    {
        const double H = std::sqrt(2.0 * n) + 3;
        roots<T>("hermite", 0, n, 0, -H, H, [n](T x) { return hermite(n, x); }, bad, cnt);
        roots<T>("legendre", 1, n, 0, -1, 1, [n](T x) { return legendre(n, x); }, bad, cnt);
        roots<T>("laguerre", 2, n, 0, 0, 4.0 * n + 6, [n](T x) { return laguerre(n, x); }, bad, cnt);
        for (unsigned mm : {1u, 2u, 5u})
        {
            roots<T>(
                "assoc_laguerre", 3, n, mm, 0, 4.0 * n + 4.0 * mm + 6,
                [n, mm](T x) { return assoc_laguerre(n, mm, x); }, bad, cnt);
            if (mm < n)
                roots<T>(
                    "assoc_legendre", 4, n, mm, -1, 1, [n, mm](T x) { return assoc_legendre(n, mm, x); }, bad, cnt);
        }
    }
    printf("poly-roots n=%ld bad=%ld\n", cnt, bad);
    fflush(stdout);
    g_bad_total += bad;
        // beta for one argument at or above 10^15 and the other below 10^4
    {
        g_precision = 1500;
        const std::vector<T> A = small_grid<T>(-6, 3, true);
        std::vector<T> L;
        for (int e : {15, 16, 20, 33, 60, 90})
        {
            add_around(L, T{1, e});
            add_around(L, T{7, e});
        }
        L.push_back(T{314159, 12});
        L.push_back(T{271828, 40});
        auto fd = [](T a, T b) { return beta(a, b); };
        check_binary<T>("beta_large", A, L, fd, large_beta_value, 1);
        check_binary<T>("beta_large", L, A, fd, large_beta_value, 1);
    }
    // beta near integers, atan2 near y = x, ellint near the edges
    {
        g_precision = 1000;
        std::vector<T> A;
        for (int n : {1, 2, 3, 10})
            add_around(A, T{n});
        add_around(A, T{1, -1});
        check_binary<T>(
            "beta_n", A, A, [](T a, T b) { return beta(a, b); },
            [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_beta(o, a, b, MPFR_RNDN); }, 1);
        check_binary<T>(
            "atan2_n", A, A, [](T a, T b) { return atan2(a, b); },
            [](mpfr_t o, mpfr_t a, mpfr_t b) { mpfr_atan2(o, a, b, MPFR_RNDN); }, 4);
    }
    {
        g_precision = 700;
        std::vector<T> K, Ph;
        add_around(K, T{1});
        add_around(K, T{9, -1});
        add_around(K, T{866, -3});
        for (int j = 1; j < std::numeric_limits<T>::digits10; j += 4)
            add_around(K, T{1} - T{1, -j});
        for (long n : {1L, 2L, 3L, 7L})
        {
            mpfr_const_pi(t.v, MPFR_RNDN);
            mpfr_mul_si(t.v, t.v, n, MPFR_RNDN);
            mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
            Ph.push_back(from_mpfr<T>(t.v));
            add_around(Ph, from_mpfr<T>(t.v));
        }
        std::vector<T> K2;
        for (size_t i = 0; i < K.size(); i += 2)
            K2.push_back(K[i]);
        std::vector<T> P2;
        for (size_t i = 0; i < Ph.size(); i += 3)
            P2.push_back(Ph[i]);
        long bF = 0, bE = 0, n = 0, bC = 0;
        big_float a, b, oF, oE;
        g_overflow = g_underflow = false;
        for (auto k : K2)
            for (auto p : P2)
            {
                to_mpfr(k, a.v);
                to_mpfr(p, b.v);
                ellint_values(a.v, b.v, oF.v, oE.v);
                std::string as = show(k) + "," + show(p);
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T f = ellint_1(k, p);
                    margin_events(oF.v, "ellint_1", as);
                    T e = ellint_2(k, p);
                    margin_events(oE.v, "ellint_2", as);
                    fesetround(g_modes[0]);
                    ++n;
                    check(f, oF.v, m, "ellint_1", as, bF);
                    check(e, oE.v, m, "ellint_2", as, bE);
                }
            }
        printf("ellint_n n=%ld badF=%ld badE=%ld\n", n, bF, bE);
        g_bad_total += bF + bE;
        n = 0;
        for (auto k : K)
        {
            to_mpfr(k, a.v);
            comp_ellint_values(a.v, oF.v, oE.v);
            std::string as = show(k);
            for (int m = 0; m < 5; ++m)
            {
                fesetround(g_modes[m]);
                T f = comp_ellint_1(k);
                margin_events(oF.v, "comp_ellint_1", as);
                T e = comp_ellint_2(k);
                margin_events(oE.v, "comp_ellint_2", as);
                fesetround(g_modes[0]);
                ++n;
                check(f, oF.v, m, "comp_ellint_1", as, bC);
                check(e, oE.v, m, "comp_ellint_2", as, bC);
            }
        }
        printf("comp_ellint_n n=%ld bad=%ld\n", n, bC);
        g_bad_total += bC;
    }
}
// The reference of a function at the arguments of a line of a list, and in ternary the sign of the reference minus
// the exact value, or 2 where it is not known
static bool list_value(const std::string &f, mpfr_t o, mpfr_t a, mpfr_t b, int &ternary)
{
    ternary = 2;
    static const std::pair<const char *, int (*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)> one[] = {
        {"exp", mpfr_exp},     {"expm1", mpfr_expm1}, {"exp2", mpfr_exp2},   {"log", mpfr_log},
        {"log2", mpfr_log2},   {"log10", mpfr_log10}, {"log1p", mpfr_log1p}, {"sin", mpfr_sin},
        {"cos", mpfr_cos},     {"tan", mpfr_tan},     {"asin", mpfr_asin},   {"acos", mpfr_acos},
        {"atan", mpfr_atan},   {"atanh", mpfr_atanh}, {"sinh", mpfr_sinh},   {"cosh", mpfr_cosh},
        {"tanh", mpfr_tanh},   {"asinh", mpfr_asinh}, {"acosh", mpfr_acosh}, {"erf", mpfr_erf},
        {"erfc", mpfr_erfc},   {"tgamma", mpfr_gamma}, {"lgamma", mpfr_lngamma}, {"zeta", mpfr_zeta}};
    for (const auto &e : one)
    {
        if (f == e.first)
        {
            const int i = e.second(o, a, MPFR_RNDN);
            ternary = (i > 0) - (i < 0);
            return true;
        }
    }
    if (f == "beta" || f == "beta_large")
    {
        if (mpfr_cmp_d(a, 1e15) >= 0 || mpfr_cmp_d(b, 1e15) >= 0)
            large_beta_value(o, a, b);
        else
        {
            const int i = mpfr_beta(o, a, b, MPFR_RNDN);
            ternary = (i > 0) - (i < 0);
        }
        return true;
    }
    if (f == "pow")
    {
        const int i = mpfr_pow(o, a, b, MPFR_RNDN);
        ternary = (i > 0) - (i < 0);
        return true;
    }
    if (f == "atan2")
    {
        const int i = mpfr_atan2(o, a, b, MPFR_RNDN);
        ternary = (i > 0) - (i < 0);
        return true;
    }
    if (f == "ellint_1" || f == "ellint_2")
    {
        big_float other;
        if (f == "ellint_1")
            ellint_values(a, b, o, other.v);
        else
            ellint_values(a, b, other.v, o);
        return true;
    }
    return false;
}
template <class T> T list_result(const std::string &f, T a, T b)
{
    using F = T (*)(T);
    static const std::pair<const char *, F> one[] = {
        {"exp", [](T x) -> T { return exp(x); }},       {"expm1", [](T x) -> T { return expm1(x); }},
        {"exp2", [](T x) -> T { return exp2(x); }},     {"log", [](T x) -> T { return log(x); }},
        {"log2", [](T x) -> T { return log2(x); }},     {"log10", [](T x) -> T { return log10(x); }},
        {"log1p", [](T x) -> T { return log1p(x); }},   {"sin", [](T x) -> T { return sin(x); }},
        {"cos", [](T x) -> T { return cos(x); }},       {"tan", [](T x) -> T { return tan(x); }},
        {"asin", [](T x) -> T { return asin(x); }},     {"acos", [](T x) -> T { return acos(x); }},
        {"atan", [](T x) -> T { return atan(x); }},     {"atanh", [](T x) -> T { return atanh(x); }},
        {"sinh", [](T x) -> T { return sinh(x); }},     {"cosh", [](T x) -> T { return cosh(x); }},
        {"tanh", [](T x) -> T { return tanh(x); }},     {"asinh", [](T x) -> T { return asinh(x); }},
        {"acosh", [](T x) -> T { return acosh(x); }},   {"erf", [](T x) -> T { return erf(x); }},
        {"erfc", [](T x) -> T { return erfc(x); }},     {"tgamma", [](T x) -> T { return tgamma(x); }},
        {"lgamma", [](T x) -> T { return lgamma(x); }}, {"zeta", [](T x) -> T { return riemann_zeta(x); }}};
    for (const auto &e : one)
    {
        if (f == e.first)
            return e.second(a);
    }
    if (f == "pow")
        return pow(a, b);
    if (f == "atan2")
        return atan2(a, b);
    if (f == "ellint_1")
        return ellint_1(a, b);
    if (f == "ellint_2")
        return ellint_2(a, b);
    return beta(a, b);
}
// The reference t against the digits of r: 0 away from a rounding boundary, 1 within the precision of one, and 2 on
// one. The boundaries are the multiples of a twentieth of the unit of r, which takes in those below a power of ten.
template <class T> int boundary_state(T r, mpfr_t t)
{
    if (!isfinite(r) || r == T{0} || !mpfr_regular_p(t))
        return 0;
    static big_float d(check_precision), u(check_precision);
    int e;
    frexp10(r, &e);
    to_mpfr(r, d.v);
    mpfr_sub(d.v, d.v, t, MPFR_RNDN);
    mpfr_ui_pow_ui(u.v, 10, e < 0 ? -e : e, MPFR_RNDN);
    if (e < 0)
        mpfr_mul(d.v, d.v, u.v, MPFR_RNDN);
    else
        mpfr_div(d.v, d.v, u.v, MPFR_RNDN);
    mpfr_mul_ui(d.v, d.v, 20, MPFR_RNDN);
    mpfr_rint(u.v, d.v, MPFR_RNDN);
    mpfr_sub(d.v, d.v, u.v, MPFR_RNDN);
    // r itself is exact only to check_precision, and t has fewer bits
    if (mpfr_zero_p(d.v) || mpfr_get_exp(d.v) < -(check_precision - 400))
        return 2;
    return mpfr_get_exp(d.v) < -(g_precision - 310) ? 1 : 0;
}
// t on a boundary of r moves past the window of check to the side of the exact value, for the sign side of t minus it
template <class T> void move_off_boundary(T r, mpfr_t t, int side)
{
    static big_float u(check_precision), s(check_precision);
    int e;
    frexp10(r, &e);
    mpfr_ui_pow_ui(s.v, 10, e < 0 ? -e : e, MPFR_RNDN);
    mpfr_set_ui_2exp(u.v, 1, -(g_precision - 320), MPFR_RNDN);
    if (e < 0)
        mpfr_div(u.v, u.v, s.v, MPFR_RNDN);
    else
        mpfr_mul(u.v, u.v, s.v, MPFR_RNDN);
    if (side > 0)
        mpfr_sub(t, t, u.v, MPFR_RNDN);
    else
        mpfr_add(t, t, u.v, MPFR_RNDN);
}
// Each line of the file is function|a|b|mode. The program prints the same line with the result and a mark if wrong.
// The precision holds a large argument exactly, and doubles up to 60000 bits next to a boundary MPFR leaves open.
template <class T> void check_list(const char *path)
{
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
        const auto p1 = line.find('|'), p2 = line.find('|', p1 + 1), p3 = line.find('|', p2 + 1);
        const std::string f = line.substr(0, p1);
        const T x = list_argument<T>(line.substr(p1 + 1, p2 - p1 - 1));
        const T y = list_argument<T>(line.substr(p2 + 1, p3 - p2 - 1));
        const int m = std::stoi(line.substr(p3 + 1));
        fesetround(g_modes[m]);
        const T r = list_result<T>(f, x, y);
        fesetround(g_modes[0]);
        int start = 1500;
        for (const T v : {x, y})
        {
            int e = 0;
            if (isfinite(v))
                frexp10(v, &e);
            start = e * 7 / 3 + 300 > start ? e * 7 / 3 + 300 : start;
        }
        for (g_precision = start;; g_precision *= 2)
        {
            big_float a, b, t;
            to_mpfr(x, a.v);
            to_mpfr(y, b.v);
            mpfr_clear_flags();
            int ternary = 2;
            if (!list_value(f, t.v, a.v, b.v, ternary))
            {
                printf("%s|unknown\n", line.c_str());
                ++g_bad_total;
                break;
            }
            const int state = boundary_state(r, t.v);
            if (state == 2 && ternary != 2)
            {
                if (ternary != 0)
                    move_off_boundary(r, t.v, ternary);
            }
            else if (state != 0 && g_precision * 2 <= check_precision)
                continue;
            g_overflow = mpfr_overflow_p();
            g_underflow = mpfr_underflow_p();
            long bad = 0;
            check(r, t.v, m, f.c_str(), line, bad);
            printf("%s|%s%s\n", line.c_str(), show(r).c_str(), bad != 0 ? "|wrong" : "");
            g_bad_total += bad;
            break;
        }
    }
}
int main(int argc, char **argv)
{
    if (argc != 3 && !(argc == 4 && std::string(argv[2]) == "list"))
    {
        std::puts("usage: rounding_grid d32|d64|d128|f32|f64|f128 main|exact|near|ext|ell|assoc|list <file>");
        return 2;
    }
    std::string t = argv[1], w = argv[2];
    const bool known_type = t == "d32" || t == "d64" || t == "d128" || t == "f32" || t == "f64" || t == "f128";
    const bool known_part =
        w == "main" || w == "exact" || w == "near" || w == "ext" || w == "ell" || w == "assoc" || w == "list";
    if (!known_type || !known_part)
    {
        std::puts("usage: rounding_grid d32|d64|d128|f32|f64|f128 main|exact|near|ext|ell|assoc|list <file>");
        return 2;
    }
#ifdef BOUND_MARGIN
    boost_mut::record_events = true;
#endif
#define ALL(T, n)                                                                                                      \
    if (t == n)                                                                                                        \
    {                                                                                                                  \
        if (w == "ell")                                                                                                \
            check_ellint<T>(n);                                                                                        \
        if (w == "assoc")                                                                                              \
            check_assoc<T>(n);                                                                                         \
        if (w == "ext")                                                                                                \
            check_extremes<T>(n);                                                                                      \
        if (w == "main")                                                                                               \
            check_main<T>(n);                                                                                          \
        if (w == "exact")                                                                                              \
            check_exact<T>();                                                                                          \
        if (w == "near")                                                                                               \
            check_near<T>();                                                                                           \
        if (w == "list")                                                                                               \
            check_list<T>(argv[3]);                                                                                    \
    }
    ALL(decimal32_t, "d32")
    ALL(decimal64_t, "d64") ALL(decimal128_t, "d128") ALL(decimal_fast32_t, "f32") ALL(decimal_fast64_t, "f64")
        ALL(decimal_fast128_t, "f128") printf("%s %s total bad %ld\n", t.c_str(), w.c_str(), g_bad_total);
    return g_bad_total == 0 ? 0 : 1;
}
