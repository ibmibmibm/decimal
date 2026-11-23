// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Compares the <cmath> functions against MPFR at fixed grids of arguments, for one type and in all five rounding modes.
// See README.md in this directory.

#include <boost/decimal.hpp>
#include "record_field.hpp"
#include <mpfr.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <array>
#include <type_traits>
#include <utility>
#ifdef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
#error "The rounding modes need a compiler that lets the library detect constant evaluation"
#endif
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
template <class T> void to_mpfr(T x, mpfr_t o, mpfr_rnd_t rnd = MPFR_RNDN)
{
    if (isnan(x))
    {
        mpfr_set_nan(o);
        return;
    }
    if (isinf(x))
    {
        mpfr_set_inf(o, signbit(x) ? -1 : 1);
        return;
    }
    int e;
    auto s = frexp10(x, &e);
    if (s == 0)
    {
        mpfr_set_zero(o, signbit(x) ? -1 : 1);
        return;
    }
    // s 10^e with one rounding from exact integers, as mpfr_set_str, which is slow at a large precision
    static mpz_t z, p;
    static big_float n(128);
    static const bool init = (mpz_inits(z, p, nullptr), true);
    static_cast<void>(init);
    mpz_set_str(z, ustr(s).c_str(), 10);
    if (signbit(x))
        mpz_neg(z, z);
    mpz_ui_pow_ui(p, 10, static_cast<unsigned long>(e < 0 ? -e : e));
    if (e >= 0)
    {
        mpz_mul(z, z, p);
        mpfr_set_z(o, z, rnd);
        return;
    }
    mpfr_set_z(n.v, z, MPFR_RNDN);
    mpfr_div_z(o, n.v, p, rnd);
}
// The exponent of the smallest subnormal value, and the exponent of the unit of r, which is that one below min()
template <class T> int tiny_exponent()
{
    int e;
    auto q = frexp10(std::numeric_limits<T>::denorm_min(), &e);
    for (q /= 10; q != 0; q /= 10)
        ++e;
    return e;
}
// The unit of an infinite r is that of the largest value, which has the overflow threshold as a boundary
template <class T> int unit_exponent(T r)
{
    int e = tiny_exponent<T>();
    if (abs(r) >= std::numeric_limits<T>::min())
        frexp10(isinf(r) ? std::numeric_limits<T>::max() : r, &e);
    return e;
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
// The digits of a field [-]d[.d][e[-]d] with no zeros at the two ends, and the exponent of the last digit
static bool field_digits(const std::string &s, std::string &d, long &e)
{
    std::size_t i = !s.empty() && s[0] == '-' ? 1 : 0;
    long frac = 0;
    bool dot = false, any = false;
    for (; i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || (s[i] == '.' && !dot)); ++i)
    {
        if (s[i] == '.')
        {
            dot = true;
            continue;
        }
        any = true;
        if (!d.empty() || s[i] != '0')
            d += s[i];
        frac += dot ? 1 : 0;
    }
    e = 0;
    if (any && i < s.size() && (s[i] == 'e' || s[i] == 'E'))
    {
        // An exponent needs digits, and one of more than 9 digits, which a long of 32 bits cannot hold, is not valid
        const bool eneg = i + 1 < s.size() && s[i + 1] == '-';
        const std::size_t j = i + 1 + (i + 1 < s.size() && (s[i + 1] == '-' || s[i + 1] == '+') ? 1 : 0);
        for (i = j; i < s.size() && i < j + 9 && s[i] >= '0' && s[i] <= '9'; ++i)
            e = e * 10 + (s[i] - '0');
        e = eneg ? -e : e;
        if (i == j)
            return false;
    }
    e -= frac;
    for (; !d.empty() && d.back() == '0'; ++e)
        d.pop_back();
    return any && i == s.size();
}
template <class T> bool list_argument(const std::string &s, T &x)
{
    if (parse_special(s, x))
        return true;
    const bool neg = !s.empty() && s[0] == '-';
    const char *last = s.data() + s.size();
    const auto r = from_chars(s.data() + (neg ? 1 : 0), last, x);
    if (r.ec != std::errc() || r.ptr != last)
        return false;
    x = neg ? -x : x;
    // The decimal must hold the value of the field: its digits and exponent are compared, with no binary conversion.
    // A second sign, too many digits or an overflow is not valid.
    std::string want, have;
    long ew = 0, eh = 0;
    return field_digits(s, want, ew) && isfinite(x) && field_digits(show(x), have, eh) && want == have &&
           (want.empty() || ew == eh);
}
#define MPFR1(f) [](mpfr_t o, mpfr_t a) { return f(o, a, MPFR_RNDN); }
#define DEC1(f) [](T x) { return f(x); }
#ifdef BOUND_MARGIN
#include "bound_margin.hpp"
#else
static void margin_events(mpfr_t, const char *, const std::string &)
{
}
static bool margin_report()
{
    return true;
}
#endif
// With GRID_REC set to a file, the file gets a line for each check: name|arguments|mode|result
static FILE *g_rec = std::getenv("GRID_REC") ? std::fopen(std::getenv("GRID_REC"), "w") : nullptr;
static long g_bad_total;
// g_rec_failed says that the file of GRID_REC did not open or did not get all records
static bool g_rec_failed;
const rounding_mode g_modes[] = {rounding_mode::fe_dec_to_nearest, rounding_mode::fe_dec_upward,
                                 rounding_mode::fe_dec_downward, rounding_mode::fe_dec_to_nearest_from_zero,
                                 rounding_mode::fe_dec_toward_zero};
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
// 10^|e| in p with the bits that hold it exactly, as a product or a quotient costs the limbs of each argument
static void pow10_exact(mpfr_t p, int e)
{
    const unsigned long k = static_cast<unsigned long>(e < 0 ? -e : e);
    mpfr_set_prec(p, static_cast<mpfr_prec_t>(k * 10 / 3 + 64));
    mpfr_ui_pow_ui(p, 10, k, MPFR_RNDN);
}
// Checks r against the value t in mode index mode, with zero, subnormal and overflow results.
// The direction dir is 0 for nearest, 1 for up and 2 for down.
static bool g_overflow, g_underflow;
// g_skipped says that check_value took a result of a fast type without a comparison
static bool g_skipped;
// The caller of check sets g_undecided where no precision gave the side of a boundary: such a result gets no record
static bool g_undecided;
static long g_inconclusive, g_skipped_total;
// g_say_undecided is false only for the test of self, which makes undecided results on purpose
static bool g_say_undecided = true;
template <class T> bool check_value(T r, mpfr_t t, int mode, const char *name, const std::string &arg, long &bad)
{
    margin_events(t, name, arg);
    static big_float d(check_precision), u(check_precision), mx(check_precision), mn(check_precision),
        a0(check_precision);
    // The largest and the smallest normal value, set one time: a conversion at this precision is slow
    static const bool limits_set = (to_mpfr(std::numeric_limits<T>::max(), mx.v),
                                    to_mpfr(std::numeric_limits<T>::min(), mn.v), true);
    static_cast<void>(limits_set);
    auto say = [&](const char *w, double du)
    {
        if (bad++ < 8)
            printf("  %s(%s) mode %d %s err %.6f got %s\n", name, arg.c_str(), mode, w, du, show(r).c_str());
        return false;
    };
    if (mpfr_nan_p(t))
        return isnan(r) ? true : say("NaN expected", 0);
    // At its pole riemann_zeta gives NaN, where MPFR gives infinity without an overflow
    if (mpfr_inf_p(t) && !g_overflow && std::strcmp(name, "zeta") == 0)
        return isnan(r) ? true : say("NaN at the pole expected", 0);
    if (isnan(r))
        return say("unexpected NaN", 0);
    const bool neg = mpfr_signbit(t);
    // effective direction: 0 nearest, 1 up, 2 down
    const int dir = mode == 1 ? 1 : mode == 2 ? 2 : mode == 4 ? (neg ? 1 : 2) : 0;
    if (mpfr_zero_p(t) && g_underflow)
    {
        // A fast type has no subnormal results: its result far below the range is a zero
        if (is_fast<T>::v)
            return (g_skipped = true) && r == T{0} ? true : say("fast underflow", 0);
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
    mpfr_abs(u.v, t, MPFR_RNDN);
    if (mpfr_inf_p(t) && !g_overflow)
    {
        return (isinf(r) && signbit(r) == neg) ? true : say("pole expected", 0);
    }
    if (mpfr_inf_p(t) || mpfr_cmp(u.v, mx.v) > 0)
    {
        // A fast type overflows as the other types do, but its result gets no record
        g_skipped = is_fast<T>::v;
        const T mxv = std::numeric_limits<T>::max();
        const bool toinf = dir == 0 || (dir == 1 && !neg) || (dir == 2 && neg);
        if (toinf)
        {
            // Nearest gives the largest value below the overflow threshold, which is half a unit above it
            int em;
            frexp10(mxv, &em);
            mpfr_set_ui(d.v, 10, MPFR_RNDN);
            mpfr_pow_si(d.v, d.v, em, MPFR_RNDN);
            mpfr_div_ui(d.v, d.v, 2, MPFR_RNDN);
            mpfr_add(d.v, d.v, mx.v, MPFR_RNDN);
            if (!mpfr_inf_p(t) && dir == 0 && mpfr_cmp(u.v, d.v) < 0)
                return r == (neg ? -mxv : mxv) ? true : say("max expected", 0);
            return isinf(r) && signbit(r) == neg ? true : say("inf expected", 0);
        }
        return r == (neg ? -mxv : mxv) ? true : say("max expected", 0);
    }
    if (isinf(r))
        return say("inf unexpected", 0);
    int e;
    auto s = frexp10(r, &e);
    const int etiny = tiny_exponent<T>();
    const bool sub = mpfr_cmp(u.v, mn.v) < 0 || r == T{0};
    const int er = e;
    if (sub)
    {
        // A fast type has no subnormal results, but a zero for a value in the normal range is a wrong result
        if (is_fast<T>::v && mpfr_cmp(u.v, mn.v) < 0)
            return (g_skipped = true) && abs(r) <= std::numeric_limits<T>::min() ? true : say("fast underflow", 0);
        e = etiny;
    }
    // A value that is not zero rounds to zero only toward zero
    if (r == T{0} && ((dir == 1 && !neg) || (dir == 2 && neg)))
        return say("zero away from zero", 0);
    to_mpfr(r, d.v);
    // Below a power of ten, the unit is 10 times smaller
    auto pw = s;
    pw = 1;
    for (int i = 1; i < std::numeric_limits<T>::digits10; ++i)
        pw *= 10;
    if (!sub && s == pw && e > etiny && mpfr_cmpabs(t, d.v) < 0)
        --e;
    mpfr_sub(d.v, d.v, t, MPFR_RNDN);
    static big_float p10;
    pow10_exact(p10.v, e);
    if (e < 0)
        mpfr_mul(d.v, d.v, p10.v, MPFR_RNDN);
    else
        mpfr_div(d.v, d.v, p10.v, MPFR_RNDN);
    double du = mpfr_get_d(d.v, MPFR_RNDN);
    static big_float tl(check_precision), w1(check_precision), w2(check_precision);
    bool ok;
    mpfr_set_ui_2exp(tl.v, 1, -(g_precision - 300), MPFR_RNDN);
    // A value far below the smallest subnormal value rounds to it only away from zero
    if (sub && abs(r) == std::numeric_limits<T>::denorm_min() && std::fabs(std::fabs(du) - 1) < 1e-9 &&
        signbit(r) == neg && (neg ? dir == 2 && du < 0 : dir == 1 && du > 0))
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
    // A zero from a value that is not zero has the sign of the value
    if (r == T{0} && signbit(r) != neg)
        return say("zero sign", du);
    if (dir == 0)
    {
        mpfr_abs(w1.v, d.v, MPFR_RNDN);
        mpfr_sub_d(w1.v, w1.v, 0.5, MPFR_RNDN);
        mpfr_abs(w1.v, w1.v, MPFR_RNDN);
        if (mpfr_cmp(w1.v, tl.v) < 0)
        { // an exact tie
            if (mode == 0)
            {
                // A subnormal result has fewer digits
                for (int i = er; sub && i < etiny; ++i)
                    s /= 10;
                if (static_cast<int>(s % 2U) != 0)
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
template <class T> void move_off_boundary(T r, mpfr_t t, int side);
// check_value, and a record in the file of GRID_REC for a result that the comparison accepts
template <class T> bool check(T r, mpfr_t t, int mode, const char *name, const std::string &arg, long &bad)
{
    g_skipped = false;
    bool ok;
    if (g_undecided)
    {
        // The value has an unknown side of the boundary t: the result is wrong only if it is wrong on each side
        static big_float v(check_precision);
        long quiet = 8;
        ok = check_value(r, t, mode, name, arg, quiet);
        for (const int side : {1, -1})
        {
            mpfr_set(v.v, t, MPFR_RNDN);
            if (!ok)
                move_off_boundary(r, v.v, side);
            ok = ok || check_value(r, v.v, mode, name, arg, quiet);
        }
        if (!ok)
            check_value(r, t, mode, name, arg, bad);
    }
    else
        ok = check_value(r, t, mode, name, arg, bad);
    if (g_rec && ok && !g_skipped && !g_undecided)
        g_rec_failed |= std::fprintf(g_rec, "%s|%s|%d|%s\n", name, arg.c_str(), mode, show(r).c_str()) < 0;
    if (g_undecided && g_say_undecided && g_inconclusive < 8)
        printf("  %s(%s) mode %d inconclusive got %s\n", name, arg.c_str(), mode, show(r).c_str());
    g_inconclusive += g_undecided;
    g_skipped_total += g_skipped;
    g_undecided = false;
    return ok;
}
template <class T> int boundary_state(T r, mpfr_t t);
static bool field_digits(const std::string &s, std::string &d, long &e);
// The decimal x as a rational in o
template <class T> void decimal_q(T x, mpq_t o)
{
    std::string d;
    long k = 0;
    field_digits(show(x), d, k);
    mpz_t p;
    mpz_init(p);
    mpz_ui_pow_ui(p, 10, static_cast<unsigned long>(k < 0 ? -k : k));
    mpq_set_ui(o, 0, 1);
    mpz_set_str(mpq_numref(o), d.empty() ? "0" : d.c_str(), 10);
    if (k < 0)
        mpz_set(mpq_denref(o), p);
    else
        mpz_mul(mpq_numref(o), mpq_numref(o), p);
    mpq_canonicalize(o);
    if (signbit(x))
        mpq_neg(o, o);
    mpz_clear(p);
}
// check for a reference t that fr computes at the precision g_precision. Next to a rounding boundary, fr runs again
// with four times the bits. Where they do not give the side of an inexact reference, ex can give the exact value as a
// rational, which decides the side, else the result is inconclusive.
template <class T, class FR, class EX>
bool check_ref(T r, mpfr_t t, int mode, const char *name, const std::string &arg, long &bad, FR fr, EX ex)
{
    if (boundary_state(r, t) == 0)
        return check(r, t, mode, name, arg, bad);
    const int prec = g_precision;
    g_precision *= 4;
    big_float t2;
    mpfr_clear_flags();
    fr(t2.v);
    g_undecided = mpfr_inexflag_p() && boundary_state(r, t2.v) != 0;
    mpq_t q;
    mpq_init(q);
    if (g_undecided && ex(q))
    {
        // The boundary is c 10^g, a multiple of a twentieth of the unit of r, as in list_exact
        static big_float u(check_precision), p10;
        const int e = unit_exponent(r);
        pow10_exact(p10.v, e);
        if (e < 0)
            mpfr_mul(u.v, t2.v, p10.v, MPFR_RNDN);
        else
            mpfr_div(u.v, t2.v, p10.v, MPFR_RNDN);
        mpfr_mul_ui(u.v, u.v, 20, MPFR_RNDN);
        mpq_t b;
        mpq_init(b);
        mpfr_get_z(mpq_numref(b), u.v, MPFR_RNDN);
        mpz_mul_ui(mpq_numref(b), mpq_numref(b), 5);
        mpz_ui_pow_ui(mpq_denref(b), 10, static_cast<unsigned long>(e - 2 < 0 ? 2 - e : e - 2));
        if (e - 2 >= 0)
        {
            mpz_mul(mpq_numref(b), mpq_numref(b), mpq_denref(b));
            mpz_set_ui(mpq_denref(b), 1);
        }
        mpq_canonicalize(b);
        const int side = mpq_cmp(b, q);
        mpfr_set_prec(t2.v, check_precision);
        mpfr_set_q(t2.v, b, MPFR_RNDN);
        if (side != 0)
            move_off_boundary(r, t2.v, side > 0 ? 1 : -1);
        g_undecided = false;
        mpq_clear(b);
    }
    mpq_clear(q);
    const bool ok = check(r, t2.v, mode, name, arg, bad);
    g_precision = prec;
    return ok;
}
template <class T, class FR> bool check_ref(T r, mpfr_t t, int mode, const char *name, const std::string &arg, long &bad, FR fr)
{
    return check_ref(r, t, mode, name, arg, bad, fr, [](mpq_t) { return false; });
}
// P(x) of poly_value as a rational in o, for a rational x, false for assoc_legendre where sqrt(1 - x^2) is not one
static bool poly_value_q(int kind, unsigned n, unsigned mm, mpq_t x, mpq_t o)
{
    mpq_t p0, p1, p2, t1, k;
    mpq_inits(p0, p1, p2, t1, k, nullptr);
    bool known = true;
    auto num = [&](unsigned long a, unsigned long b = 1) { mpq_set_ui(k, a, b); mpq_canonicalize(k); return k; };
    mpq_set_ui(p0, 1, 1);
    if (kind == 4)
    {
        // P_mm^mm = (2mm - 1)!! (1 - x^2)^(mm/2), then the recurrence in l
        mpq_mul(t1, x, x);
        mpq_sub(t1, num(1), t1);
        const bool square = mpq_sgn(t1) >= 0 && mpz_perfect_square_p(mpq_numref(t1)) && mpz_perfect_square_p(mpq_denref(t1));
        if (mm > n)
            mpq_set_ui(o, 0, 1);
        else if (mm % 2 != 0 && !square)
            known = false;
        else
        {
            mpq_t c;
            mpq_init(c);
            mpz_sqrt(mpq_numref(c), mpq_numref(t1));
            mpz_sqrt(mpq_denref(c), mpq_denref(t1));
            for (unsigned i = 0; i < mm; ++i)
            {
                mpq_mul(p0, p0, num(2 * i + 1));
                mpq_mul(p0, p0, mm % 2 != 0 ? c : t1);
                if (mm % 2 == 0 && ++i < mm)
                    mpq_mul(p0, p0, num(2 * i + 1));
            }
            mpq_clear(c);
            mpq_mul(p1, x, p0);
            mpq_mul(p1, p1, num(2 * mm + 1));
            for (unsigned l = mm + 1; l < n; ++l)
            {
                mpq_mul(p2, x, p1);
                mpq_mul(p2, p2, num(2 * l + 1));
                mpq_mul(t1, p0, num(l + mm));
                mpq_sub(p2, p2, t1);
                mpq_div(p2, p2, num(l - mm + 1));
                mpq_swap(p0, p1);
                mpq_swap(p1, p2);
            }
            mpq_set(o, n == mm ? p0 : p1);
        }
    }
    else if (n == 0)
        mpq_set_ui(o, 1, 1);
    else
    {
        if (kind == 0)
            mpq_add(p1, x, x);
        else if (kind == 1)
            mpq_set(p1, x);
        else
            mpq_sub(p1, num(1 + mm), x);
        for (unsigned j = 1; j < n; ++j)
        {
            if (kind == 0)
            {
                mpq_mul(p2, x, p1);
                mpq_add(p2, p2, p2);
                mpq_mul(t1, p0, num(2 * j));
                mpq_sub(p2, p2, t1);
            }
            else if (kind == 1)
            {
                mpq_mul(p2, x, p1);
                mpq_mul(p2, p2, num(2 * j + 1));
                mpq_mul(t1, p0, num(j));
                mpq_sub(p2, p2, t1);
                mpq_div(p2, p2, num(j + 1));
            }
            else
            {
                mpq_sub(t1, num(2 * j + 1 + mm), x);
                mpq_mul(p2, t1, p1);
                mpq_mul(t1, p0, num(j + mm));
                mpq_sub(p2, p2, t1);
                mpq_div(p2, p2, num(j + 1));
            }
            mpq_swap(p0, p1);
            mpq_swap(p1, p2);
        }
        mpq_set(o, p1);
    }
    mpq_clears(p0, p1, p2, t1, k, nullptr);
    return known;
}
// The polynomial of poly_value_q at the decimal x
template <class T> bool poly_exact_q(int kind, unsigned n, unsigned mm, T x, mpq_t o)
{
    mpq_t q;
    mpq_init(q);
    decimal_q(x, q);
    const bool known = poly_value_q(kind, n, mm, q, o);
    mpq_clear(q);
    return known;
}
// x^k for the decimal x, false for x = 0 and k < 0
template <class T> bool pow_exact_q(T x, int k, mpq_t o)
{
    decimal_q(x, o);
    if (k < 0 && mpq_sgn(o) == 0)
        return false;
    const unsigned long a = static_cast<unsigned long>(k < 0 ? -static_cast<long>(k) : k);
    mpz_pow_ui(mpq_numref(o), mpq_numref(o), a);
    mpz_pow_ui(mpq_denref(o), mpq_denref(o), a);
    if (k < 0)
        mpq_inv(o, o);
    return true;
}
// hypot(x, y, z) where x^2 + y^2 + z^2 is the square of a rational
template <class T> bool hypot3_exact_q(T x, T y, T z, mpq_t o)
{
    mpq_t a, s;
    mpq_inits(a, s, nullptr);
    mpq_set_ui(s, 0, 1);
    for (const T v : {x, y, z})
    {
        decimal_q(v, a);
        mpq_mul(a, a, a);
        mpq_add(s, s, a);
    }
    const bool square = mpz_perfect_square_p(mpq_numref(s)) && mpz_perfect_square_p(mpq_denref(s));
    if (square)
    {
        mpz_sqrt(mpq_numref(o), mpq_numref(s));
        mpz_sqrt(mpq_denref(o), mpq_denref(s));
    }
    mpq_clears(a, s, nullptr);
    return square;
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
    if (mpfr_nan_p(k) || mpfr_nan_p(phi))
    {
        mpfr_set_nan(oF);
        mpfr_set_nan(oE);
        return;
    }
    int c = mpfr_cmp_ui(ak.v, 1);
    // A zero phase keeps its sign, and the library gives an infinity for an infinite phase
    if (c <= 0 && (mpfr_zero_p(phi) || mpfr_inf_p(phi)))
    {
        mpfr_set(oF, phi, MPFR_RNDN);
        mpfr_set(oE, phi, MPFR_RNDN);
        return;
    }
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
    if (r > 0 || mpfr_nan_p(k))
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
                    check_ref(f, oF.v, m, "ellint_1", as, bF,
                              [&](mpfr_t o) { big_float a2, b2, e2; to_mpfr(k, a2.v); to_mpfr(p, b2.v); ellint_values(a2.v, b2.v, o, e2.v); },
                              [&](mpq_t o) { return k == T{0} && (decimal_q(p, o), true); });
                    check_ref(e, oE.v, m, "ellint_2", as, bE,
                              [&](mpfr_t o) { big_float a2, b2, f2; to_mpfr(k, a2.v); to_mpfr(p, b2.v); ellint_values(a2.v, b2.v, f2.v, o); },
                              [&](mpq_t o) { return k == T{0} && (decimal_q(p, o), true); });
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
                check_ref(e, oE.v, m, "comp_ellint_2", as, bC, [&](mpfr_t o) { big_float a2; to_mpfr(k, a2.v); comp_ellint_e(a2.v, o); });
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
                        check_ref(r, res.v, m, "assoc_laguerre", as, bL,
                                  [&](mpfr_t o) { big_float a2; to_mpfr(x, a2.v); poly_value(3, nn, mm, a2.v, o); },
                                  [&](mpq_t o) { return poly_exact_q(3, nn, mm, x, o); });
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
                        check_ref(r, res.v, m, "assoc_legendre", as, bP,
                                  [&](mpfr_t o) { big_float a2; to_mpfr(x, a2.v); poly_value(4, nn, mm, a2.v, o); },
                                  [&](mpq_t o) { return poly_exact_q(4, nn, mm, x, o); });
                    }
                }
    printf("%s assoc_legendre n=%ld bad=%ld\n", tn, n, bP);
    fflush(stdout);
    g_bad_total += bL + bP;
}
template <class T> int boundary_state(T r, mpfr_t t);
template <class T> bool list_exact(const std::string &f, T x, T y, T r, mpfr_t t);
// The reference fm of x in o, at the precision of a and o. The result is the sign of o minus the exact value, or 2
// where it is not known: fm does not give it, or a is not x.
template <class T, class FM> int reference(FM fm, mpfr_t o, mpfr_t a, T x)
{
    mpfr_clear_flags();
    to_mpfr(x, a);
    const bool exact = !mpfr_inexflag_p();
    mpfr_clear_flags();
    if constexpr (std::is_void<decltype(fm(o, a))>::value)
    {
        fm(o, a);
        return !exact || mpfr_inexflag_p() ? 2 : 0;
    }
    else
    {
        const int i = fm(o, a);
        return exact ? (i > 0) - (i < 0) : 2;
    }
}
// Whether f is monotonic in each argument on an interval that holds no exact binary value where it turns, such as 0 for
// cosh and 1 for x in pow. Then the values at the ends of the binary intervals of the arguments bound the exact value.
static bool monotonic(const std::string &name)
{
    const bool x = name.size() > 2 && name.compare(name.size() - 2, 2, "_x") == 0;
    const std::string f = x ? name.substr(0, name.size() - 2) : name;
    for (const char *m : {"exp", "exp2", "expm1", "log", "log2", "log10", "log1p", "sqrt", "cbrt", "sinh", "cosh", "tanh",
                          "asinh", "acosh", "atanh", "asin", "acos", "atan", "erf", "erfc", "comp_ellint_2", "pow",
                          "atan2", "hypot", "beta"})
        if (f == m)
            return true;
    return false;
}
// riemann_zeta falls on (-2.717, 1) and on (1, inf), where tiny and large x give values next to -1/2 and 1
template <class T> static bool monotonic(const std::string &name, T x)
{
    return monotonic(name) || (name == "zeta" && x > T{-27, -1});
}
// For an x with a binary error, the same sign where fm at both ends of the binary interval of x rounds to t from one
// side, for an fm that is monotonic in the interval
template <class T, class FM> int reference_side(FM fm, mpfr_t t, T x)
{
    int side = 2;
    if constexpr (!std::is_void<decltype(fm(t, t))>::value)
    {
        big_float a, o;
        for (const mpfr_rnd_t rnd : {MPFR_RNDD, MPFR_RNDU})
        {
            to_mpfr(x, a.v, rnd);
            const int i = fm(o.v, a.v);
            const int s = (i > 0) - (i < 0);
            if (!mpfr_equal_p(o.v, t) || s == 0 || (rnd == MPFR_RNDU && s != side))
                return 2;
            side = s;
        }
    }
    return side;
}
// Whether t decides the rounding of r. On a boundary, t moves to the side of the exact value from the ternary value.
template <class T> bool decides(T r, mpfr_t t, int ternary)
{
    const int state = boundary_state(r, t);
    if (state == 2 && ternary != 2)
    {
        if (ternary != 0)
            move_off_boundary(r, t, ternary);
        return true;
    }
    return state == 0;
}
// A reference that does not decide a result is taken again with two times the bits, up to check_precision. A result
// that no precision decides is on a boundary for check(), as an exact value of an argument with a binary error.
template <class T, class FD, class FM>
void check_unary(const char *name, const std::vector<T> &G, FD fd, FM fm, int prec, double lim = HUGE_VAL)
{
    long bad = 0, n = 0, skipped = 0;
    for (auto x0 : G)
        for (int sg = 0; sg < 2; ++sg)
        {
            T x = sg ? -x0 : x0;
            g_precision = prec;
            big_float a, t;
            to_mpfr(x, a.v);
            if (mpfr_cmp_d(a.v, lim) > 0 || mpfr_cmp_d(a.v, -lim) < 0)
            {
                ++skipped;
                continue;
            }
            int ternary = reference(fm, t.v, a.v, x);
            bool overflow = mpfr_overflow_p(), underflow = mpfr_underflow_p();
            std::string as = show(x);
            for (int m = 0; m < 5; ++m)
            {
                fesetround(g_modes[m]);
                T r = fd(x);
                fesetround(g_modes[0]);
                ++n;
                for (;;)
                {
                    big_float v;
                    mpfr_set(v.v, t.v, MPFR_RNDN);
                    g_overflow = overflow;
                    g_underflow = underflow;
                    if (ternary == 2 && monotonic(name, x) && boundary_state(r, v.v) == 2)
                        ternary = reference_side(fm, t.v, x);
                    bool decided = decides(r, v.v, ternary);
                    if (!decided && boundary_state(r, v.v) != 0 && list_exact(name, x, T{0}, r, v.v))
                        decided = true;
                    if (decided || g_precision * 2 > check_precision)
                    {
                        g_undecided = !decided;
                        check(r, v.v, m, name, as, bad);
                        break;
                    }
                    g_precision *= 2;
                    mpfr_set_prec(a.v, g_precision);
                    mpfr_set_prec(t.v, g_precision);
                    ternary = reference(fm, t.v, a.v, x);
                    overflow = mpfr_overflow_p();
                    underflow = mpfr_underflow_p();
                }
            }
        }
    g_precision = prec;
    printf("%s n=%ld bad=%ld skipped=%ld\n", name, n, bad, skipped);
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
    // Significands with all the digits of a type that has 34
    int count = 0;
    for (unsigned long long a : {10000000000000000ULL, 31415926535897932ULL, 49999999999999999ULL, 99999999999999999ULL})
        for (unsigned long long b : {1ULL, 38462643383279502ULL, 99999999999999999ULL})
            for (int e = lo; P == 34 && e <= hi; e += 3)
                if (!small || count++ % 5 == 0)
                    v.push_back(T{boost::int128::uint128_t{a} * 100000000000000000ULL + b, e - 33});
    return v;
}
template <class T, class FD, class FM>
void check_binary(const char *name, const std::vector<T> &H, const std::vector<T> &K, FD fd, FM fm, int sgs = 4)
{
    long bad = 0, n = 0;
    big_float a, b, t, a2(4 * g_precision), b2(4 * g_precision), t2(4 * g_precision);
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
                // The modes share the reference with more bits, and list_exact while the result stays the same
                bool wide = false, inexact = false, exact = false;
                std::string exact_r;
                for (int m = 0; m < 5; ++m)
                {
                    fesetround(g_modes[m]);
                    T r = fd(x, y);
                    fesetround(g_modes[0]);
                    ++n;
                    if (boundary_state(r, t.v) == 0)
                    {
                        check(r, t.v, m, name, as, bad);
                        continue;
                    }
                    // A reference on a boundary is exact, or it has too few bits to give the side: take more bits
                    g_precision *= 4;
                    if (!wide)
                    {
                        mpfr_clear_flags();
                        to_mpfr(x, a2.v);
                        to_mpfr(y, b2.v);
                        fm(t2.v, a2.v, b2.v);
                        // A reference with no rounding is exact, so a boundary there is a true one
                        inexact = mpfr_inexflag_p();
                        wide = true;
                    }
                    g_undecided = false;
                    if (inexact && boundary_state(r, t2.v) != 0)
                    {
                        if (exact_r != show(r))
                        {
                            exact_r = show(r);
                            exact = list_exact(name, x, y, r, t2.v);
                        }
                        g_undecided = !exact;
                    }
                    check(r, t2.v, m, name, as, bad);
                    g_precision /= 4;
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
                    check_ref(r, res.v, m, name, as, bad, [&](mpfr_t o) { big_float a2; to_mpfr(x, a2.v); poly_value(kind, l, 0, a2.v, o); },
                              [&](mpq_t o) { return poly_exact_q(kind, l, 0, x, o); });
                }
            }
    printf("%s n=%ld bad=%ld\n", name, n, bad);
    fflush(stdout);
    g_bad_total += bad;
}
// Each k-th value of a grid, for the short run
template <class T> std::vector<T> thin(const std::vector<T> &v, size_t k)
{
    std::vector<T> o;
    for (size_t i = 0; i < v.size(); i += k)
        o.push_back(v[i]);
    return o;
}
template <class T> void check_main(const char *tn, size_t k = 1)
{
    g_precision = 700;
    auto G = thin(wide_grid<T>(false), k);
    auto H = thin(wide_grid<T>(true), k);
    printf("== %s G=%zu H=%zu\n", tn, G.size(), H.size());
    std::vector<T> H2;
    for (size_t i = 0; i < H.size(); i += 4)
        H2.push_back(H[i]);
    for (int m : {1, 2, 3, 5, 7})
        for (int e = -3; e <= 2; e += static_cast<int>(k))
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
                        check_ref(r, t.v, m, "pow_int", as, bad, [&](mpfr_t o) { big_float a2; to_mpfr(x, a2.v); mpfr_pow_si(o, a2.v, k, MPFR_RNDN); },
                                  [&](mpq_t o) { return pow_exact_q(x, k, o); });
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
                    check_ref(r, t.v, m, "hypot3_x", as, bad, [&](mpfr_t o) {
                        big_float a2, b2, c2, s2;
                        to_mpfr(x, a2.v);
                        to_mpfr(y, b2.v);
                        to_mpfr(z, c2.v);
                        mpfr_sqr(o, a2.v, MPFR_RNDN);
                        mpfr_sqr(s2.v, b2.v, MPFR_RNDN);
                        mpfr_add(o, o, s2.v, MPFR_RNDN);
                        mpfr_sqr(s2.v, c2.v, MPFR_RNDN);
                        mpfr_add(o, o, s2.v, MPFR_RNDN);
                        mpfr_sqrt(o, o, MPFR_RNDN);
                    }, [&](mpq_t o) { return hypot3_exact_q(x, y, z, o); });
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
            check_ref(res, y.v, m, name, as, bad, [&](mpfr_t o) { big_float x2; to_mpfr(q, x2.v); poly_value(kind, n, mm, x2.v, o); },
                      [&](mpq_t o) { return poly_exact_q(kind, n, mm, q, o); });
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
    // ln(gamma(m)) goes into the exponent: gamma(m) alone and the other factor alone go out of the range of MPFR
    mpfr_lngamma(t, m, MPFR_RNDN);
    mpfr_sub(t, t, sum, MPFR_RNDN);
    mpfr_exp(t, t, MPFR_RNDN);
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
    // A double holds each n exactly, and a long of 32 bits does not hold the last one
    for (double n : {1.0,  2.0,  3.0,  4.0,   5.0,   6.0,    7.0,    8.0,      9.0,       10.0,  11.0,  12.0, 13.0,
                     16.0, 31.0, 32.0, 100.0, 101.0, 1000.0, 1001.0, 100000.0, 1000001.0, 1e9,   1e9 + 1, 1e12})
    {
        mpfr_mul_d(t.v, pi.v, n, MPFR_RNDN);
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
                    check_ref(f, oF.v, m, "ellint_1", as, bF,
                              [&](mpfr_t o) { big_float a2, b2, e2; to_mpfr(k, a2.v); to_mpfr(p, b2.v); ellint_values(a2.v, b2.v, o, e2.v); },
                              [&](mpq_t o) { return k == T{0} && (decimal_q(p, o), true); });
                    check_ref(e, oE.v, m, "ellint_2", as, bE,
                              [&](mpfr_t o) { big_float a2, b2, f2; to_mpfr(k, a2.v); to_mpfr(p, b2.v); ellint_values(a2.v, b2.v, f2.v, o); },
                              [&](mpq_t o) { return k == T{0} && (decimal_q(p, o), true); });
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
                check_ref(f, oF.v, m, "comp_ellint_1", as, bC, [&](mpfr_t o) {
                    big_float a2, e2; to_mpfr(k, a2.v); comp_ellint_values(a2.v, o, e2.v); });
                check_ref(e, oE.v, m, "comp_ellint_2", as, bC, [&](mpfr_t o) {
                    big_float a2, f2; to_mpfr(k, a2.v); comp_ellint_values(a2.v, f2.v, o); });
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
        {"erfc", mpfr_erfc},   {"tgamma", mpfr_gamma}, {"zeta", mpfr_zeta},   {"sqrt", mpfr_sqrt},
        {"cbrt", mpfr_cbrt}};
    for (const auto &e : one)
    {
        if (f == e.first)
        {
            const int i = e.second(o, a, MPFR_RNDN);
            ternary = (i > 0) - (i < 0);
            return true;
        }
    }
    if (f == "lgamma")
    {
        int sign;
        const int i = mpfr_lgamma(o, &sign, a, MPFR_RNDN);
        ternary = (i > 0) - (i < 0);
        return true;
    }
    if (f == "beta" || f == "beta_large")
    {
        // The library gives a NaN for an argument that is not above 0
        if (mpfr_nan_p(a) || mpfr_nan_p(b) || mpfr_sgn(a) <= 0 || mpfr_sgn(b) <= 0)
            mpfr_set_nan(o);
        else if (mpfr_inf_p(a) || mpfr_inf_p(b))
            mpfr_set_zero(o, 1);
        else if (mpfr_cmp_d(a, 1e15) >= 0 || mpfr_cmp_d(b, 1e15) >= 0)
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
    if (f == "atan2" || f == "hypot")
    {
        const int i = (f == "atan2" ? mpfr_atan2 : mpfr_hypot)(o, a, b, MPFR_RNDN);
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
// The result of the library in r, or false for a name that has no function here
template <class T> bool list_result(const std::string &f, T a, T b, T &r)
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
        {"lgamma", [](T x) -> T { return lgamma(x); }}, {"zeta", [](T x) -> T { return riemann_zeta(x); }},
        {"sqrt", [](T x) -> T { return sqrt(x); }},     {"cbrt", [](T x) -> T { return cbrt(x); }}};
    using G = T (*)(T, T);
    static const std::pair<const char *, G> two[] = {
        {"pow", [](T x, T y) -> T { return pow(x, y); }},           {"atan2", [](T x, T y) -> T { return atan2(x, y); }},
        {"hypot", [](T x, T y) -> T { return hypot(x, y); }},       {"beta", [](T x, T y) -> T { return beta(x, y); }},
        {"beta_large", [](T x, T y) -> T { return beta(x, y); }},   {"ellint_1", [](T x, T y) -> T { return ellint_1(x, y); }},
        {"ellint_2", [](T x, T y) -> T { return ellint_2(x, y); }}};
    for (const auto &e : one)
    {
        if (f == e.first)
        {
            r = e.second(a);
            return true;
        }
    }
    for (const auto &e : two)
    {
        if (f == e.first)
        {
            r = e.second(a, b);
            return true;
        }
    }
    return false;
}
// The reference t against the digits of r: 0 away from a rounding boundary, 1 within the precision of one, and 2 on
// one. The boundaries are the multiples of a twentieth of the unit of r, which takes in those below a power of ten.
template <class T> int boundary_state(T r, mpfr_t t)
{
    if (isnan(r) || !mpfr_regular_p(t))
        return 0;
    // An infinite r is next to the overflow threshold, a boundary of the largest value, and not decided by the result
    if (isinf(r))
        r = signbit(r) ? -std::numeric_limits<T>::max() : std::numeric_limits<T>::max();
    static big_float d(check_precision), u(check_precision);
    // Above two times the largest value only the overflow rules of check apply
    to_mpfr(std::numeric_limits<T>::max(), d.v);
    mpfr_mul_ui(d.v, d.v, 2, MPFR_RNDN);
    if (mpfr_cmpabs(t, d.v) > 0)
        return 0;
    const int e = unit_exponent(r);
    static big_float p10;
    pow10_exact(p10.v, e);
    // 0 is no boundary: a regular t far below the unit of r has the sign of the exact value
    for (const bool difference : {false, true})
    {
        mpfr_set(d.v, t, MPFR_RNDN);
        if (difference)
        {
            to_mpfr(r, d.v);
            mpfr_sub(d.v, d.v, t, MPFR_RNDN);
        }
        if (e < 0)
            mpfr_mul(d.v, d.v, p10.v, MPFR_RNDN);
        else
            mpfr_div(d.v, d.v, p10.v, MPFR_RNDN);
        mpfr_mul_ui(d.v, d.v, 20, MPFR_RNDN);
        if (!difference && mpfr_cmpabs_ui(d.v, 1) < 0 && mpfr_cmp_ui(d.v, 0) != 0)
            return 0;
    }
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
    static big_float u(check_precision), s;
    const int e = unit_exponent(r);
    pow10_exact(s.v, e);
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
// For arguments with a binary error, the sign of list_value where the four corners of their binary intervals round to
// t from one side, for a function that is monotonic in each argument there. Else 2.
template <class T> int list_side(const std::string &f, mpfr_t t, T x, T y)
{
    big_float a, b, o;
    int side = 2;
    for (int corner = 0; corner < 4; ++corner)
    {
        to_mpfr(x, a.v, corner & 1 ? MPFR_RNDU : MPFR_RNDD);
        to_mpfr(y, b.v, corner & 2 ? MPFR_RNDU : MPFR_RNDD);
        int s = 2;
        list_value(f, o.v, a.v, b.v, s);
        if (!mpfr_equal_p(o.v, t) || s == 0 || s == 2 || (corner != 0 && s != side))
            return 2;
        side = s;
    }
    return side;
}
// A 10^a == B 10^b for the integers A and B
static bool scaled_equal(mpz_srcptr A, long a, mpz_srcptr B, long b)
{
    if (a - b > 100000 || b - a > 100000)
        return false;
    mpz_t p;
    mpz_init(p);
    mpz_ui_pow_ui(p, 10, static_cast<unsigned long>(a > b ? a - b : b - a));
    mpz_mul(p, p, a > b ? A : B);
    const bool same = mpz_cmp(p, a > b ? B : A) == 0;
    mpz_clear(p);
    return same;
}
// o 10^e, for e >= 0
static void scale10(mpz_ptr o, long e)
{
    mpz_t p;
    mpz_init(p);
    mpz_ui_pow_ui(p, 10, static_cast<unsigned long>(e));
    mpz_mul(o, o, p);
    mpz_clear(p);
}
// d 10^e as an integer, for e >= 0
static void decimal_integer(mpz_ptr o, const std::string &d, long e)
{
    mpz_set_str(o, d.empty() ? "0" : d.c_str(), 10);
    scale10(o, e);
}
// If the value of f at the decimal arguments is the boundary of r next to t, sets t to that boundary at the full
// precision. f is log10, sqrt, cbrt, hypot, pow with x = 1 or y = P/Q for |P| up to 10000 and Q up to 1000, or beta
// with an argument in 1 to 64. A name can have a suffix after _.
template <class T> bool list_exact(const std::string &name, T x, T y, T r, mpfr_t t)
{
    const std::string f = name.substr(0, name.find('_'));
    // hypot is symmetric: a zero argument goes second, where the test below takes it
    if (f == "hypot" && x == T{0})
        std::swap(x, y);
    std::string xd, yd;
    long xe = 0, ye = 0;
    if ((f != "log10" && f != "pow" && f != "sqrt" && f != "cbrt" && f != "hypot" && f != "beta") ||
        !field_digits(show(x), xd, xe) || xd.empty() || !field_digits(show(y), yd, ye))
        return false;
    // The boundary is c 10^g: a multiple of a twentieth of the unit of r, as in boundary_state
    static big_float u(check_precision), p10;
    const int e = unit_exponent(r);
    pow10_exact(p10.v, e);
    if (e < 0)
        mpfr_mul(u.v, t, p10.v, MPFR_RNDN);
    else
        mpfr_div(u.v, t, p10.v, MPFR_RNDN);
    mpfr_mul_ui(u.v, u.v, 20, MPFR_RNDN);
    mpz_t c, m, A, B;
    mpz_inits(c, m, A, B, nullptr);
    mpfr_get_z(c, u.v, MPFR_RNDN);
    mpz_mul_ui(c, c, 5);
    const long g = e - 2;
    mpz_set_str(m, xd.c_str(), 10);
    // d 10^k as an integer below 100, or -1
    auto small = [](const std::string &d, long k)
    {
        const bool fits = !d.empty() && k >= 0 && d.size() + static_cast<std::size_t>(k) <= 2;
        return fits ? std::stol(d) * (k > 0 ? 10 : 1) : -1L;
    };
    bool same = false;
    if (f == "log10")
    {
        // x = 10^xe only if its digits are 1
        mpz_set_si(B, xe);
        same = !signbit(x) && xd == "1" && scaled_equal(c, g, B, 0);
    }
    else if (f == "hypot")
    {
        // c^2 10^(2 g) == (mx^2 10^sx + my^2 10^sy) 10^e0
        const long e0 = yd.empty() || xe < ye ? 2 * xe : 2 * ye;
        const long sx = 2 * xe - e0, sy = yd.empty() ? 0 : 2 * ye - e0;
        if (mpz_sgn(c) > 0 && sx <= 100000 && sy <= 100000)
        {
            mpz_mul(A, c, c);
            mpz_mul(B, m, m);
            scale10(B, sx);
            decimal_integer(m, yd, 0);
            mpz_mul(m, m, m);
            scale10(m, sy);
            mpz_add(B, B, m);
            same = scaled_equal(A, 2 * g, B, e0);
        }
    }
    else if (f == "beta")
    {
        // For an integer n: beta(x, n) (x)(x + 1)...(x + n - 1) == (n - 1)!, with x = X / 10^s
        const long yn = small(yd, ye);
        const bool swap = !(yn >= 1 && yn <= 64);
        const long n = swap ? small(xd, xe) : yn;
        const std::string &od = swap ? yd : xd;
        const long oe = swap ? ye : xe;
        if (n >= 1 && n <= 64 && !od.empty() && !signbit(x) && !signbit(y) && mpz_sgn(c) > 0)
        {
            const long sc = oe < 0 ? -oe : 0;
            decimal_integer(m, od, oe + sc);
            mpz_set(A, c);
            mpz_t d, p;
            mpz_inits(d, p, nullptr);
            mpz_ui_pow_ui(d, 10, static_cast<unsigned long>(sc));
            for (long i = 0; i < n; ++i)
            {
                mpz_mul_si(p, d, i);
                mpz_add(p, p, m);
                mpz_mul(A, A, p);
            }
            mpz_fac_ui(B, static_cast<unsigned long>(n - 1));
            mpz_clears(d, p, nullptr);
            same = scaled_equal(A, g, B, sc * n);
        }
    }
    else if (f == "pow" && xd == "1" && xe == 0 && !signbit(x))
    {
        // pow(1, y) is 1 for each y
        mpz_set_ui(B, 1);
        same = scaled_equal(c, g, B, 0);
    }
    else if (f != "pow" || (yd.size() <= 4 && ye <= 3 && ye >= -3))
    {
        long P = f == "pow" ? std::stol(yd.empty() ? "0" : yd) : 1, Q = f == "sqrt" ? 2 : f == "cbrt" ? 3 : 1;
        for (long k = f == "pow" ? ye : 0; k > 0; --k)
            P *= 10;
        for (long k = f == "pow" ? ye : 0; k < 0; ++k)
            Q *= 10;
        P = f == "pow" && signbit(y) ? -P : P;
        for (const long k : {2L, 5L})
            for (; P % k == 0 && Q % k == 0; P /= k, Q /= k)
                ;
        // c^Q 10^(g Q) == m^P 10^(xe P), with m^|P| on the left for P < 0, and the sign of x^P for x < 0
        const int sign = signbit(x) && P % 2 != 0 ? -1 : 1;
        if (P >= -10000 && P <= 10000 && Q <= 1000 && (Q % 2 != 0 || !signbit(x)) && mpz_sgn(c) == sign)
        {
            mpz_abs(A, c);
            mpz_pow_ui(A, A, static_cast<unsigned long>(Q));
            mpz_pow_ui(B, m, static_cast<unsigned long>(P < 0 ? -P : P));
            if (P < 0)
            {
                mpz_mul(A, A, B);
                mpz_set_ui(B, 1);
            }
            same = scaled_equal(A, g * Q, B, xe * P);
        }
    }
    if (same)
    {
        mpfr_set_prec(t, check_precision);
        mpfr_set_z(t, c, MPFR_RNDN);
        pow10_exact(p10.v, static_cast<int>(g));
        if (g < 0)
            mpfr_div(t, t, p10.v, MPFR_RNDN);
        else
            mpfr_mul(t, t, p10.v, MPFR_RNDN);
    }
    mpz_clears(c, m, A, B, nullptr);
    return same;
}
template <class T> long check_lines(std::istream &in)
{
    std::string line;
    long count = 0;
    while (std::getline(in, line))
    {
        ++count;
        const auto p1 = line.find('|');
        const auto p2 = p1 == std::string::npos ? p1 : line.find('|', p1 + 1);
        const auto p3 = p2 == std::string::npos ? p2 : line.find('|', p2 + 1);
        T x{}, y{};
        // A line with a missing field, an argument that is not a decimal, or a mode not in 0 to 4 counts as wrong
        if (p3 == std::string::npos || !list_argument(line.substr(p1 + 1, p2 - p1 - 1), x) ||
            !list_argument(line.substr(p2 + 1, p3 - p2 - 1), y) || line.size() != p3 + 2 || line[p3 + 1] < '0' ||
            line[p3 + 1] > '4')
        {
            printf("%s|invalid line %ld\n", line.c_str(), count);
            ++g_bad_total;
            continue;
        }
        const std::string f = line.substr(0, p1);
        const int m = line[p3 + 1] - '0';
        fesetround(g_modes[m]);
        T r{};
        const bool known = list_result<T>(f, x, y, r);
        fesetround(g_modes[0]);
        if (!known)
        {
            printf("%s|unknown\n", line.c_str());
            ++g_bad_total;
            continue;
        }
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
            mpfr_clear_flags();
            to_mpfr(x, a.v);
            to_mpfr(y, b.v);
            const bool exact = !mpfr_inexflag_p();
            mpfr_clear_flags();
            int ternary = 2;
            if (!list_value(f, t.v, a.v, b.v, ternary))
            {
                printf("%s|unknown\n", line.c_str());
                ++g_bad_total;
                break;
            }
            const bool overflow = mpfr_overflow_p(), underflow = mpfr_underflow_p();
            int state = boundary_state(r, t.v);
            // An exact decimal value on the boundary decides first, as the sign of MPFR is that of a binary rounding.
            // The sign is that of the binary arguments: for the decimal ones, it comes from both ends of each.
            if (state != 0 && list_exact(f, x, y, r, t.v))
            {
                state = 2;
                ternary = 0;
            }
            else if (!exact && ternary != 2)
                ternary = state == 2 && monotonic(f, x) ? list_side(f, t.v, x, y) : 2;
            if (state == 2 && ternary != 2)
            {
                if (ternary != 0)
                    move_off_boundary(r, t.v, ternary);
            }
            else if (state != 0)
            {
                if (g_precision * 2 <= check_precision)
                    continue;
                g_undecided = true;
            }
            g_overflow = overflow;
            g_underflow = underflow;
            long bad = 0;
            const bool undecided = g_undecided;
            check(r, t.v, m, f.c_str(), show(x) + "," + show(y), bad);
            const char *mark = bad != 0 ? "|wrong" : undecided ? "|inconclusive" : "";
            printf("%s|%s%s\n", line.c_str(), show(r).c_str(), mark);
            g_bad_total += bad;
            break;
        }
    }
    return count;
}
template <class T> void check_list(const char *path)
{
    std::ifstream in(path);
    if (check_lines<T>(in) == 0)
    {
        printf("%s: no line to check\n", path);
        ++g_bad_total;
    }
}
// Gives wrong results to check(), which must reject each of them, and the right results, which it must accept
template <class T> void check_self_cases()
{
    static big_float t(check_precision), v(check_precision), h(check_precision);
    const T dm = std::numeric_limits<T>::denorm_min(), mx = std::numeric_limits<T>::max(), one{1};
    const T inf = std::numeric_limits<T>::infinity(), qnan = std::numeric_limits<T>::quiet_NaN();
    int e;
    frexp10(one, &e);
    const T below = one - T{1, e - 1};
    long n = 0, missed = 0;
    // The value t and the result r with both signs: up and down change places for a negative value
    auto probe = [&](T r, int mode, bool right, const char *name)
    {
        for (int sg = 0; sg < 2; ++sg)
        {
            // From 8, check() does not print
            long bad = 8;
            mpfr_set(v.v, t.v, MPFR_RNDN);
            if (sg)
                mpfr_neg(v.v, v.v, MPFR_RNDN);
            const int m = sg && (mode == 1 || mode == 2) ? 3 - mode : mode;
            ++n;
            if (check<T>(sg ? -r : r, v.v, m, name, sg ? "-" : "+", bad) != right)
            {
                ++missed;
                printf("  self: %s mode %d sign %d has the wrong verdict\n", name, m, sg);
            }
        }
    };
    if (is_fast<T>::v)
    {
        // A fast type has no subnormal results: only the results past the range, which had no comparison
        to_mpfr(mx, t.v);
        mpfr_mul_ui(t.v, t.v, 2, MPFR_RNDN);
        probe(one, 0, false, "fast overflow");
        probe(mx, 0, false, "fast overflow");
        probe(inf, 0, true, "fast overflow");
        probe(mx, 4, true, "fast overflow");
        to_mpfr(std::numeric_limits<T>::min(), t.v);
        mpfr_div_ui(t.v, t.v, 1000, MPFR_RNDN);
        probe(one, 0, false, "fast underflow");
        probe(T{0}, 0, true, "fast underflow");
        printf("self n=%ld bad=%ld\n", n, missed);
        g_bad_total += missed;
        return;
    }
    // A value far below the smallest subnormal value
    to_mpfr(dm, t.v);
    mpfr_div_ui(t.v, t.v, 1000000, MPFR_RNDN);
    for (int mode = 0; mode < 5; ++mode)
    {
        probe(dm, mode, mode == 1, "tiny");
        probe(T{0}, mode, mode != 1, "tiny");
    }
    // The same value far below the precision of the comparison
    mpfr_div_2ui(t.v, t.v, 8000, MPFR_RNDN);
    probe(dm, 1, true, "deep tiny");
    probe(T{0}, 1, false, "deep tiny");
    probe(T{0}, 0, true, "deep tiny");
    // A quarter of the smallest subnormal value: a zero result has the sign of the value
    to_mpfr(dm, t.v);
    mpfr_div_ui(t.v, t.v, 4, MPFR_RNDN);
    for (int sg = 0; sg < 2; ++sg)
    {
        long bad = 8;
        mpfr_set(v.v, t.v, MPFR_RNDN);
        if (sg)
            mpfr_neg(v.v, v.v, MPFR_RNDN);
        n += 2;
        if (!check<T>(sg ? -T{0} : T{0}, v.v, 0, "zero sign", sg ? "-" : "+", bad) ||
            check<T>(sg ? T{0} : -T{0}, v.v, 0, "zero sign", sg ? "-" : "+", bad))
        {
            ++missed;
            printf("  self: zero sign %d has the wrong verdict\n", sg);
        }
    }
    // The two directed conversions of a negative argument enclose it
    to_mpfr(-T{1, -1}, v.v, MPFR_RNDD);
    to_mpfr(-T{1, -1}, h.v, MPFR_RNDU);
    mpfr_mul_ui(t.v, v.v, 10, MPFR_RNDD);
    ++n;
    if (mpfr_cmp(v.v, h.v) >= 0 || mpfr_cmp_si(t.v, -1) >= 0)
    {
        ++missed;
        printf("  self: the directed conversions of -0.1 do not enclose it\n");
    }
    // A tie between 0 and the smallest subnormal value
    to_mpfr(dm, t.v);
    mpfr_div_ui(t.v, t.v, 2, MPFR_RNDN);
    probe(T{0}, 0, true, "subnormal tie");
    probe(dm, 0, false, "subnormal tie");
    probe(dm, 3, true, "subnormal tie");
    probe(T{0}, 3, false, "subnormal tie");
    // Two times the largest value
    to_mpfr(mx, t.v);
    mpfr_mul_ui(t.v, t.v, 2, MPFR_RNDN);
    probe(mx, 0, false, "overflow");
    probe(mx, 3, false, "overflow");
    probe(inf, 0, true, "overflow");
    probe(mx, 2, true, "overflow");
    probe(mx, 4, true, "overflow");
    probe(inf, 1, true, "overflow");
    // A quarter of a unit above the largest value, and the overflow threshold at half a unit
    int em;
    frexp10(mx, &em);
    mpfr_set_ui(h.v, 10, MPFR_RNDN);
    mpfr_pow_si(h.v, h.v, em, MPFR_RNDN);
    mpfr_div_ui(h.v, h.v, 4, MPFR_RNDN);
    to_mpfr(mx, t.v);
    mpfr_add(t.v, t.v, h.v, MPFR_RNDN);
    probe(mx, 0, true, "below threshold");
    probe(inf, 0, false, "below threshold");
    mpfr_add(t.v, t.v, h.v, MPFR_RNDN);
    probe(mx, 0, false, "threshold");
    probe(inf, 0, true, "threshold");
    probe(inf, 3, true, "threshold");
    // Just below the threshold nearest gives the largest value. At 1500 bits the reference is on the threshold, and
    // an infinite result must not decide it, so that the precision goes up.
    mpfr_set_ui_2exp(h.v, 1, -1600, MPFR_RNDN);
    mpfr_sub(t.v, t.v, h.v, MPFR_RNDN);
    probe(mx, 0, true, "below threshold");
    probe(inf, 0, false, "below threshold");
    {
        big_float lo(1500);
        mpfr_set(lo.v, t.v, MPFR_RNDN);
        const int keep = g_precision;
        g_precision = 1500;
        for (const T r : {inf, mx})
        {
            ++n;
            if (decides(r, lo.v, 2))
            {
                ++missed;
                printf("  self: %s decides a reference on the overflow threshold\n", show(r).c_str());
            }
        }
        g_precision = keep;
    }
    // 0.6 of a unit below 1, where the unit is that of the values below 1
    to_mpfr(one - below, h.v);
    mpfr_mul_ui(h.v, h.v, 6, MPFR_RNDN);
    mpfr_div_ui(h.v, h.v, 10, MPFR_RNDN);
    to_mpfr(one, t.v);
    mpfr_sub(t.v, t.v, h.v, MPFR_RNDN);
    probe(one, 0, false, "power of ten");
    probe(below, 0, true, "power of ten");
    probe(one, 1, true, "power of ten");
    probe(below, 2, true, "power of ten");
    probe(one, 2, false, "power of ten");
    probe(below, 1, false, "power of ten");
    // An infinite zeta without an overflow is the pole, which gives NaN; with an overflow it rounds as other results
    mpfr_set_inf(t.v, 1);
    probe(qnan, 0, true, "zeta");
    probe(inf, 0, false, "zeta");
    g_overflow = true;
    probe(qnan, 0, false, "zeta");
    probe(inf, 0, true, "zeta");
    probe(mx, 2, true, "zeta");
    probe(inf, 2, false, "zeta");
    g_overflow = false;
    // A reference on a boundary with no sign of its error is not decided, and it counts as inconclusive
    const long inconclusive = g_inconclusive;
    g_say_undecided = false;
    check_unary<T>(
        "self undecided", std::vector<T>{one}, [](T x) { return x; },
        [](mpfr_t o, mpfr_t a)
        {
            mpfr_set(o, a, MPFR_RNDN);
            mpfr_set_inexflag();
        },
        500);
    n += 10;
    if (g_inconclusive != inconclusive + 10)
    {
        ++missed;
        printf("  self: an undecided reference is not inconclusive\n");
    }
    // On an undecided boundary, a result that is right on one side is not wrong, and one that is right on none is
    struct side_case
    {
        T r;
        int mode;
        bool right;
    };
    const T above = one + T{1, e};
    mpfr_set_ui(t.v, 1, MPFR_RNDN);
    for (const side_case &c : {side_case{above, 1, true}, side_case{below, 2, true}, side_case{one, 4, true},
                               side_case{above, 0, false}, side_case{above + T{1, e}, 1, false}})
    {
        long bad = 8;
        g_undecided = true;
        ++n;
        if (check<T>(c.r, t.v, c.mode, "self side", "+", bad) != c.right)
        {
            ++missed;
            printf("  self: an undecided reference has the wrong verdict for %s mode %d\n", show(c.r).c_str(), c.mode);
        }
    }
    g_inconclusive = inconclusive;
    g_say_undecided = true;
    // Exact results of arguments that have a binary error, in the directed modes, and a beta far below the range
    const long before = g_bad_total;
    std::istringstream exact_lines(
        "pow|1e-1|-1e0|1\npow|1e-1|-1e0|2\npow|1e-1|-1e0|4\npow|2e-1|-1e0|2\npow|1e-2|5e-1|1\npow|1e-2|5e-1|2\n"
        "log10|1e-3|0|1\nbeta|1e15|1e15|0\n");
    n += check_lines<T>(exact_lines);
    // Zeros after the last digit keep the value of an argument, and a digit that the type cannot hold is not valid
    const std::string zeros(200, '0');
    std::istringstream zero_lines("log10|1." + zeros + "e-3|0|1\n");
    n += check_lines<T>(zero_lines);
    if (g_inconclusive != inconclusive)
    {
        ++missed;
        printf("  self: an exact result is inconclusive\n");
    }
    // list_exact accepts a boundary only if it is the exact value
    const struct
    {
        const char *f;
        T x, y, r;
        bool exact;
    } cases[] = {{"pow", T{1, -2}, T{5, -1}, T{1, -1}, true},  {"pow", T{2, -1}, T{-2}, T{25}, true},
                 {"pow", T{1, -1}, T{3}, T{1, -3}, true},     {"pow", T{2, -2}, T{5, -1}, T{1, -1}, false},
                 {"log10", T{1, -3}, T{0}, T{-3}, true},      {"log10", T{2, -3}, T{0}, T{-3}, false},
                 {"pow_x", T{-2, -1}, T{-3}, T{-125}, true},  {"pow_x", T{-2, -1}, T{-2}, T{25}, true},
                 {"pow", T{-1, -2}, T{-5, -1}, T{-10}, false}, {"log10", T{-1, -3}, T{0}, T{-3}, false},
                 {"sqrt", T{1, -2}, T{0}, T{1, -1}, true},    {"sqrt", T{2, -2}, T{0}, T{1, -1}, false},
                 {"cbrt", T{-8, -3}, T{0}, T{-2, -1}, true},  {"pow", T{1, -1}, T{70}, T{1, -70}, true},
                 {"hypot_x", T{3, -2}, T{4, -2}, T{5, -2}, true}, {"hypot", T{3, -2}, T{5, -2}, T{6, -2}, false},
                 {"hypot_x", T{11, -1}, T{6}, T{61, -1}, true}, {"hypot", T{0}, T{1, -50}, T{1, -50}, true},
                 {"beta_n", T{1, -5}, T{1}, T{1, 5}, true},   {"beta", T{2}, T{5, -1}, T{4, -1}, false},
                 {"beta", T{3}, T{2}, T{8333333, -8}, false}, {"beta", T{1, -1}, T{3}, T{50, -1}, false},
                 {"beta_large", T{1}, T{1, 15}, T{1, -15}, true}, {"pow", T{1, -1}, T{12345, -1}, T{1}, false},
                 {"pow", T{1}, T{11, -11}, T{1}, true},        {"pow", T{1}, T{-75, -10}, T{1}, true}};
    for (const auto &k : cases)
    {
        big_float a, b, o;
        to_mpfr(k.x, a.v);
        to_mpfr(k.y, b.v);
        int side = 2;
        // MPFR has no list_value for the names with a suffix, or for sqrt, cbrt and hypot
        const std::string base = std::string(k.f).substr(0, std::string(k.f).find('_'));
        if (base == "sqrt")
            mpfr_sqrt(o.v, a.v, MPFR_RNDN);
        else if (base == "cbrt")
            mpfr_cbrt(o.v, a.v, MPFR_RNDN);
        else if (base == "hypot")
            mpfr_hypot(o.v, a.v, b.v, MPFR_RNDN);
        else
            list_value(base, o.v, a.v, b.v, side);
        ++n;
        if (list_exact(k.f, k.x, k.y, k.r, o.v) != k.exact)
        {
            ++missed;
            printf("  self: list_exact(%s, %s, %s) is not %d\n", k.f, show(k.x).c_str(), show(k.y).c_str(), k.exact);
        }
    }
    const long valid = g_bad_total;
    std::istringstream tail_lines("log|1." + zeros + "1|0|0\nlog|1e-9223372036854775808|0|0\nlog|1e|0|0\n"
                                  "exp|--inf|0e0|0\nexp|--nan|0e0|0\n");
    n += check_lines<T>(tail_lines);
    if (g_bad_total != valid + 5)
    {
        ++missed;
        printf("  self: an argument that is not valid passes\n");
    }
    g_bad_total = valid;
    printf("self n=%ld bad=%ld\n", n, missed + g_bad_total - before);
    g_bad_total += missed;
}
// The results of the self test are no records for a replay
template <class T> void check_self()
{
    FILE *const rec = g_rec;
    g_rec = nullptr;
    check_self_cases<T>();
    g_rec = rec;
}
// The part short takes each k-th value of the grids of main
static constexpr size_t short_stride = 16;
int main(int argc, char **argv)
{
    if (argc < 3 || argc != (std::string(argv[2]) == "list" ? 4 : 3))
    {
        std::puts("usage: rounding_grid d32|d64|d128|f32|f64|f128 "
                  "main|short|exact|near|ext|ell|assoc|self|list <file>");
        return 2;
    }
    std::string t = argv[1], w = argv[2];
    if (std::getenv("GRID_REC") && !g_rec)
    {
        std::printf("cannot open the file of GRID_REC: %s\n", std::getenv("GRID_REC"));
        return 2;
    }
    const bool known_type = t == "d32" || t == "d64" || t == "d128" || t == "f32" || t == "f64" || t == "f128";
    const bool known_part =
        w == "main" || w == "short" || w == "exact" || w == "near" || w == "ext" || w == "ell" || w == "assoc" ||
        w == "self" || w == "list";
    if (!known_type || !known_part)
    {
        std::puts("usage: rounding_grid d32|d64|d128|f32|f64|f128 "
                  "main|short|exact|near|ext|ell|assoc|self|list <file>");
        return 2;
    }
#ifdef BOUND_MARGIN
    margin_self_test();
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
        if (w == "short")                                                                                              \
            check_main<T>(n, short_stride);                                                                            \
        if (w == "exact")                                                                                              \
            check_exact<T>();                                                                                          \
        if (w == "near")                                                                                               \
            check_near<T>();                                                                                           \
        if (w == "self")                                                                                               \
            check_self<T>();                                                                                           \
        if (w == "list")                                                                                               \
            check_list<T>(argv[3]);                                                                                    \
    }
    ALL(decimal32_t, "d32")
    ALL(decimal64_t, "d64") ALL(decimal128_t, "d128") ALL(decimal_fast32_t, "f32") ALL(decimal_fast64_t, "f64")
        ALL(decimal_fast128_t, "f128") printf("%s %s total bad %ld inconclusive %ld skipped %ld\n", t.c_str(),
                                              w.c_str(), g_bad_total, g_inconclusive, g_skipped_total);
    if (g_rec && std::fclose(g_rec) != 0)
        g_rec_failed = true;
    if (g_rec_failed)
        std::puts("the file of GRID_REC did not get all records");
    // With GRID_STRICT set, an inconclusive result also gives the status 1
    const bool strict_failed = std::getenv("GRID_STRICT") && g_inconclusive != 0;
    const bool margin_failed = !margin_report();
    return g_bad_total == 0 && !g_rec_failed && !strict_failed && !margin_failed ? 0 : 1;
}
