// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// With BOUND_MARGIN, rounding_grid measures how much of each error bound the digits use.
// See README.md in this directory.

#ifndef BOOST_DECIMAL_TOOLS_BOUND_MARGIN_HPP
#define BOOST_DECIMAL_TOOLS_BOUND_MARGIN_HPP

#include <bound_runtime.hpp>
#include <mpfr.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

struct site_margin
{
    unsigned long checks = 0;
    unsigned long skipped = 0;
    double largest = 0;
    std::string function;
    std::string arguments;
};

static site_margin g_margins[boost_mut::max_sites];

// GMP takes all the bits of each word, as unsigned long can have 32 bits
static void words_to_mpfr(const std::uint64_t high, const std::uint64_t low, mpfr_t o)
{
    const std::uint64_t w[2] = {low, high};
    mpz_t z;
    mpz_init(z);
    mpz_import(z, 2, -1, sizeof(std::uint64_t), 0, 0, w);
    mpfr_set_z(o, z, MPFR_RNDN);
    mpz_clear(z);
}

// |digits - value| / bound for each check since the last call, where the digits scale the value by a power of ten
static void margin_events(mpfr_t value, const char *name, const std::string &arg)
{
    static mpfr_t d, w, s, e;
    static bool ready = false;
    if (!ready)
    {
        mpfr_inits2(2000, d, w, s, e, static_cast<mpfr_ptr>(nullptr));
        ready = true;
        // Both words with bits above bit 31
        words_to_mpfr(UINT64_C(0x8000000100000003), UINT64_C(0x8000000200000005), d);
        mpfr_set_str(w, "80000001000000038000000200000005", 16, MPFR_RNDN);
        if (!mpfr_equal_p(d, w))
        {
            std::printf("words_to_mpfr loses bits\n");
            std::exit(2);
        }
    }
    for (int i = 0; i < boost_mut::event_count; ++i)
    {
        const auto &ev = boost_mut::events[i];
        auto &g = g_margins[ev.site];
        if (mpfr_nan_p(value) || mpfr_inf_p(value) || mpfr_zero_p(value))
            continue;
        if ((ev.digits_high == 0U && ev.digits_low == 0U) || (ev.bound_high == 0U && ev.bound_low == 0U))
            continue;
        words_to_mpfr(ev.digits_high, ev.digits_low, d);
        words_to_mpfr(ev.bound_high, ev.bound_low, w);
        mpfr_abs(s, value, MPFR_RNDN);
        mpfr_div(e, d, s, MPFR_RNDN);
        mpfr_log10(e, e, MPFR_RNDN);
        mpfr_round(e, e);
        mpfr_ui_pow(e, 10, e, MPFR_RNDN);
        mpfr_mul(s, s, e, MPFR_RNDN);
        mpfr_sub(e, d, s, MPFR_RNDN);
        mpfr_abs(e, e, MPFR_RNDN);
        mpfr_div(s, e, d, MPFR_RNDN);
        ++g.checks;
        // A check of another value, such as a part of the function, is not a measure of this bound
        if (mpfr_cmp_d(s, 1e-5) > 0)
        {
            ++g.skipped;
            continue;
        }
        mpfr_div(e, e, w, MPFR_RNDN);
        const double ratio = mpfr_get_d(e, MPFR_RNDN);
        if (ratio > g.largest)
        {
            g.largest = ratio;
            g.function = name;
            g.arguments = arg;
        }
    }
    boost_mut::event_count = 0;
}

// Feeds margin_events two checks of site 0: digits 10^18 + 3 with a bound of 10 against the value 1, which gives the
// ratio 0.3, and digits 10^18 with a bound of 1 against 1.00002, an error of 2 * 10^13 units that it must skip
static void margin_self_test()
{
    mpfr_t v;
    mpfr_init2(v, 200);
    const std::uint64_t e18 {UINT64_C(1000000000000000000)};
    boost_mut::events[0] = boost_mut::bound_event {0, 0U, e18 + 3U, 0U, 10U};
    boost_mut::event_count = 1;
    mpfr_set_ui(v, 1U, MPFR_RNDN);
    margin_events(v, "self", "+");
    boost_mut::events[0] = boost_mut::bound_event {0, 0U, e18, 0U, 1U};
    boost_mut::event_count = 1;
    mpfr_set_str(v, "1.00002", 10, MPFR_RNDN);
    margin_events(v, "self", "+");
    mpfr_clear(v);
    const auto &g = g_margins[0];
    if (g.checks != 2U || g.skipped != 1U || g.largest < 0.2999 || g.largest > 0.3001)
    {
        std::printf("margin_events failed its self-test\n");
        std::exit(2);
    }
    g_margins[0] = site_margin {};
}

// With MARGIN_OUT set to a file, the file gets a line for each site: site checks skipped ratio function arguments. The
// ratio is "unmeasured" where the program skipped each check. false, with a message, if the file does not get all lines.
static bool margin_report()
{
    const char *path = std::getenv("MARGIN_OUT");
    if (!path)
        return true;
    FILE *f = std::fopen(path, "w");
    bool ok = f != nullptr;
    for (int i = 0; ok && i < boost_mut::max_sites; ++i)
    {
        const auto &g = g_margins[i];
        if (g.checks == 0U)
            continue;
        if (g.checks == g.skipped)
            ok = std::fprintf(f, "%d %lu %lu unmeasured\n", i, g.checks, g.skipped) >= 0;
        else
            ok = std::fprintf(f, "%d %lu %lu %.4f %s %s\n", i, g.checks, g.skipped, g.largest, g.function.c_str(),
                              g.arguments.c_str()) >= 0;
    }
    ok = f != nullptr && std::fclose(f) == 0 && ok;
    if (!ok)
        std::printf("the file of MARGIN_OUT did not get all lines: %s\n", path);
    return ok;
}

#endif
