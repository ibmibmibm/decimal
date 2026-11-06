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

static void words_to_mpfr(const std::uint64_t high, const std::uint64_t low, mpfr_t o)
{
    mpfr_set_ui_2exp(o, high, 64, MPFR_RNDN);
    mpfr_add_ui(o, o, low, MPFR_RNDN);
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

// With MARGIN_OUT set to a file, the file gets a line for each site: site checks skipped ratio function arguments
struct margin_report
{
    ~margin_report()
    {
        const char *path = std::getenv("MARGIN_OUT");
        FILE *f = path ? std::fopen(path, "w") : nullptr;
        if (!f)
            return;
        for (int i = 0; i < boost_mut::max_sites; ++i)
        {
            const auto &g = g_margins[i];
            if (g.checks)
                std::fprintf(f, "%d %lu %lu %.4f %s %s\n", i, g.checks, g.skipped, g.largest, g.function.c_str(),
                             g.arguments.c_str());
        }
        std::fclose(f);
    }
};
static margin_report g_margin_report;

#endif
