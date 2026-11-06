// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The checks of the error bounds in the headers that bound_mutate.py writes.
// See README.md in this directory.

#ifndef BOOST_DECIMAL_TOOLS_BOUND_RUNTIME_HPP
#define BOOST_DECIMAL_TOOLS_BOUND_RUNTIME_HPP

#include <cstdint>

namespace boost_mut {

struct site_count
{
    unsigned long long calls;
    unsigned long long nears;
};

// A check of a site with its digits and its bound, for bound_margin.hpp
struct bound_event
{
    int site;
    std::uint64_t digits_high;
    std::uint64_t digits_low;
    std::uint64_t bound_high;
    std::uint64_t bound_low;
};

constexpr int max_sites {512};
constexpr int max_hits {32};
constexpr int max_events {64};

inline int site {0};
inline int target {-1};
inline unsigned shift {0};
inline site_count counts[max_sites] {};
inline int hits[max_hits];
inline int hit_count {0};
inline bound_event events[max_events];
inline int event_count {0};
inline bool record_events {false};

constexpr auto set(const int id) noexcept -> bool
{
    if (!__builtin_is_constant_evaluated())
    {
        site = id;
    }
    return true;
}

// The bound of the selected site gets smaller by 2^shift. A shift of 99 or more removes the check.
template <typename W>
constexpr auto adjust(const W w) noexcept -> W
{
    if (__builtin_is_constant_evaluated() || site != target)
    {
        return w;
    }
    return shift >= 99U ? W {0U} : static_cast<W>(w >> shift);
}

template <typename U>
inline auto high_word(const U& u) noexcept -> std::uint64_t
{
    if constexpr (sizeof(U) > 8U)
    {
        return u.high;
    }
    else
    {
        return 0U;
    }
}

template <typename U>
inline auto low_word(const U& u) noexcept -> std::uint64_t
{
    if constexpr (sizeof(U) > 8U)
    {
        return u.low;
    }
    else
    {
        return static_cast<std::uint64_t>(u);
    }
}

template <typename D, typename W>
inline void note(const bool is_near, const D& d, const W& w) noexcept
{
    ++counts[site].calls;
    if (is_near)
    {
        ++counts[site].nears;
        if (hit_count < max_hits)
        {
            hits[hit_count++] = site;
        }
    }
    if (record_events && event_count < max_events)
    {
        events[event_count++] = bound_event {site, high_word(d), low_word(d), high_word(w), low_word(w)};
    }
}

template <typename D, typename R, typename W, typename M>
constexpr auto near(const D& d, const R& r, W w, const M& m) noexcept -> bool
{
    if (!__builtin_is_constant_evaluated())
    {
        note(r < w || r > m - w, d, w);
    }
    w = adjust(w);
    return w != W {0U} && (r < w || r > m - w);
}

template <typename D, typename R, typename U, typename W>
constexpr auto near_trig(const D& d, const R& r, const U& rem, W w, const U& unit) noexcept -> bool
{
    if (!__builtin_is_constant_evaluated())
    {
        note((r == 0U && rem < w) || (r == 499U && rem >= unit - w), d, w);
    }
    w = adjust(w);
    return w != W {0U} && ((r == 0U && rem < w) || (r == 499U && rem >= unit - w));
}

template <typename D, typename U>
constexpr auto poly_ambiguous(const D& d, U units, const U& h, const U& r) noexcept -> bool
{
    if (!__builtin_is_constant_evaluated())
    {
        note(units >= h || r <= units || r + units + 1U >= h, d, units);
    }
    units = adjust(units);
    return units >= h || r <= units || r + units + 1U >= h;
}

} // namespace boost_mut

#endif
