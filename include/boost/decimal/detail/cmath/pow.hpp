// Copyright 2023 - 2026 Matt Borland
// Copyright 2023 - 2026 Christopher Kormanyos
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_DECIMAL_DETAIL_CMATH_POW_HPP
#define BOOST_DECIMAL_DETAIL_CMATH_POW_HPP

#include <boost/decimal/fwd.hpp> // NOLINT(llvm-include-order)
#include <boost/decimal/detail/type_traits.hpp>
#include <boost/decimal/detail/concepts.hpp>
#include <boost/decimal/detail/config.hpp>
#include <boost/decimal/detail/remove_trailing_zeros.hpp>
#include <boost/decimal/detail/cmath/log.hpp>
#include <boost/decimal/detail/cmath/impl/exp_impl.hpp>

#ifndef BOOST_DECIMAL_BUILD_MODULE
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#endif

namespace boost {
namespace decimal {

namespace detail {

namespace pow_detail {

// s^n in r, or false where it does not fit in 128 bits
constexpr auto pow_power(const boost::int128::uint128_t& s, unsigned n, boost::int128::uint128_t& r) noexcept -> bool
{
    r = 1U;
    const auto limit {(std::numeric_limits<boost::int128::uint128_t>::max)() / s};
    while (n-- > 0U)
    {
        if (r > limit)
        {
            return false;
        }
        r *= s;
    }
    return true;
}

// An exact power sig * 10^exp, where ok
struct pow_exact_value
{
    boost::int128::uint128_t sig;
    int exp;
    bool ok;
};

// x^n or x^-n for x = s * 10^ex with s not a multiple of 10, where it is exact in 128 bits. s^-n is a finite decimal
// only for s = 2^i or 5^j, as 5^(in) * 10^(-in) or 2^(jn) * 10^(-jn), and s^n of more than 38 digits is not exact.
constexpr auto pow_exact(boost::int128::uint128_t s, const int ex, unsigned n, const bool neg_n) noexcept -> pow_exact_value
{
    int scale {ex * static_cast<int>(n)};
    if (neg_n)
    {
        scale = -scale;
        const int i {boost::int128::countr_zero(s)};
        unsigned j {};
        if ((s >> static_cast<unsigned>(i)) == 1U)
        {
            j = static_cast<unsigned>(i);
            s = 5U;
        }
        else
        {
            while (s % 5U == 0U)
            {
                s /= 5U;
                ++j;
            }
            if (s != 1U)
            {
                return {};
            }
            s = 2U;
        }
        if (j * n >= 128U)
        {
            return {};
        }
        n *= j;
        scale -= static_cast<int>(n);
    }
    boost::int128::uint128_t v {};
    if (n >= 128U || !pow_power(s, n, v))
    {
        return {};
    }
    return {v, scale, true};
}

// x^y where it is exact, for x = sx * 10^ex with sx > 1 and y = m / d in lowest terms with d = 2^a * 5^b > 1:
// x^m = (x^y)^d needs sx = q^d with d dividing ex, and then x^y = q^m * 10^(ex m / d).
constexpr auto pow_exact_root(const boost::int128::uint128_t& sx, const int ex, const boost::int128::uint128_t& sy, const int ey, const bool neg_y) noexcept -> pow_exact_value
{
    const auto k {static_cast<unsigned>(-ey)};
    const auto a {k - (std::min)(k, static_cast<unsigned>(boost::int128::countr_zero(sy)))};
    unsigned fives {};
    auto t {sy};
    while (fives < k && t % 5U == 0U)
    {
        t /= 5U;
        ++fives;
    }
    const unsigned b {k - fives};
    const auto bits {static_cast<unsigned>(128 - countl_zero(sx))};
    if (a >= 7U || b >= 3U)
    {
        return {};
    }
    const unsigned d {(1U << a) * (b == 0U ? 1U : b == 1U ? 5U : 25U)};
    if (d >= bits || ex % static_cast<int>(d) != 0)
    {
        return {};
    }

    // q = floor(sx^(1/d)) by bisection, with q^d < 2^bits
    boost::int128::uint128_t lo {1U};
    boost::int128::uint128_t hi {boost::int128::uint128_t {1U} << ((bits + d - 1U) / d)};
    while (hi - lo > 1U)
    {
        const auto mid {(lo + hi) / 2U};
        boost::int128::uint128_t v {};
        if (!pow_power(mid, d, v) || v > sx)
        {
            hi = mid;
        }
        else
        {
            lo = mid;
        }
    }
    boost::int128::uint128_t v {};
    pow_power(lo, d, v);
    const auto m {sy / (pow10(static_cast<boost::int128::uint128_t>(k)) / d)};
    if (v != sx || m >= 128U)
    {
        return {};
    }
    return pow_exact(lo, ex / static_cast<int>(d), static_cast<unsigned>(m), neg_y);
}

// |log10(x)| * 10^-k in N + 1 words as fx_wide, and the sign of log10(x)
template <int N>
struct pow_log
{
    fx_wide<N> v;
    int k;
    bool neg;
};

// log10(x) for x = sig * 10^e != 1 with p digits. Near 1, w = x - 1 = m * 10^k is exact with m in [0.1, 1), and
// log10(1 + w) = w * ln(1 + w) / w * log10(e) with the ratio from log_ratio keeps all digits. Else |log10(x)| > 0.0016.
template <int N>
constexpr auto pow_log10(const boost::int128::uint128_t& sig, const int e, const int p) noexcept -> pow_log<N>
{
    if (e == -p || e == 1 - p)
    {
        const auto one {pow10(static_cast<boost::int128::uint128_t>(-e))};
        const bool below {sig < one};
        const auto d {below ? one - sig : sig - one};
        if (decimal_tiny(d, e))
        {
            const auto g {log_detail::log_ratio(fx_from_decimal<N>(d, -e), below)};
            const int nd {num_digits(d)};
            return {fx_widen(fx_mul(fx_from_decimal<N>(d, nd), fx_mul(g, fx_top<N>(fx_table::log10_e)))), e + nd, below};
        }
    }
    auto r {log_detail::log_wide_of<N>(sig, p, e, log_detail::log_base::ten)};
    bool neg {};
    r = fx_wide_abs(r, neg);
    return {r, 0, neg};
}

// Whether digits with an error below 8, 10^12 or 10^9 units for three, two or one words can round another way
// than the exact value: the digits below the p + 1 kept ones are near ...0 or ...5, as for an exact value.
template <int words, typename digits_type>
constexpr auto pow_near_boundary(const digits_type& sig, const int p) noexcept -> bool
{
    BOOST_DECIMAL_IF_CONSTEXPR (words == 3)
    {
        return fx_near_boundary(sig);
    }
    constexpr int scale10 {fx_scale10(words)};
    const int len {sig >= pow10(static_cast<digits_type>(scale10)) ? scale10 + 1 : scale10};
    const auto m {pow10(static_cast<digits_type>(len - p - 1)) * 5U};
    const auto r {sig % m};
    const auto w {pow10(static_cast<digits_type>(words == 1 ? 9 : 12))};
    return r < w || r > m - w;
}

enum class pow_kind { digits, overflow, underflow, above_one, below_one };

// The result as digits * 10^exp10 with scale10 + 1 digits, or the kind of a result beyond them
template <int words>
struct pow_parts
{
    fx_digits_t<words> digits;
    int exp10;
    pow_kind kind;
};

// |x|^y for x = sig * 10^e with p digits and y = sy * 10^ey in words: y log10(x) = whole + f gives 10^whole * 10^f,
// or 10^(-whole - 1) * 10^(1 - f) for a negative product. The result overflows or underflows from whole >= bound.
template <int words>
constexpr auto pow_digits(const boost::int128::uint128_t& sig, const int e, const int p, const boost::int128::uint128_t& sy, const int ey,
                          const bool neg_y, const std::uint64_t bound) noexcept -> pow_parts<words>
{
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    // |y log10(x)| = u * |log10(x)| * 10^k > 10^(k - 4) with u in [0.1, 1). For k < 0 it is below 1, where an error of
    // 10^k in fixed point is small enough. Below 2^-(64 words - 20), it cannot change the rounding of 1.
    const auto lg {pow_log10<words>(sig, e, p)};
    const int nds {num_digits(sy)};
    auto t {fx_wide_mul(lg.v, fx_from_decimal<words>(sy, nds))};
    const int k {ey + nds + lg.k};
    const bool neg_t {neg_y != lg.neg};
    if (k > 0 && k < 9)
    {
        std::uint64_t carry {};
        for (int i {}; i <= words; ++i)
        {
            const auto v {static_cast<boost::int128::uint128_t>(t.w[i]) * pow10(static_cast<std::uint64_t>(k)) + carry};
            t.w[i] = v.low;
            carry = v.high;
        }
    }
    else if (k < 0)
    {
        t = fx_wide_mul(t, fx_from_decimal<words>(1U, -k));
    }
    const std::uint64_t whole {(t.w[words] << 2U) | (t.w[words - 1] >> 62U)};
    if (k >= 9 || whole >= bound)
    {
        return {0U, 0, neg_t ? pow_kind::underflow : pow_kind::overflow};
    }
    fx<words> f {};
    bool tiny {whole == 0U && t.w[0] < (UINT64_C(1) << 18U)};
    for (int i {}; i < words; ++i)
    {
        f.w[i] = t.w[i];
        tiny = tiny && (i == 0 || t.w[i] == 0U);
    }
    f.w[words - 1] &= (UINT64_C(1) << 62U) - 1U;
    if (tiny)
    {
        return {0U, 0, neg_t ? pow_kind::below_one : pow_kind::above_one};
    }

    auto exp10 {static_cast<int>(whole)};
    if (neg_t)
    {
        exp10 = -exp10 - 1;
        f = fx_sub(fx_one<words>(), f);
    }
    const auto r10 {exp_detail::exp10_frac(f)};
    return {static_cast<digits_type>(fx_scale(pow10(static_cast<boost::int128::uint128_t>(scale10)) << (r10.a + 1U), fx_shr(r10.m, 1U))),
            exp10 - scale10, pow_kind::digits};
}

// remove_trailing_zeros, which takes six products, with a quick return where the last digit is not zero
template <typename Sig>
constexpr auto pow_trim(const Sig s) noexcept -> remove_trailing_zeros_return<Sig>
{
    if ((s & 1U) != 0U || s % 5U != 0U)
    {
        return {s, 0U};
    }
    return remove_trailing_zeros(s);
}

// As 2^64 = 1 (mod 5), s = high + low (mod 5), without the division of uint128_t
constexpr auto pow_trim(const boost::int128::uint128_t& s) noexcept -> remove_trailing_zeros_return<boost::int128::uint128_t>
{
    if ((s.low & 1U) != 0U || (s.high % 5U + s.low % 5U) % 5U != 0U)
    {
        return {s, 0U};
    }
    return remove_trailing_zeros(s);
}

// x^y for x = sig * 10^e and y = sy * 10^ey in words, rounded once with the sign neg. One word computes it again in
// two near a boundary, and two or three words test for an exact value there.
template <typename T, int words>
constexpr auto pow_words(const boost::int128::uint128_t& sig, const int e, const boost::int128::uint128_t& sy, const int ey,
                         const bool neg_y, const bool neg) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int scale10 {fx_scale10(words)};
    using digits_type = fx_digits_t<words>;

    const auto r {pow_digits<words>(sig, e, p, sy, ey, neg_y, static_cast<std::uint64_t>(std::numeric_limits<T>::max_exponent10 + p + 1))};
    switch (r.kind)
    {
        case pow_kind::overflow:
            return T {1U, 2 * std::numeric_limits<T>::max_exponent10, neg};
        case pow_kind::underflow:
            return T {1U, -2 * std::numeric_limits<T>::max_exponent10, neg};
        case pow_kind::above_one:
            return fx_round<T>(static_cast<digits_type>(pow10(static_cast<boost::int128::uint128_t>(scale10))), -scale10, neg);
        case pow_kind::below_one:
            return exp_detail::below_one<T, digits_type>(scale10, neg);
        default:
            break;
    }
    if (pow_near_boundary<words>(r.digits, p))
    {
        BOOST_DECIMAL_IF_CONSTEXPR (words == 1)
        {
            return pow_words<T, 2>(sig, e, sy, ey, neg_y, neg);
        }
        const auto tx {pow_trim(sig)};
        const auto exact {ey < 0 ? pow_exact_root(tx.trimmed_number, e + static_cast<int>(tx.number_of_removed_zeros), sy, ey, neg_y) : pow_exact_value {}};
        if (exact.ok)
        {
            return T {exact.sig, exact.exp, neg};
        }
    }
    return fx_round<T>(r.digits, r.exp10, neg);
}

// x^y for a finite x != 0 and y = sy * 10^ey != 0 with sy not a multiple of 10 and the sign neg_y, rounded once in
// the current mode
template <typename T>
constexpr auto pow_trimmed(const T x, const boost::int128::uint128_t& sy, const int ey, const bool neg_y) noexcept -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    constexpr int words {fx_words(p)};

    // A negative x has a real power only for an integer y
    const bool neg_x {signbit(x)};
    if (neg_x && ey < 0)
    {
        #ifndef BOOST_DECIMAL_FAST_MATH
        return std::numeric_limits<T>::quiet_NaN();
        #else
        return T {0};
        #endif
    }
    const bool neg {neg_x && ey == 0 && (sy & 1U) != 0U};

    // x = sig * 10^e, a power of 10 for sig = 10^(p - 1) with the exponent e10
    int e {};
    const auto sig_x {frexp10(x, &e)};
    const auto sig {static_cast<boost::int128::uint128_t>(sig_x)};
    const bool power10 {sig == pow10(static_cast<boost::int128::uint128_t>(p - 1))};
    const int e10 {e + p - 1};

    // |x| = 1 gives 1. An integer y = n gives 10^(ex n) for sx = 1, and x * x and 1 / x round once. Else x^n is
    // exact only where it fits in 128 bits, for n < 128.
    if (power10 && e10 == 0)
    {
        return T {1U, 0, neg};
    }
    if (ey >= 0)
    {
        const auto tx {pow_trim(sig_x)};
        const auto sx {static_cast<boost::int128::uint128_t>(tx.trimmed_number)};
        const int ex {e + static_cast<int>(tx.number_of_removed_zeros)};
        const auto n {ey < 6 && sy < 100000U ? sy * pow10(static_cast<boost::int128::uint128_t>(ey)) : boost::int128::uint128_t {UINT64_MAX}};
        if (sx == 1U)
        {
            // Past 10^5, 4 max_exponent10 is far enough out to overflow or underflow
            const int scale {n > 100000U ? (ex < 0 ? -4 : 4) * std::numeric_limits<T>::max_exponent10 : ex * static_cast<int>(n)};
            return T {1U, neg_y ? -scale : scale, neg};
        }
        if (n == 1U)
        {
            return neg_y ? T {1} / x : x;
        }
        if (n == 2U && !neg_y)
        {
            return x * x;
        }
        const auto exact {n < 128U ? pow_exact(sx, ex, static_cast<unsigned>(n), neg_y) : pow_exact_value {}};
        if (exact.ok)
        {
            return T {exact.sig, exact.exp, neg};
        }
    }
    else if (power10 && ey > -39)
    {
        // 10^(ex y) is exact for an integer ex * y
        const auto prod10 {static_cast<boost::int128::uint128_t>(e10 < 0 ? -e10 : e10) * sy};
        const auto den {pow10(static_cast<boost::int128::uint128_t>(-ey))};
        if (prod10 % den == 0U)
        {
            const auto n {prod10 / den};
            const int scale {n > 100000U ? 4 * std::numeric_limits<T>::max_exponent10 : static_cast<int>(n)};
            return T {1, (e10 < 0) != neg_y ? -scale : scale};
        }
    }
    return pow_words<T, words>(sig, e, sy, ey, neg_y, neg);
}

// pow_trimmed for y = sig_y * 10^e_y, with the trailing zeros of sig_y removed in its own type
template <typename T, typename Sig>
constexpr auto pow_fixed(const T x, const Sig sig_y, const int e_y, const bool neg_y) noexcept -> T
{
    const auto ty {pow_trim(sig_y)};
    return pow_trimmed(x, static_cast<boost::int128::uint128_t>(ty.trimmed_number), e_y + static_cast<int>(ty.number_of_removed_zeros), neg_y);
}

} // namespace pow_detail

} // namespace detail

BOOST_DECIMAL_EXPORT template <typename T, typename IntegralType>
constexpr auto pow(const T b, const IntegralType p) noexcept
    BOOST_DECIMAL_REQUIRES_TWO_RETURN(detail::is_decimal_floating_point_v, T, detail::is_integral_v, IntegralType, T)
{
    using local_integral_type = IntegralType;

    constexpr T zero { 0, 0 };
    constexpr T one  { 1, 0 };

    T result { };

    const auto fpc_x = fpclassify(b);

    const auto p_is_odd = (static_cast<local_integral_type>(p & 1) != static_cast<local_integral_type>(0));

    if (p == static_cast<local_integral_type>(0))
    {
        // pow(base, +/-0) returns 1 for any base, even when base is NaN.

        // Excluded from LCOV since it's apparently optimized away or otherwise
        // missing from LCOV. Verified this line is covered in the unit tests.

        result = one; // LCOV_EXCL_LINE
    }
    else if (fpc_x == FP_ZERO)
    {
        // pow(  +0, exp), where exp is a negative odd  integer, returns +infinity.
        // pow(  -0, exp), where exp is a negative odd  integer, returns +infinity.
        // pow(+/-0, exp), where exp is a negative even integer, returns +infinity.

        // pow(  +0, exp), where exp is a positive odd  integer, returns +0.
        // pow(  -0, exp), where exp is a positive odd  integer, returns -0.
        // pow(+/-0, exp), where exp is a positive even integer, returns +0.

        result = ((p < static_cast<local_integral_type>(0)) ? std::numeric_limits<T>::infinity() : ((p_is_odd && signbit(b)) ? -zero : zero));
    }
    #ifndef BOOST_DECIMAL_FAST_MATH
    else if (fpc_x == FP_INFINITE)
    {
        if (signbit(b))
        {
            if (p < static_cast<local_integral_type>(0))
            {
                // pow(-infinity, exp) returns -0 if exp is a negative odd integer.
                // pow(-infinity, exp) returns +0 if exp is a negative even integer.

                result = (p_is_odd ? -zero : zero);
            }
            else
            {
                // pow(-infinity, exp) returns -infinity if exp is a positive odd integer.
                // pow(-infinity, exp) returns +infinity if exp is a positive even integer.

                result = (p_is_odd ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity());
            }
        }
        else
        {
            // pow(+infinity, exp) returns +0 for any negative exp.
            // pow(+infinity, exp) returns +infinity for any positive exp.

            result = ((p < static_cast<local_integral_type>(0)) ? zero : std::numeric_limits<T>::infinity());
        }
    }
    else if (fpc_x == FP_NAN)
    {
        // Excluded from LCOV since it's apparently optimized away or otherwise
        // missing from LCOV. Verified this line is covered in the unit tests.

        result = std::numeric_limits<T>::quiet_NaN(); // LCOV_EXCL_LINE
    }
    #endif
    else
    {
        using local_unsigned_integral_type = std::make_unsigned_t<IntegralType>;

        const auto up {static_cast<local_unsigned_integral_type>(p < static_cast<local_integral_type>(0) ? static_cast<local_unsigned_integral_type>(~p) + 1U : static_cast<local_unsigned_integral_type>(p))};

        using magnitude_type = std::conditional_t<(sizeof(IntegralType) > sizeof(std::uint64_t)), boost::int128::uint128_t, std::uint64_t>;

        result = detail::pow_detail::pow_fixed(b, static_cast<magnitude_type>(up), 0, p < static_cast<local_integral_type>(0));
    }

    return result;
}

BOOST_DECIMAL_EXPORT template <typename T>
constexpr auto pow(const T x, const T a) noexcept
    BOOST_DECIMAL_REQUIRES(detail::is_decimal_floating_point_v, T)
{
    constexpr T zero { 0, 0 };
    constexpr T one  { 1, 0 };

    auto result = zero;

    const auto fpc_a = fpclassify(a);
    const auto fpc_x = fpclassify(x);

    if (fpc_a == FP_ZERO)
    {
        // pow(base, +/-0) returns 1 for any base, even when base is NaN.

        result = one;
    }
    else if ((fpc_x == FP_NORMAL || fpc_x == FP_SUBNORMAL) && (fpc_a == FP_NORMAL || fpc_a == FP_SUBNORMAL))
    {
        int e_a {};
        const auto sig_a {frexp10(a, &e_a)};

        result = detail::pow_detail::pow_fixed(x, sig_a, e_a, signbit(a));
    }
    else if (static_cast<int>(a) == a)
    {
        result = pow(x, static_cast<int>(a));
    }
    else
    {
        if ((fpc_x == FP_NAN) || (fpc_a == FP_NAN))
        {
            // This line is known to be covered by tests.
            result = std::numeric_limits<T>::quiet_NaN(); // LCOV_EXCL_LINE
        }
        else if (fpc_x == FP_ZERO)
        {
            #ifndef BOOST_DECIMAL_FAST_MATH
            if ((fpc_a == FP_NORMAL) || (fpc_a == FP_INFINITE))
            {
                // pow(+/-0, exp), where exp is negative and finite, returns +infinity.
                // pow(+/-0, exp), where exp is positive non-integer, returns +0.

                // pow(+/-0, -infinity) returns +infinity.
                // pow(+/-0, +infinity) returns +0.

                result = (signbit(a) ? std::numeric_limits<T>::infinity() : zero);
            }
            #else
            if (fpc_a == FP_NORMAL)
            {
                result = T{0};
            }
            #endif
        }
        #ifndef BOOST_DECIMAL_FAST_MATH
        else if (fpc_x == FP_INFINITE)
        {
            if ((fpc_a == FP_NORMAL) || (fpc_a == FP_INFINITE))
            {
                // pow(+infinity, exp) returns +0 for any negative exp.
                // pow(-infinity, exp) returns +infinity for any positive exp.

                result = (signbit(a) ? zero : std::numeric_limits<T>::infinity());
            }
        }
        #endif
        else
        {
            #ifndef BOOST_DECIMAL_FAST_MATH
            if (fpc_a == FP_INFINITE)
            {
                result =
                    (
                          (fabs(x) < one) ? (signbit(a) ? std::numeric_limits<T>::infinity() : zero)
                        : (fabs(x) > one) ? (signbit(a) ? zero : std::numeric_limits<T>::infinity())
                        : one
                    );
            }
            #endif
        }
    }

    return result;
}

} // namespace decimal
} // namespace boost

#endif // BOOST_DECIMAL_DETAIL_CMATH_POW_HPP
