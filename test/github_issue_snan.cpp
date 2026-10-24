// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <limits>

using namespace boost::decimal;

// x is the quiet form of the signaling NaN s, with the same sign and payload
template <typename T>
bool is_quiet_of(const T x, const T s)
{
    return isnan(x) && !issignaling(x) && signbit(x) == signbit(s) && read_payload(x) == read_payload(s);
}

// A signaling NaN operand gives a quiet NaN (IEEE 754 6.2)
template <typename T>
void test_quiet_result()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const T one {1};

    for (const auto s : {snan<T>("7"), -snan<T>("7")})
    {
        BOOST_TEST(is_quiet_of(sqrt(s), s));
        BOOST_TEST(is_quiet_of(cbrt(s), s));
        BOOST_TEST(is_quiet_of(logb(s), s));
        BOOST_TEST(is_quiet_of(cosh(s), s));
        BOOST_TEST(is_quiet_of(ceil(s), s));
        BOOST_TEST(is_quiet_of(floor(s), s));
        BOOST_TEST(is_quiet_of(trunc(s), s));
        BOOST_TEST(is_quiet_of(round(s), s));
        BOOST_TEST(is_quiet_of(rint(s), s));
        BOOST_TEST(is_quiet_of(nearbyint(s), s));
        BOOST_TEST(is_quiet_of(roundeven(s), s));
        BOOST_TEST(is_quiet_of(ldexp(s, 1), s));
        BOOST_TEST(is_quiet_of(rescale(s, 2), s));
        BOOST_TEST(is_quiet_of(remainder(s, one), s));
        BOOST_TEST(is_quiet_of(remainder(one, s), s));
        BOOST_TEST(is_quiet_of(nextafter(s, one), s));
        BOOST_TEST(is_quiet_of(nextafter(one, s), s));

        int quo {};
        BOOST_TEST(is_quiet_of(remquo(s, one, &quo), s));
        BOOST_TEST(is_quiet_of(remquo(one, s, &quo), s));

        // The NaN goes before the infinite dividend
        const T inf {std::numeric_limits<T>::infinity()};
        BOOST_TEST(is_quiet_of(remainder(inf, s), s));
        BOOST_TEST(is_quiet_of(remquo(inf, s, &quo), s));
        BOOST_TEST(is_quiet_of(remainder(s, inf), s));

        T ip {};
        BOOST_TEST(is_quiet_of(modf(s, &ip), s));
        BOOST_TEST(is_quiet_of(ip, s));
    }
    #endif
}

// abs only clears the sign bit, also of a NaN, and keeps a signaling NaN (IEEE 754 5.5.1)
template <typename T>
void test_abs()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    const T q {std::numeric_limits<T>::quiet_NaN()};
    const T s {snan<T>("7")};

    BOOST_TEST(isnan(abs(-q)) && !signbit(abs(-q)));
    BOOST_TEST(isnan(fabs(-q)) && !signbit(fabs(-q)));
    BOOST_TEST(issignaling(abs(-s)) && !signbit(abs(-s)) && read_payload(abs(-s)) == 7U);
    BOOST_TEST(issignaling(fabs(-s)) && !signbit(fabs(-s)));
    #endif
}

// The result has no reserved bits (IEEE 754 6.2.3). 0x7E10007B is a signaling NaN with the payload 123 and one bit
// set between the signaling bit and the payload field. read_payload ignores that bit, so the test reads the encoding
void test_canonical_result()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto bits : {UINT32_C(0x7E10007B), UINT32_C(0x7C10007B), UINT32_C(0xFE10007B)})
    {
        const auto x {from_bits(bits)};
        const auto expected {(bits & UINT32_C(0x80000000)) | UINT32_C(0x7C00007B)};
        BOOST_TEST_EQ(to_bits(sqrt(x)), expected);
        BOOST_TEST_EQ(to_bits(ceil(x)), expected);

        decimal32_t ip {};
        BOOST_TEST_EQ(to_bits(modf(x, &ip)), expected);
        BOOST_TEST_EQ(to_bits(ip), expected);
    }
    #endif
}

template <typename T>
void test()
{
    test_quiet_result<T>();
    test_abs<T>();
}

int main()
{
    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();
    test_canonical_result();

    return boost::report_errors();
}
