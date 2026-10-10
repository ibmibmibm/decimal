// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// See: https://github.com/boostorg/decimal/issues/NNNN

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include "quiet_nan.hpp"

using namespace boost::decimal;
using namespace boost::decimal::literals;

// The expected values are erf(x) and erfc(x) rounded once, from mpmath
void test_nearest()
{
    BOOST_TEST_EQ(erfc(0.6054432_DF), 0.3918728_DF);
    BOOST_TEST_EQ(erfc(8_DD), 1.122429717298293E-29_DD);
    BOOST_TEST_EQ(erfc(8_DL), "1.122429717298292707996788844317028E-29"_DL);
    BOOST_TEST_EQ(erfc(7.237657_DF), 1.373516E-24_DF);
    BOOST_TEST_EQ(erf(0.7374917078929143_DD), 0.7030381718352665_DD);
    BOOST_TEST_EQ(erf("0.6909990869580467250525630679297074"_DL), "0.6715398901256361307984119041402326"_DL);
    BOOST_TEST_EQ(erf(0.42_DL), "0.4474676184260253198194105938973533"_DL);
    BOOST_TEST_EQ(erfc(0.3_DL), "0.6713732405408725723610859521332434"_DL);
    BOOST_TEST_EQ(erfc(-0.09_DL), "1.101280593914626883352498163244121"_DL);
    BOOST_TEST_EQ(erfc("1.466714220113811455610092649481663"_DL), "0.03805636406474010858014756910353931"_DL);
    BOOST_TEST_EQ(erfc("5.286644996632435042372313904104612"_DL), "7.636004811898490972481158413652983E-14"_DL);
    BOOST_TEST_EQ(erfc("9.717086581047046900731686691090881"_DL), "5.685445805411021947837466094592123E-43"_DL);
    BOOST_TEST_EQ(erfc("14.44551254136481821536462822834949"_DL), "9.229920158564956817042147395508213E-93"_DL);
    BOOST_TEST_EQ(erfc("18.40693879138063372752601449365155"_DL), "2.188508402769297332913995696907663E-149"_DL);
    BOOST_TEST_EQ(erfc("75.04431114407240626427448017378449"_DL), "1.208210210358260564889915095324682E-2448"_DL);
    BOOST_TEST_EQ(erfc(7.237657_DFF), 1.373516E-24_DFF);
    BOOST_TEST_EQ(erf(0.7374917078929143_DDF), 0.7030381718352665_DDF);
    BOOST_TEST_EQ(erfc("1.466714220113811455610092649481663"_DLF), "0.03805636406474010858014756910353931"_DLF);

    // The first try in fewer words, next to a rounding boundary or not, and where erfc(|x|) cannot change 1 or 2
    BOOST_TEST_EQ(erf(0.8813637618879451_DD), 0.7873958432227577_DD);
    BOOST_TEST_EQ(erf(2.093058595003888_DD), 0.9969239276281441_DD);
    BOOST_TEST_EQ(erfc(-2.533559957247113_DD), 1.999660322696180_DD);
    BOOST_TEST_EQ(erfc(2.5_DD), 4.069520174449589E-4_DD);
    BOOST_TEST_EQ(erfc(1.1_DD), 0.1197949304259183_DD);
    BOOST_TEST_EQ(erfc(15.5_DD), 1.663201640048872E-106_DD);
    BOOST_TEST_EQ(erfc(10.19_DD), 4.423224412798637E-47_DD);
    BOOST_TEST_EQ(erf("1.538478965855649762301814093302143"_DL), "0.9704251294103920334237066174742096"_DL);
    BOOST_TEST_EQ(erfc("-1.901260189359848622879600582536023"_DL), "1.992828804211194806465840121302044"_DL);
    BOOST_TEST_EQ(erf(1.5_DL), "0.9661051464753107270669762616459479"_DL);
    BOOST_TEST_EQ(erfc(-2.5_DL), "1.999593047982555041060435784260025"_DL);
    BOOST_TEST_EQ(erf(1.3106_DL), "0.9361847335846901415609512372770475"_DL);
    BOOST_TEST_EQ(erfc(1.1947_DL), "0.09111198338967185460439260870308278"_DL);
    BOOST_TEST_EQ(erfc(12.1154_DL), "8.310159528653606278842591913656304E-66"_DL);
    BOOST_TEST_EQ(erfc(119_DL), 4.28257066177265941055824E-6153_DL);
    BOOST_TEST_EQ(erf(8_DD), 1_DD);
    BOOST_TEST_EQ(erfc(-8_DD), 2_DD);
    BOOST_TEST_EQ(erf(9.5_DL), 1_DL);
    BOOST_TEST_EQ(erfc(-9.5_DL), 2_DL);

    // Small and large arguments
    BOOST_TEST_EQ(erf(1E-300_DD), 1.128379167095513E-300_DD);
    BOOST_TEST_EQ(erfc(1E-300_DD), 1_DD);
    BOOST_TEST_EQ(erf(10_DD), 1_DD);
    BOOST_TEST_EQ(erfc(30_DD), 2.564656203756112E-393_DD);
    BOOST_TEST_EQ(erfc(-30_DD), 2_DD);
    BOOST_TEST_EQ(erfc(0_DD), 1_DD);
}

void test_directed()
{
    #ifndef BOOST_DECIMAL_NO_CONSTEVAL_DETECTION
    // fesetround has an effect only with the detection of constant evaluation
    fesetround(rounding_mode::fe_dec_downward);
    BOOST_TEST_EQ(erfc(-1.5_DD), 1.966105146475310_DD);
    BOOST_TEST_EQ(erf(-0.8463104_DF), -0.7686403_DF);
    BOOST_TEST_EQ(erf(43.36809_DF), 0.9999999_DF);
    BOOST_TEST_EQ(erfc("0.7628776023409478494943231288202300"_DL), "0.2806448585899707069063940456989828"_DL);
    BOOST_TEST_EQ(erf(8_DD), 0.9999999999999999_DD);
    BOOST_TEST_EQ(erfc(-9.5_DL), "1.999999999999999999999999999999999"_DL);

    fesetround(rounding_mode::fe_dec_upward);
    BOOST_TEST_EQ(erfc(-1.5_DD), 1.966105146475311_DD);
    BOOST_TEST_EQ(erfc(0.5183075_DF), 0.4635600_DF);
    BOOST_TEST_EQ(erf(-4.253308970344476_DD), -0.9999999982021014_DD);
    BOOST_TEST_EQ(erfc(3.791849285048550_DD), 8.208492978872465E-8_DD);
    BOOST_TEST_EQ(erf("-63.26060993105608666610375734038622"_DL), "-0.9999999999999999999999999999999999"_DL);
    BOOST_TEST_EQ(erfc("-0.6661145562372408908860666158430723"_DL), "1.653821817916521188864576677893881"_DL);
    BOOST_TEST_EQ(erfc("71.24424751345447299539316736125789"_DL), "3.400488078544209326870696980098378E-2207"_DL);

    fesetround(rounding_mode::fe_dec_toward_zero);
    BOOST_TEST_EQ(erf(-52.85656_DF), -0.9999999_DF);
    BOOST_TEST_EQ(erfc(0.9457868535668712_DD), 0.1810449356798116_DD);
    BOOST_TEST_EQ(erfc("0.7509075114865893430360491026231036"_DL), "0.2882612961234590548286220881929440"_DL);

    fesetround(rounding_mode::fe_dec_to_nearest);
    #endif
}

// A signaling NaN operand gives a quiet NaN (IEEE 754 6.2)
template <typename T>
void test_snan()
{
    #ifndef BOOST_DECIMAL_FAST_MATH
    for (const auto s : {snan<T>("7"), -snan<T>("7")})
    {
        BOOST_TEST(is_quiet_of(erf(s), s));
        BOOST_TEST(is_quiet_of(erfc(s), s));
    }
    #endif
}

namespace erf_table {

using namespace boost::decimal::detail;
using namespace boost::decimal::detail::erf_detail;

// l for c in [2^(64N - 66 - l), 2^(64N - 65 - l))
template <int N>
auto float_l(const fx<N>& c) -> int
{
    int z {};
    for (int i {N - 1}; i >= 0 && c.w[i] == 0U; --i)
    {
        z += 64;
    }
    return z + (z < 64 * N ? countl_zero(c.w[N - 1 - z / 64]) : 0) - 1;
}

// The top N - 1 words of c * 2^l
template <int N>
auto float_m(const fx<N>& c, const int l) -> erf_mant_t<N>
{
    std::uint64_t w[2] {};
    for (int i {N - 2}; i < N; ++i)
    {
        const int j {i - l / 64};
        const auto r {static_cast<unsigned>(l % 64)};
        const auto hi {j >= 0 ? c.w[j] : UINT64_C(0)};
        const auto lo {j >= 1 ? c.w[j - 1] : UINT64_C(0)};
        w[i - (N - 2)] = r == 0U ? hi : (hi << r) | (lo >> (64U - r));
    }
    return static_cast<erf_mant_t<N>>(boost::int128::uint128_t {w[1], w[0]} >> (N == 2 ? 64 : 0));
}

// The table of erf_float_imp<N> made again from the coefficients c
template <int N, std::size_t S>
void test_float_table(const fx<N> (&c)[S])
{
    const auto& t {erf_float_imp<N>::t};
    for (unsigned i {}; i <= 3U; ++i)
    {
        const auto iv {erf_interval_of(N, i)};
        int next {};
        for (int k {iv.first + iv.count - 1}; k >= iv.first; --k)
        {
            const int own {float_l(c[k])};
            const int l {k + 1 == iv.first + iv.count || own < next ? own : next - 1};
            BOOST_TEST(erf_float_imp<N>::mant.m[k] == float_m(c[k], l));
            BOOST_TEST_EQ(t.shift[k], k + 1 == iv.first + iv.count ? 0 : next - l - 1);
            next = l;
        }
        BOOST_TEST_EQ(t.l[i], next);
        for (unsigned j {}; j < 8U; ++j)
        {
            // The terms after the last sum to below 2^-130 of the scale 2^-l of the first
            fx<N> um {};
            um.w[N - 1] = static_cast<std::uint64_t>(j + 1U) << 59U;
            fx<N> pw[64] {fx_one<N>()};
            for (int k {1}; k < iv.count; ++k)
            {
                pw[k] = fx_mul(pw[k - 1], um);
            }
            fx<N> tail {};
            fx<N> bound {};
            bound.w[0] = UINT64_C(1) << (60 - next);
            int last {iv.first + iv.count - 1};
            while (last > iv.first + 1 && !fx_less(bound, tail = fx_add(tail, fx_mul(c[last], pw[last - iv.first]))))
            {
                --last;
            }
            int split {iv.first};
            while (N == 3 && split < last && float_l(fx_mul(c[split], pw[split - iv.first])) < next + 65)
            {
                ++split;
            }
            BOOST_TEST_EQ(t.last[i][j], last);
            BOOST_TEST_EQ(t.split[i][j], split);
        }
    }
}

} // namespace erf_table

int main()
{
    erf_table::test_float_table(boost::decimal::detail::fx_table::erfc_d64);
    erf_table::test_float_table(boost::decimal::detail::fx_table::erfc_d128);
    test_nearest();
    test_directed();
    test_snan<decimal32_t>();
    test_snan<decimal64_t>();
    test_snan<decimal128_t>();
    test_snan<decimal_fast32_t>();
    test_snan<decimal_fast64_t>();
    test_snan<decimal_fast128_t>();

    return boost::report_errors();
}
