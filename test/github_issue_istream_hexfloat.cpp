// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// operator>> reads the format of the float field, so that it reads what operator<< writes

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

using namespace boost::decimal;
using namespace boost::decimal::literals;

template <typename T>
auto read(const std::string& text, std::ios_base& (*manip)(std::ios_base&), bool& failed) -> T
{
    std::istringstream is {text};
    T value {-1};
    is >> manip >> value;
    failed = is.fail();
    return value;
}

template <typename T>
void test_round_trip(std::ios_base& (*manip)(std::ios_base&))
{
    for (const auto value : {T {1.5_DD}, T {255_DD}, T {-0.001_DD}, T {1e20_DD}, T {0_DD}})
    {
        std::ostringstream os;
        os << manip << value;

        bool failed {};
        BOOST_TEST_EQ(read<T>(os.str(), manip, failed), value);
        BOOST_TEST(!failed);

        std::ostringstream upper_os;
        upper_os << std::uppercase << manip << value;
        BOOST_TEST_EQ(read<T>(upper_os.str(), manip, failed), value);
        BOOST_TEST(!failed);
    }
}

template <typename T>
void test()
{
    test_round_trip<T>(std::defaultfloat);
    test_round_trip<T>(std::fixed);
    test_round_trip<T>(std::scientific);
    test_round_trip<T>(std::hexfloat);

    bool failed {};
    BOOST_TEST_EQ(read<T>("f.fp+01", std::hexfloat, failed), T {255});
    BOOST_TEST(!failed);

    // The hex flag is the base of the integers, and operator<< does not use it too
    BOOST_TEST_EQ(read<T>("10", std::hex, failed), T {10});
    BOOST_TEST(!failed);
    BOOST_TEST_EQ(read<T>("1e1", std::hex, failed), T {10});
    BOOST_TEST(!failed);

    std::istringstream is {"ff 2.5"};
    int n {};
    T value {};
    is >> std::hex >> n >> value;
    BOOST_TEST_EQ(n, 255);
    BOOST_TEST_EQ(value, T {2.5_DD});

    // scientific needs an exponent, and it is not hexfloat
    BOOST_TEST_EQ(read<T>("2.5e+00", std::scientific, failed), T {2.5_DD});
    BOOST_TEST(!failed);
    static_cast<void>(read<T>("1p+00", std::scientific, failed));
    BOOST_TEST(failed);
}

int main()
{
    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    return boost::report_errors();
}
