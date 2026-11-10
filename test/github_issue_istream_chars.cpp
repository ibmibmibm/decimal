// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// operator>> takes only the characters of the number from the stream, so that it works on a stream
// that can not put back what it read before, as a pipe, and keeps the skipws flag

#include <boost/decimal.hpp>
#include <boost/core/lightweight_test.hpp>
#include <algorithm>
#include <cstring>
#include <iterator>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>

using namespace boost::decimal;
using namespace boost::decimal::literals;

// Gives the string four characters at a time, and can not put back a character of an earlier block
class block_buffer : public std::streambuf
{
public:
    explicit block_buffer(const std::string& str) : str_ {str} {}

private:
    auto underflow() -> int_type override
    {
        if (pos_ == str_.size())
        {
            return traits_type::eof();
        }
        const auto n {std::min<std::size_t>(sizeof(block_), str_.size() - pos_)};
        std::memcpy(block_, str_.data() + pos_, n);
        pos_ += n;
        setg(block_, block_, block_ + n);
        return traits_type::to_int_type(block_[0]);
    }

    std::string str_;
    std::size_t pos_ {};
    char block_[4] {};
};

template <typename T>
void check(const std::string& str, const T value, const std::string& rest)
{
    block_buffer buffer {str};
    std::istream is {&buffer};
    T read {};
    is >> read;
    BOOST_TEST_EQ(read, value);
    BOOST_TEST(!is.fail());
    BOOST_TEST(!is.bad());
    is.clear();
    BOOST_TEST_EQ(std::string(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>()), rest);
}

template <typename T>
void test()
{
    check("1.5abcdef", T{1.5_DD}, "abcdef");
    check("-25e-1xyz", T{-2.5_DD}, "xyz");
    check("infx", std::numeric_limits<T>::infinity(), "x");
    check("12345,678", T{12345_DD}, ",678");

    // A list with commas
    block_buffer buffer {"1.25,2.5,-3,4e1"};
    std::istream is {&buffer};
    T value {};
    T sum {};
    char comma {};
    while (is >> value)
    {
        sum += value;
        is >> comma;
    }
    BOOST_TEST_EQ(sum, T{40.75_DD});
    BOOST_TEST(is.eof());
    BOOST_TEST(!is.bad());

    // A space before the number fails with noskipws, as for double
    std::istringstream spaces {" 7"};
    T read {1_DD};
    spaces >> std::noskipws >> read;
    BOOST_TEST(spaces.fail());
    BOOST_TEST_EQ(read, T{1_DD});
    spaces.clear();
    BOOST_TEST_EQ(std::string(std::istreambuf_iterator<char>(spaces), std::istreambuf_iterator<char>()), " 7");

    // The fixed format has no exponent, so the exponent stays in the stream
    block_buffer fixed_buffer {"1e+12345678901234567890"};
    std::istream fixed {&fixed_buffer};
    fixed >> std::fixed >> read;
    BOOST_TEST_EQ(read, T{1_DD});
    BOOST_TEST(fixed.good());
    BOOST_TEST_EQ(std::string(std::istreambuf_iterator<char>(fixed), std::istreambuf_iterator<char>()), "e+12345678901234567890");

    // A wide stream
    std::wistringstream wide {L" 2.5e1x"};
    wide >> read;
    BOOST_TEST_EQ(read, T{25_DD});
    BOOST_TEST_EQ(static_cast<char>(wide.get()), 'x');
}

#ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS

// Gives one character, and then throws
class throwing_buffer : public std::streambuf
{
private:
    auto underflow() -> int_type override
    {
        if (done_)
        {
            throw std::runtime_error("No more characters");
        }
        done_ = true;
        setg(&ch_, &ch_, &ch_ + 1);
        return traits_type::to_int_type(ch_);
    }

    char ch_ {'1'};
    bool done_ {false};
};

// An exception of the stream buffer sets badbit, and leaves operator>> only with badbit in the mask
void test_buffer_exception()
{
    throwing_buffer quiet_buffer;
    std::istream quiet {&quiet_buffer};
    decimal64_t value {7};
    quiet >> value;
    BOOST_TEST(quiet.bad());
    BOOST_TEST_EQ(value, decimal64_t{7});

    throwing_buffer loud_buffer;
    std::istream loud {&loud_buffer};
    loud.exceptions(std::ios_base::badbit);
    BOOST_TEST_THROWS(loud >> value, std::runtime_error);
    BOOST_TEST(loud.bad());
}

// A grouping that throws sets badbit as an exception of the stream buffer. The locale of
// libstdc++ reads the grouping once for its cache, so the facet throws from the second call
class throwing_grouping : public std::numpunct<char>
{
private:
    auto do_grouping() const -> std::string override
    {
        if (calls_++ > 0)
        {
            throw std::runtime_error("No grouping");
        }
        return "\3";
    }

    mutable int calls_ {};
};

void test_facet_exception()
{
    std::istringstream in {"1,234"};
    in.imbue(std::locale(std::locale::classic(), new throwing_grouping));
    decimal64_t value {7};
    in >> value;
    BOOST_TEST_EQ(value, decimal64_t{in.bad() ? 7 : 1234});
}

// The value is in place before an exception for eofbit
void test_eof_exception()
{
    std::istringstream in {"1.5"};
    in.exceptions(std::ios_base::eofbit);
    decimal64_t value {7};
    BOOST_TEST_THROWS(in >> value, std::ios_base::failure);
    BOOST_TEST_EQ(value, (decimal64_t{15, -1}));
}

#endif

// A locale with ',' as the point and '.' between groups of three digits
class point_comma : public std::numpunct<char>
{
private:
    auto do_decimal_point() const -> char override { return ','; }
    auto do_thousands_sep() const -> char override { return '.'; }
    auto do_grouping() const -> std::string override { return "\3"; }
};

// A '.' that is not between two digits is not the point of that locale, and stays in the stream
void test_separator_as_point()
{
    std::istringstream in {"1.e2 1.234,5"};
    in.imbue(std::locale(std::locale::classic(), new point_comma));
    decimal64_t value {7};
    in >> value;
    BOOST_TEST_EQ(value, decimal64_t{1});
    BOOST_TEST(in.good());
    std::string rest;
    in >> rest >> value;
    BOOST_TEST_EQ(rest, ".e2");
    BOOST_TEST_EQ(value, (decimal64_t{12345, -1}));
    BOOST_TEST(in.eof());
    BOOST_TEST(!in.fail());
}

// A wide stream reads the punctuation of its own character type, as operator<< writes it
class wide_point_comma : public std::numpunct<wchar_t>
{
private:
    auto do_decimal_point() const -> wchar_t override { return L','; }
    auto do_thousands_sep() const -> wchar_t override { return L'.'; }
    auto do_grouping() const -> std::string override { return "\3"; }
};

void test_wide_punctuation()
{
    std::wstringstream wide;
    wide.imbue(std::locale(std::locale::classic(), new wide_point_comma));
    wide << decimal64_t{12345, -1};
    BOOST_TEST(wide.str() == L"1.234,5");
    decimal64_t value {7};
    wide >> value;
    BOOST_TEST_EQ(value, (decimal64_t{12345, -1}));
    BOOST_TEST(!wide.fail());
}

int main()
{
    #ifndef BOOST_DECIMAL_DISABLE_EXCEPTIONS
    test_buffer_exception();
    test_facet_exception();
    test_eof_exception();
    #endif
    test_separator_as_point();
    test_wide_punctuation();

    test<decimal32_t>();
    test<decimal64_t>();
    test<decimal128_t>();
    test<decimal_fast32_t>();
    test<decimal_fast64_t>();
    test<decimal_fast128_t>();

    return boost::report_errors();
}
