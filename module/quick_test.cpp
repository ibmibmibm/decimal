// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A module does not export the macro stdout
#include <cstdio>

// Import the STL if the build supplies the std module
#ifdef BOOST_DECIMAL_USE_STD_MODULE
import std;
#else
#include <iostream>
#include <string_view>
#endif

import boost.decimal;

int main()
{
    using namespace boost::decimal;

    decimal32_t a {2, 0};
    a += 2;
    std::cout << a << std::endl;

    // Every use below needs a name that the module exports
    constexpr decimal32_t negative {2U, 0, construction_sign::negative};

    int quotient {};
    const auto remainder_value {remquo(a, negative, &quotient)};

    int failures {};
    const auto check = [&failures](const bool ok, const char* what)
    {
        if (!ok)
        {
            std::cout << "Failed: " << what << std::endl;
            ++failures;
        }
    };

    char buffer[64] {};
    check(boost::decimal::snprintf(buffer, sizeof(buffer), "%H", remainder_value) == 1, "snprintf");
    check(std::string_view(buffer) == "0" && quotient == -2, "remquo");
    check(boost::decimal::printf("%H\n", normalize(a)) == 2, "printf");
    check(boost::decimal::fprintf(stdout, "%H\n", remainder_value) == 2, "fprintf");

    // The fixed format of 10^300 is longer than the local buffer of snprintf
    char wide_buffer[512] {};
    const auto wide_length {boost::decimal::snprintf(wide_buffer, sizeof(wide_buffer), "%.0Df", decimal64_t {1, 300})};
    check(wide_length == 301 && wide_buffer[0] == '1' && wide_buffer[300] == '0', "snprintf of a long fixed value");

    std::cout << buffer << ' ' << decompose(a).full_significand() << ' '
              << formatting_limits<decimal32_t>::scientific_format_max_chars << std::endl;

    // A qualified name does not find a function that only has a friend declaration
    const boost::decimal::uint128_t wide {1U};
    std::cout << boost::decimal::isnan(a) << boost::decimal::signbit(negative) << boost::decimal::quantexpd32(a)
              << boost::decimal::samequantumd32(a, negative) << wide.low << std::endl;
    check(a == 4 && !boost::decimal::isnan(a) && boost::decimal::signbit(negative), "values");

    // The importer instantiates the math templates, and sin reduces a large argument
    constexpr decimal64_t tolerance {1, -15};
    check(fabs(exp(decimal64_t {1}) - decimal64_t {2718281828459045, -15}) < tolerance, "exp");
    check(fabs(sin(decimal64_t {1, 22}) - decimal64_t {-8522008497671888, -16}) < tolerance, "sin");

    return failures;
}
