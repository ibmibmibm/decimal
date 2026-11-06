// Copyright 2026 Shen-Ta Hsieh
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Cuts each error bound in the headers that bound_mutate.py writes, and shows which inputs then give another result.
// See README.md in this directory.

#include <boost/decimal.hpp>
#include <bound_runtime.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using namespace boost::decimal;

namespace {

const rounding_mode modes[] {rounding_mode::fe_dec_to_nearest, rounding_mode::fe_dec_upward,
                             rounding_mode::fe_dec_downward, rounding_mode::fe_dec_to_nearest_from_zero,
                             rounding_mode::fe_dec_toward_zero};
const unsigned cuts[] {1U, 2U, 4U, 99U};

template <typename U>
auto digits_of(U s) -> std::string
{
    if (s == 0U)
    {
        return "0";
    }
    std::string r;
    for (; s != 0U; s /= 10U)
    {
        r.insert(r.begin(), static_cast<char>('0' + static_cast<int>(s % 10U)));
    }
    return r;
}

template <typename T>
auto show(const T x) -> std::string
{
    if (isnan(x))
    {
        return "nan";
    }
    if (isinf(x))
    {
        return signbit(x) ? "-inf" : "inf";
    }
    int e {};
    const auto s {frexp10(x, &e)};
    return std::string(signbit(x) ? "-" : "") + digits_of(s) + "e" + std::to_string(e);
}

template <typename T>
auto parse(const std::string& s) -> T
{
    if (s == "nan")
    {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (s == "inf" || s == "-inf")
    {
        return s[0] == '-' ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::infinity();
    }
    T x {};
    const bool neg {s[0] == '-'};
    from_chars(s.data() + (neg ? 1 : 0), s.data() + s.size(), x);
    return neg ? -x : x;
}

template <typename T>
struct arguments
{
    T a;
    T b;
    T c;
    unsigned n;
    unsigned m;
};

template <typename T>
using function = T (*)(const arguments<T>&);

// The functions by the names in the records of rounding_grid. A name can have a suffix, such as beta_x.
template <typename T>
auto functions() -> const std::map<std::string, function<T>>&
{
    using A = arguments<T>;
    static const std::map<std::string, function<T>> table {
        {"exp", [](const A& x) -> T { return exp(x.a); }},
        {"expm1", [](const A& x) -> T { return expm1(x.a); }},
        {"exp2", [](const A& x) -> T { return exp2(x.a); }},
        {"log", [](const A& x) -> T { return log(x.a); }},
        {"log2", [](const A& x) -> T { return log2(x.a); }},
        {"log10", [](const A& x) -> T { return log10(x.a); }},
        {"log1p", [](const A& x) -> T { return log1p(x.a); }},
        {"sin", [](const A& x) -> T { return sin(x.a); }},
        {"cos", [](const A& x) -> T { return cos(x.a); }},
        {"tan", [](const A& x) -> T { return tan(x.a); }},
        {"asin", [](const A& x) -> T { return asin(x.a); }},
        {"acos", [](const A& x) -> T { return acos(x.a); }},
        {"atan", [](const A& x) -> T { return atan(x.a); }},
        {"sinh", [](const A& x) -> T { return sinh(x.a); }},
        {"cosh", [](const A& x) -> T { return cosh(x.a); }},
        {"tanh", [](const A& x) -> T { return tanh(x.a); }},
        {"asinh", [](const A& x) -> T { return asinh(x.a); }},
        {"acosh", [](const A& x) -> T { return acosh(x.a); }},
        {"atanh", [](const A& x) -> T { return atanh(x.a); }},
        {"erf", [](const A& x) -> T { return erf(x.a); }},
        {"erfc", [](const A& x) -> T { return erfc(x.a); }},
        {"tgamma", [](const A& x) -> T { return tgamma(x.a); }},
        {"lgamma", [](const A& x) -> T { return lgamma(x.a); }},
        {"zeta", [](const A& x) -> T { return riemann_zeta(x.a); }},
        {"sqrt", [](const A& x) -> T { return sqrt(x.a); }},
        {"cbrt", [](const A& x) -> T { return cbrt(x.a); }},
        {"comp_ellint_1", [](const A& x) -> T { return comp_ellint_1(x.a); }},
        {"comp_ellint_2", [](const A& x) -> T { return comp_ellint_2(x.a); }},
        {"ellint_1", [](const A& x) -> T { return ellint_1(x.a, x.b); }},
        {"ellint_2", [](const A& x) -> T { return ellint_2(x.a, x.b); }},
        {"atan2", [](const A& x) -> T { return atan2(x.a, x.b); }},
        {"pow", [](const A& x) -> T { return pow(x.a, x.b); }},
        {"hypot", [](const A& x) -> T { return hypot(x.a, x.b); }},
        {"beta", [](const A& x) -> T { return beta(x.a, x.b); }},
        {"hermite", [](const A& x) -> T { return hermite(x.n, x.a); }},
        {"legendre", [](const A& x) -> T { return legendre(x.n, x.a); }},
        {"laguerre", [](const A& x) -> T { return laguerre(x.n, x.a); }},
        {"assoc_laguerre", [](const A& x) -> T { return assoc_laguerre(x.n, x.m, x.a); }},
        {"assoc_legendre", [](const A& x) -> T { return assoc_legendre(x.n, x.m, x.a); }},
        {"pow_int", [](const A& x) -> T { return pow(x.a, static_cast<int>(x.n)); }},
        {"hypot3", [](const A& x) -> T { return hypot(x.a, x.b, x.c); }}};
    return table;
}

template <typename T>
auto find_function(const std::string& name) -> function<T>
{
    const auto& table {functions<T>()};
    auto found {table.find(name)};
    if (found == table.end())
    {
        const auto last {name.rfind('_')};
        found = last == std::string::npos ? table.end() : table.find(name.substr(0, last));
    }
    return found == table.end() ? nullptr : found->second;
}

template <typename T>
auto run(const function<T> f, const arguments<T>& x, const int mode) -> T
{
    fesetround(modes[mode]);
    const T v {f(x)};
    fesetround(modes[0]);
    return v;
}

auto split(const std::string& s, const char mark) -> std::vector<std::string>
{
    std::vector<std::string> parts;
    std::size_t start {0};
    for (auto next {s.find(mark)}; next != std::string::npos; next = s.find(mark, start))
    {
        parts.push_back(s.substr(start, next - start));
        start = next + 1;
    }
    parts.push_back(s.substr(start));
    return parts;
}

template <typename T>
struct record
{
    function<T> f;
    std::string name;
    arguments<T> x;
    int mode;
    std::string want;
};

// A record of rounding_grid is name|arguments|mode|result. Polynomials have n and m first, and pow_int has "int n".
template <typename T>
auto parse_record(const std::string& line, record<T>& r) -> bool
{
    const auto parts {split(line, '|')};
    if (parts.size() != 4U)
    {
        return false;
    }
    r.name = parts[0];
    r.f = find_function<T>(r.name);
    if (r.f == nullptr)
    {
        std::printf("unknown function %s\n", r.name.c_str());
        return false;
    }
    auto list {split(parts[1], ',')};
    r.x = arguments<T> {};
    const auto base {r.name.substr(0, r.name.find('_'))};
    if (base == "hermite" || base == "legendre" || base == "laguerre" || base == "assoc")
    {
        r.x.n = static_cast<unsigned>(std::stoi(list[0]));
        r.x.m = list.size() == 3U ? static_cast<unsigned>(std::stoi(list[1])) : 0U;
        r.x.a = parse<T>(list.back());
    }
    else if (r.name.rfind("pow_int", 0) == 0)
    {
        r.x.a = parse<T>(list[0]);
        r.x.n = static_cast<unsigned>(std::stoi(list[1].substr(4)));
    }
    else
    {
        r.x.a = parse<T>(list[0]);
        r.x.b = list.size() > 1U ? parse<T>(list[1]) : T {};
        r.x.c = list.size() > 2U ? parse<T>(list[2]) : T {};
    }
    r.mode = std::stoi(parts[2]);
    r.want = parts[3];
    return true;
}

void print_site(const int s)
{
    std::printf("site %d calls %llu nears %llu", s, boost_mut::counts[s].calls, boost_mut::counts[s].nears);
}

// Runs the records again with each bound cut, only where the site was next to a rounding boundary. No MPFR is needed.
template <typename T>
auto replay(const int argc, char** argv) -> int
{
    using namespace boost_mut;
    std::vector<record<T>> records;
    std::unordered_set<std::string> seen;
    for (int i {3}; i < argc; ++i)
    {
        std::ifstream in(argv[i]);
        std::string line;
        while (std::getline(in, line))
        {
            record<T> r {};
            if (seen.insert(line).second && parse_record(line, r))
            {
                records.push_back(r);
            }
        }
    }
    seen.clear();

    std::map<int, std::vector<std::size_t>> near_records;
    long differs {0};
    target = -1;
    for (std::size_t i {0}; i < records.size(); ++i)
    {
        const auto& r {records[i]};
        hit_count = 0;
        differs += show(run(r.f, r.x, r.mode)) != r.want ? 1 : 0;
        for (const int s : std::set<int>(hits, hits + hit_count))
        {
            near_records[s].push_back(i);
        }
    }
    std::printf("records %zu, differ without a cut %ld\n", records.size(), differs);

    for (int s {1}; s < max_sites; ++s)
    {
        if (counts[s].calls == 0U)
        {
            continue;
        }
        print_site(s);
        for (const unsigned cut : cuts)
        {
            target = s;
            shift = cut;
            long kills {0};
            for (const auto i : near_records[s])
            {
                const auto& r {records[i]};
                kills += show(run(r.f, r.x, r.mode)) != r.want ? 1 : 0;
            }
            target = -1;
            std::printf(" cut%u=%ld", cut, kills);
        }
        std::puts("");
    }
    return differs == 0 ? 0 : 1;
}

// A range of inputs: the decimal exponents of the first digit of a and of b. Without a range for b, b is zero.
struct input_range
{
    const char* name;
    const char* function_name;
    int a_low;
    int a_high;
    bool a_negative;
    int b_low;
    int b_high;
};

constexpr input_range ranges[] {
    {"exp", "exp", -3, 2, true, 0, 0},           {"expm1", "expm1", -60, 2, true, 0, 0},
    {"exp2", "exp2", -3, 2, true, 0, 0},         {"log", "log", -5, 5, false, 0, 0},
    {"log1p", "log1p", -6, 3, false, 0, 0},      {"sin", "sin", -3, 3, true, 0, 0},
    {"cos", "cos", -3, 3, true, 0, 0},           {"tan", "tan", -3, 3, true, 0, 0},
    {"asin", "asin", -3, -1, true, 0, 0},        {"acos", "acos", -3, -1, true, 0, 0},
    {"atan", "atan", -4, 4, true, 0, 0},         {"atanh", "atanh", -3, -1, true, 0, 0},
    {"sinh", "sinh", -3, 2, true, 0, 0},         {"cosh", "cosh", -3, 2, true, 0, 0},
    {"tanh", "tanh", -3, 2, true, 0, 0},         {"asinh", "asinh", -3, 6, true, 0, 0},
    {"acosh", "acosh", 0, 6, false, 0, 0},       {"erf", "erf", -3, 1, true, 0, 0},
    {"erfc", "erfc", -3, 2, true, 0, 0},         {"lgamma", "lgamma", -3, 5, false, 0, 0},
    {"tgamma", "tgamma", -3, 2, false, 0, 0},    {"zeta", "zeta", -2, 2, false, 0, 0},
    {"ellint_1", "ellint_1", -3, -1, true, -2, 0}, {"ellint_2", "ellint_2", -3, -1, true, -2, 0},
    {"beta", "beta", -2, 3, false, -2, 3},       {"beta_large", "beta", -6, 3, false, 15, 60},
    {"pow", "pow", -3, 3, false, -2, 2}};

template <typename T, typename G>
auto random_input(G& g, const int low, const int high, const bool negative) -> T
{
    constexpr int p {std::numeric_limits<T>::digits10};
    // Half of the inputs have few digits, as results next to a boundary often come from such inputs
    const int used {(g() & 1U) != 0U ? p : 1 + static_cast<int>(g() % 8U)};
    boost::int128::uint128_t s {g() % 9U + 1U};
    for (int i {1}; i < p; ++i)
    {
        s = s * 10U + (i < used ? g() % 10U : 0U);
    }
    const int e {low + static_cast<int>(g() % static_cast<unsigned>(high - low + 1))};
    const T x {s, e - p + 1};
    return negative && (g() & 1U) != 0U ? -x : x;
}

// For each input and mode, prints the smallest cut that changes the result: site|cut|range|a|b|mode|new result|result
template <typename T>
auto hunt(const input_range& range, const double seconds, const unsigned seed) -> int
{
    using namespace boost_mut;
    const auto f {find_function<T>(range.function_name)};
    std::mt19937_64 g(seed);
    std::map<std::pair<int, unsigned>, long> kills;
    const auto start {std::chrono::steady_clock::now()};
    unsigned long n {0};
    const bool two {range.b_low != 0 || range.b_high != 0};
    for (;; ++n)
    {
        const std::chrono::duration<double> used {std::chrono::steady_clock::now() - start};
        if ((n & 255U) == 0U && used.count() > seconds)
        {
            break;
        }
        arguments<T> x {};
        x.a = random_input<T>(g, range.a_low, range.a_high, range.a_negative);
        x.b = two ? random_input<T>(g, range.b_low, range.b_high, false) : T {};
        target = -1;
        hit_count = 0;
        (void)run(f, x, 0);
        for (const int s : std::set<int>(hits, hits + hit_count))
        {
            for (int mode {0}; mode < 5; ++mode)
            {
                target = -1;
                const auto base {show(run(f, x, mode))};
                for (const unsigned cut : cuts)
                {
                    target = s;
                    shift = cut;
                    const auto v {show(run(f, x, mode))};
                    target = -1;
                    if (v != base)
                    {
                        // A few inputs for each pair are enough for a table of tests
                        if (kills[{s, cut}]++ < 12)
                        {
                            std::printf("%d|%u|%s|%s|%s|%d|%s|%s\n", s, cut, range.name, show(x.a).c_str(),
                                        show(x.b).c_str(), mode, v.c_str(), base.c_str());
                        }
                        break;
                    }
                }
            }
        }
    }
    std::fprintf(stderr, "%s n=%lu", range.name, n);
    for (int s {1}; s < max_sites; ++s)
    {
        if (counts[s].calls != 0U)
        {
            std::fprintf(stderr, " | %d calls=%llu nears=%llu", s, counts[s].calls, counts[s].nears);
            for (const unsigned cut : cuts)
            {
                const auto found {kills.find({s, cut})};
                if (found != kills.end())
                {
                    std::fprintf(stderr, " cut%u=%ld", cut, found->second);
                }
            }
        }
    }
    std::fputs("\n", stderr);
    return 0;
}

template <typename T>
auto start(const int argc, char** argv) -> int
{
    const std::string part {argv[2]};
    if (part == "replay")
    {
        return replay<T>(argc, argv);
    }
    if (part == "hunt" && argc == 6)
    {
        for (const auto& range : ranges)
        {
            if (range.name == std::string(argv[3]))
            {
                return hunt<T>(range, std::atof(argv[4]), static_cast<unsigned>(std::atoi(argv[5])));
            }
        }
    }
    std::puts("usage: bound_mutants <type> replay <record>...");
    std::puts("       bound_mutants <type> hunt <range> <seconds> <seed>");
    return 2;
}

} // namespace

int main(const int argc, char** argv)
{
    const std::string type {argc > 2 ? argv[1] : ""};
    if (type == "d32")
    {
        return start<decimal32_t>(argc, argv);
    }
    if (type == "d64")
    {
        return start<decimal64_t>(argc, argv);
    }
    if (type == "d128")
    {
        return start<decimal128_t>(argc, argv);
    }
    if (type == "f32")
    {
        return start<decimal_fast32_t>(argc, argv);
    }
    if (type == "f64")
    {
        return start<decimal_fast64_t>(argc, argv);
    }
    if (type == "f128")
    {
        return start<decimal_fast128_t>(argc, argv);
    }
    std::puts("usage: bound_mutants d32|d64|d128|f32|f64|f128 replay|hunt ...");
    return 2;
}
