# Copyright 2026 Shen-Ta Hsieh
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

# Checks each row of the power tables of the conversions with exact integers.
# Run it from the root of the repository: python3 tools/check_power_tables.py

import re
import sys

FAST_FLOAT = 'include/boost/decimal/detail/fast_float/'


def read(path):
    with open(path) as f:
        return f.read()


def array(text, name):
    body = re.search(name + r'\[\] = \{(.*?)\};', text, re.S).group(1)
    return [int(v, 0) for v in re.findall(r'(-?(?:0x[0-9a-fA-F]+|\d+))\)?,?', re.sub(r'UINT64_C\(', '', body))]


def constant(text, name):
    return int(re.search(r'int ' + name + r' \{(-?\d+)\};', text).group(1))


def top_bits(numerator, denominator, up):
    # The top 128 bits of numerator / denominator, and the power of two that scales them back
    shift = 127 - (numerator.bit_length() - denominator.bit_length())
    while True:
        if shift >= 0:
            value, rest = divmod(numerator << shift, denominator)
        else:
            value, rest = divmod(numerator, denominator << -shift)
        if value >> 128:
            shift -= 1
        elif not value >> 127:
            shift += 1
        else:
            return value + (1 if up and rest else 0), rest == 0, -shift


def main():
    failures = 0

    tables = read(FAST_FLOAT + 'significand_tables.hpp')
    high = array(tables, 'significand_64')
    low = array(tables, 'significand_128')
    min_power = constant(tables, 'significand_table_min_power')
    max_power = constant(tables, 'significand_table_max_power')
    min_up_power = constant(tables, 'significand_table_min_up_power')
    max_exact_power = constant(tables, 'significand_table_max_exact_power')
    if len(high) != max_power - min_power + 1 or len(low) != len(high):
        print('significand tables: wrong number of rows')
        failures += 1
    for i, (h, l) in enumerate(zip(high, low)):
        k = min_power + i
        up = min_up_power <= k < 0
        value, exact, _ = top_bits(5 ** k, 1, up) if k >= 0 else top_bits(1, 5 ** -k, up)
        if (h << 64 | l) != value or exact != (0 <= k <= max_exact_power):
            print('significand tables: wrong row for 10^%d' % k)
            failures += 1

    coarse = read(FAST_FLOAT + 'compute_float80_128.hpp')
    coarse_high = array(coarse, 'high')
    coarse_low = array(coarse, 'low')
    coarse_exponent = array(coarse, 'exponent')
    if len(coarse_high) != 34 or len(coarse_low) != 34 or len(coarse_exponent) != 34:
        print('coarse table: wrong number of rows')
        failures += 1
    for i, (h, l, e) in enumerate(zip(coarse_high, coarse_low, coarse_exponent)):
        k = 300 * (i - 17)
        value, _, power = top_bits(5 ** k, 1, False) if k >= 0 else top_bits(1, 5 ** -k, False)
        if (h << 64 | l) != value or e != power + k:
            print('coarse table: wrong row for 10^%d' % k)
            failures += 1

    print('%d rows, %d failures' % (len(high) + len(coarse_high), failures))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
