#!/usr/bin/env python3
# Copyright 2026 Shen-Ta Hsieh
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Writes lists for rounding_grid with arguments whose results are next to a rounding boundary, for the functions that
# are a constant plus a small term: there the arguments next to each other give results that are very close.
# Each case takes a random argument, solves for the argument of the nearest boundary and rounds that to the digits.
# Needs mpmath. Usage: near_ties.py <digits> <count> <seed> [case ...], which writes near_<case>.txt for each case.

import random
import sys
from mpmath import (mp, mpf, floor, log10, findroot, power, erf, erfc, acos, atan, atan2, tanh, expm1, zeta, cos, cosh,
                    exp, ellipf, ellipe)

mp.dps = 150


def decimal_text(x, digits):
    if x == 0:
        return '0e0'
    sign = '-' if x < 0 else ''
    x = abs(x)
    e = int(floor(log10(x))) - (digits - 1)
    m = int(floor(x / mpf(10) ** e + mpf(1) / 2))
    if m >= 10 ** digits:
        m //= 10
        e += 1
    return f'{sign}{m}e{e}'


# The ranges below are those of 34 digits, where the small term is between a unit of the result and 10^-10. Fewer
# digits take a part of each range in proportion, or else all arguments of a case solve for the same boundary.
SCALE = 1.0


def log_range(low, high):
    return lambda: mpf(10) ** (SCALE * random.uniform(low, high))


def linear_range(low, high, power=1.0):
    return lambda: mpf(SCALE ** power * random.uniform(low, high))


def first_kind(phi):
    return lambda k: ellipf(phi, k * k)


def second_kind(phi):
    return lambda k: ellipe(phi, k * k)


# The name, the function of the argument that changes, the line of the list, the range, and whether both signs occur
CASES = [
    ('exp', exp, 'exp|{}|0e0', log_range(-36, -10), True),
    ('exp2', lambda x: power(2, x), 'exp2|{}|0e0', log_range(-36, -10), True),
    ('expm1', expm1, 'expm1|{}|0e0', linear_range(-80, -22), False),
    ('cos', cos, 'cos|{}|0e0', log_range(-18, -5), True),
    ('cosh', cosh, 'cosh|{}|0e0', log_range(-18, -5), True),
    ('tanh', tanh, 'tanh|{}|0e0', linear_range(11, 40), True),
    ('pow_2', lambda y: power(2, y), 'pow|2e0|{}', log_range(-36, -10), True),
    ('pow_07', lambda y: power(mpf(7) / 10, y), 'pow|7e-1|{}', log_range(-36, -10), True),
    ('pow_1e50', lambda y: power(mpf(10) ** 50, y), 'pow|1e50|{}', log_range(-38, -12), True),
    ('pow_1001', lambda y: power(mpf(1001) / 1000, y), 'pow|1001e-3|{}', log_range(-33, -8), True),
    ('erfc', erfc, 'erfc|{}|0e0', log_range(-36, -10), True),
    ('erf', erf, 'erf|{}|0e0', linear_range(5, 8.9, 0.5), True),
    ('erfc_negative', erfc, 'erfc|{}|0e0', linear_range(-8.9, -5, 0.5), False),
    ('acos', acos, 'acos|{}|0e0', log_range(-36, -10), True),
    ('atan', atan, 'atan|{}|0e0', log_range(10, 36), True),
    ('atan2_steep', lambda x: atan2(1, x), 'atan2|1e0|{}', log_range(-36, -10), True),
    ('atan2_flat', lambda y: atan2(y, -1), 'atan2|{}|-1e0', log_range(-36, -10), True),
    ('zeta', zeta, 'zeta|{}|0e0', linear_range(34, 116), False),
    ('zeta_0', zeta, 'zeta|{}|0e0', log_range(-36, -9), True),
    ('ellint_1', first_kind(1), 'ellint_1|{}|1e0', log_range(-18, -5), True),
    ('ellint_2', second_kind(1), 'ellint_2|{}|1e0', log_range(-18, -5), True),
    ('ellint_1_small', first_kind(mpf(3) / 1000), 'ellint_1|{}|3e-3', log_range(-15, -5), True),
    ('ellint_2_small', second_kind(mpf(3) / 1000), 'ellint_2|{}|3e-3', log_range(-15, -5), True),
    ('ellint_1_large', first_kind(mpf(12345)), 'ellint_1|{}|12345e0', log_range(-18, -5), True),
    ('ellint_2_large', second_kind(mpf(12345)), 'ellint_2|{}|12345e0', log_range(-18, -5), True),
]


def main():
    global SCALE
    if len(sys.argv) < 4:
        sys.exit('usage: near_ties.py <digits> <count> <seed> [case ...]')
    digits, count, seed = (int(a) for a in sys.argv[1:4])
    random.seed(seed)
    SCALE = digits / 34
    for name, function, line, argument, both in CASES:
        if sys.argv[4:] and name not in sys.argv[4:]:
            continue
        with open(f'near_{name}.txt', 'w') as out:
            for i in range(count):
                x0 = argument()
                if both and random.random() < 0.5:
                    x0 = -x0
                v = function(x0)
                unit = mpf(10) ** (int(floor(log10(abs(v)))) - (digits - 1))
                # A tie of the nearest modes, then a boundary of the directed modes
                half = i % 2 == 0
                boundary = (floor(v / unit) + (mpf(1) / 2 if half else 0)) * unit
                try:
                    x = findroot(lambda z: function(z) - boundary, x0, tol=mpf(10) ** -130, maxsteps=60)
                except (ValueError, ZeroDivisionError):
                    continue
                # Where the boundary is the constant itself, the solver leaves a zero as noise far below the type
                if x != 0 and abs(x) < mpf(10) ** (-12 * digits):
                    continue
                for mode in ((0, 3) if half else (1, 2, 4)):
                    out.write(f'{line.format(decimal_text(x, digits))}|{mode}\n')


if __name__ == '__main__':
    main()
