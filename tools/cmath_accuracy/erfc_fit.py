#!/usr/bin/env python3
# Copyright 2026 Shen-Ta Hsieh
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Prints the table erfc_d32, erfc_d64 or erfc_d128 of impl/fixed_point_tables.hpp. On each interval [c - h, c + h],
# g(x) = erfc(x) exp(x^2) gets the Chebyshev interpolation of the lowest degree from 4 up whose relative error on 201
# points is below the limit of the type. The table has |q_k| of sum q_k u^k with u = (x - c) / h, rounded to 2^(2 - 64N)
# for N words. Needs mpmath. Usage: erfc_fit.py d32|d64|d128 > table.txt, with the degrees on the standard error.
# With --float for d64 or d128, it prints the bodies of erf_float_mant_for and erf_float_table_for in erf.hpp instead,
# which test/github_issue_erf.cpp makes again from the table in the same steps.

import sys
from mpmath import mp, mpf, erfc, exp, cos, pi, fsum, nint

mp.dps = 150
# The words, the limit of the relative error, and the intervals of each type
TYPES = {
    'd32': (1, mpf('1e-21'), [(1, 2), (2, 4), (4, 8)]),
    'd64': (2, mpf('1e-34'), [(1, 2), (2, 4), (4, 8), (8, 10)]),
    'd128': (3, mpf('1e-50'), [(1, 2), (2, 4), (4, 8), (8, 12), (12, 16)]),
}


def g(x):
    return erfc(x) * exp(x * x)


def monomials(c):
    """sum c_j T_j(u) as the coefficients of u^k."""
    t = [[mpf(1)], [mpf(0), mpf(1)]]
    for j in range(2, len(c)):
        a = [mpf(0)] + [2 * v for v in t[j - 1]]
        b = t[j - 2] + [mpf(0)] * (len(a) - len(t[j - 2]))
        t.append([a[k] - b[k] for k in range(len(a))])
    m = [mpf(0)] * len(c)
    for j in range(len(c)):
        for k, v in enumerate(t[j]):
            m[k] += c[j] * v
    return m


def fit(a, b, n):
    """The interpolation of degree n at the Chebyshev nodes of [a, b]."""
    c0, h = mpf(a + b) / 2, mpf(b - a) / 2
    count = n + 1
    xs = [cos(pi * (k + mpf(1) / 2) / count) for k in range(count)]
    fs = [g(c0 + h * x) for x in xs]
    c = [2 * fsum(fs[k] * cos(pi * j * (k + mpf(1) / 2) / count) for k in range(count)) / count for j in range(count)]
    c[0] /= 2
    return monomials(c)


def error(a, b, m):
    c0, h = mpf(a + b) / 2, mpf(b - a) / 2
    worst = mpf(0)
    for i in range(201):
        u = -1 + mpf(i) / 100
        worst = max(worst, abs(fsum(m[k] * u ** k for k in range(len(m))) / g(c0 + h * u) - 1))
    return worst


# The first coefficient and the count of the intervals [1, 2], [2, 4], [4, 8] and the next one, as in erf_interval_of
FLOAT_INTERVALS = {2: [(0, 31), (31, 36), (67, 40), (107, 27)], 3: [(0, 43), (43, 50), (93, 57), (150, 49)]}


def fx_mul(a, b, n):
    """fx_mul of N words: the products of words i + k >= N - 1, truncated to the scale 2^(64N - 2)."""
    mask = 2 ** 64 - 1
    total = 0
    for i in range(n):
        for k in range(n):
            if i + k >= n - 1:
                total += ((a >> (64 * i)) & mask) * ((b >> (64 * k)) & mask) << (64 * (i + k))
    return (total >> (64 * n - 2)) & (2 ** (64 * n) - 1)


def float_l(c, n):
    """l for c in [2^(64N - 66 - l), 2^(64N - 65 - l)), as in github_issue_erf.cpp."""
    return 64 * n - c.bit_length() - 1


def print_float(words, q):
    """The mantissas m, the shifts, l, and the last and split terms for |u| < (j + 1) / 8 of each interval."""
    mant, shift = [0] * sum(c for _, c in FLOAT_INTERVALS[words]), [0] * sum(c for _, c in FLOAT_INTERVALS[words])
    ls, lasts, splits = [], [], []
    for first, count in FLOAT_INTERVALS[words]:
        nxt = 0
        for k in range(first + count - 1, first - 1, -1):
            own = float_l(q[k], words)
            l = own if k + 1 == first + count or own < nxt else nxt - 1
            mant[k] = ((q[k] << l) >> 64) & (2 ** (64 * (words - 1)) - 1)
            shift[k] = 0 if k + 1 == first + count else nxt - l - 1
            nxt = l
        ls.append(nxt)
        last_row, split_row = [], []
        for j in range(8):
            um = (j + 1) << (59 + 64 * (words - 1))
            pw = [2 ** (64 * words - 2)]
            for _ in range(1, count):
                pw.append(fx_mul(pw[-1], um, words))
            tail, bound, last = 0, 1 << (60 - nxt), first + count - 1
            while last > first + 1:
                tail = (tail + fx_mul(q[last], pw[last - first], words)) % 2 ** (64 * words)
                if bound < tail:
                    break
                last -= 1
            split = first
            while words == 3 and split < last and float_l(fx_mul(q[split], pw[split - first], words), words) < nxt + 65:
                split += 1
            last_row.append(last)
            split_row.append(split)
        lasts.append(last_row)
        splits.append(split_row)
    print('        {')
    if words == 2:
        cells = [f'UINT64_C(0x{m:016X})' for m in mant]
        for i in range(0, len(cells), 3):
            print('            ' + ', '.join(cells[i:i + 3]) + (',' if i + 3 < len(cells) else ''))
    else:
        for i, m in enumerate(mant):
            print(f'            boost::int128::uint128_t {{UINT64_C(0x{m >> 64:016X}), UINT64_C(0x{m & (2 ** 64 - 1):016X})}}'
                  + (',' if i + 1 < len(mant) else ''))
    print('        }')
    print('        {' + ', '.join(str(v) for v in ls) + '},')
    for table in (lasts, splits):
        text = ['{' + ', '.join(str(v) for v in row) + '}' for row in table]
        print('        {' + text[0] + ', ' + text[1] + ',')
        print('         ' + text[2] + ', ' + text[3] + '},')
    print('        {')
    cells = [str(v) for v in shift]
    for i in range(0, len(cells), 36):
        print('            ' + ', '.join(cells[i:i + 36]) + (',' if i + 36 < len(cells) else ''))
    print('        }')


def main():
    flt = sys.argv[2:] == ['--float']
    if len(sys.argv) != (3 if flt else 2) or sys.argv[1] not in TYPES or (flt and sys.argv[1] == 'd32'):
        sys.exit('usage: erfc_fit.py d32|d64|d128 [--float]')
    words, limit, intervals = TYPES[sys.argv[1]]
    table = []
    for a, b in intervals:
        n = 4
        m = fit(a, b, n)
        while error(a, b, m) >= limit:
            n += 1
            m = fit(a, b, n)
        # The library takes q_k as alternate in sign, and each partial sum of Horner's rule for u in [0, 1] as positive
        assert all((m[k] > 0) == (k % 2 == 0) for k in range(len(m))), (a, b)
        for i in range(101):
            s = mpf(0)
            for k in range(len(m) - 1, -1, -1):
                s = abs(m[k]) - mpf(i) / 100 * s
                assert s > 0, (a, b, k, i)
        print(f'[{a}, {b}]: degree {n}', file=sys.stderr)
        table += [abs(v) for v in m]
    if flt:
        print_float(words, [int(nint(v * mpf(2) ** (64 * words - 2))) for v in table])
        return
    print(f'    static constexpr fx<{words}> erfc_{sys.argv[1]}[{len(table)}] =')
    print('    {')
    for i, v in enumerate(table):
        q = int(nint(v * mpf(2) ** (64 * words - 2)))
        cells = ', '.join(f'UINT64_C(0x{(q >> (64 * w)) & (2 ** 64 - 1):016X})' for w in range(words))
        print(f'        {{{{{cells}}}}}' + (',' if i + 1 < len(table) else ''))
    print('    };')


if __name__ == '__main__':
    main()
