#!/usr/bin/env python3
# Copyright 2026 Shen-Ta Hsieh
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Prints the tables fit1, fit2_head and fit2_tail of impl/riemann_zeta_impl.hpp and the error of each interval.
# S(x) = zeta(x) - 1 / (x - 1) gets a Chebyshev fit, then each integer evaluation of the library runs at a grid of
# 2001 points against mpmath. Needs mpmath. Usage: zeta_fit.py > tables.txt, with the errors on the standard error.

import sys
from mpmath import mp, mpf, zeta, chebyfit, floor

mp.dps = 90
INTERVALS = [(0, 1), (1, 1), (2, 1), (3, 1), (4, 2), (6, 2), (8, 4), (12, 4), (16, 8), (24, 8), (32, 16), (48, 16)]
# The size of the one-word polynomial of each interval
SIZES_ONE = [12, 12, 12, 12, 14, 14, 16, 14, 18, 16, 18, 16]
MASK = 2 ** 64 - 1


def s_of(x):
    return zeta(x) - 1 / (x - 1) if abs(x - 1) > mpf(10) ** -70 else mp.euler


def rounded(v, bits):
    return int(floor(v * mpf(2) ** bits + mpf(1) / 2))


def one_word_value(c, u):
    """fx_poly_signed: Horner's rule on u / 2^64 with rounded products."""
    acc = c[-1]
    for k in reversed(c[:-1]):
        acc = k + ((acc * u + 2 ** 63) >> 64)
    return acc


def two_word_value(head, tail, h, t, below):
    """zeta_fit in two words: t = |v| times 2^126, with v < 0 in below."""
    u = t >> 62
    acc = tail[-1]
    for c in reversed(tail[:-1]):
        m = (acc * u + 2 ** 63) >> 64
        acc = c - m if below else c + m
    acc <<= h
    for c in reversed(head):
        m = (acc * t) >> 126
        acc = c - m if below else c + m
    return acc


def fit_one(a, w, n):
    p = chebyfit(lambda u: s_of(mpf(a) + w * u), [0, 1], n)[::-1]
    c = [rounded(v, 62) for v in p]
    worst = mpf(0)
    for i in range(2001):
        u = (i * 2 ** 64) // 2000 if i < 2000 else 2 ** 64 - 1
        worst = max(worst, abs(one_word_value(c, u) - s_of(mpf(a) + w * mpf(u) / 2 ** 64) * 2 ** 62))
    return c, worst


def fit_two(a, w):
    """The smallest size with a fit error below 0.3 units of 2^-126, and the smallest head where the tail fits a word."""
    half = mpf(1) / 2
    for n in range(16, 56, 2):
        p, e = chebyfit(lambda v: s_of(mpf(a) + w * (v + half)), [-half, half], n, error=True)
        p = p[::-1]
        if e * 2 ** 126 < mpf(3) / 10:
            break
    else:
        raise RuntimeError('no size for the interval from %d' % a)
    for h in range(4, n + 1):
        limit = mpf(2) ** (h - 63) * mpf(9) / 10
        partial = mpf(0)
        for c in reversed(p[h:]):
            partial = abs(c) + partial / 2
            if partial >= limit:
                break
        else:
            break
    head = [rounded(c, 126) for c in p[:h]]
    tail = [rounded(c, 126 - h) for c in p[h:]] or [0]
    worst = mpf(0)
    for j in range(-1000, 1001):
        t = abs(int(floor(mpf(j) / 2000 * 2 ** 126)))
        v = mpf(-t if j < 0 else t) / 2 ** 126
        worst = max(worst, abs(two_word_value(head, tail, h, t, j < 0) - s_of(mpf(a) + w * (half + v)) * mpf(2) ** 126))
    return head, tail, worst


def word(v):
    return 'UINT64_C(0x%016X)' % (v & MASK)


def main():
    fit1, head, tail = [], [], []
    first_one, first_head, first_tail = [0], [0], [0]
    for (a, w), n in zip(INTERVALS, SIZES_ONE):
        c, e1 = fit_one(a, w, n)
        h, t, e2 = fit_two(a, w)
        fit1 += c
        head += h
        tail += t
        first_one.append(len(fit1))
        first_head.append(len(head))
        first_tail.append(len(tail))
        print('[%d, %d): one word %.2f units of 2^-62, two words %.2f units of 2^-126' % (a, a + w, e1, e2),
              file=sys.stderr, flush=True)
    print('first', first_one)
    print('first_head', first_head)
    print('first_tail', first_tail)
    print('fit1[%d]' % len(fit1))
    print(',\n'.join('        ' + word(v) for v in fit1))
    print('fit2_head[%d]' % len(head))
    print(',\n'.join('        {{%s, %s}}' % (word(v), word(v >> 64)) for v in head))
    print('fit2_tail[%d]' % len(tail))
    print(',\n'.join('        ' + word(v) for v in tail))


if __name__ == '__main__':
    main()
