#!/usr/bin/env python3
# Copyright 2026 Shen-Ta Hsieh
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#
# Copies the headers and marks each check of an error bound, so that bound_mutants can cut each bound.
# See README.md in this directory.

import os
import re
import shutil
import sys

CHECKS = ('fx_near_boundary fx_try_round fx_try_round_word exp_try_round exp_try_scale fx_wide_try_round fx_try_scale2 '
          'trig_try_round trig_near_boundary trig_near_boundary_word pow_near_boundary beta_near_boundary gamma_near_boundary '
          'gamma_try_round ellint_near_boundary atan2_near_boundary fx_unsigned_try_round exp_near_boundary').split()
# These functions only pass a bound to another check, so the callers of them are the sites
WRAPPERS = re.compile(r'constexpr auto (fx_try_round_word|exp_try_round|exp_try_scale|fx_wide_try_round|fx_try_scale2|'
                      r'trig_try_round|trig_near_boundary|pow_near_boundary|gamma_try_round|fx_unsigned_try_round|exp_near_boundary)\b')
CALL = re.compile(r'(?<![\w:])((?:\w+::)*(?:' + '|'.join(CHECKS) + r'))(<[^()]*?>)?\(')
FUNCTION = re.compile(r'\n(?:\w[^\n;]*?)?constexpr auto (\w+)|\n\s*auto (\w+)\(')


def closing(text, start, opening, closer):
    """The index of the bracket that closes the bracket at start."""
    depth = 0
    i = start
    while True:
        if text[i] == opening:
            depth += 1
        elif text[i] == closer:
            depth -= 1
            if depth == 0:
                return i
        i += 1


def wrapper_bodies(text):
    bodies = []
    for match in WRAPPERS.finditer(text):
        brace = text.find('{', match.end())
        semicolon = text.find(';', match.end())
        if semicolon == -1 or brace < semicolon:
            bodies.append((brace, closing(text, brace, '{', '}')))
    return bodies


def mark_file(path, name, sites):
    text = open(path).read()
    if not any(check in text for check in CHECKS):
        return
    bodies = wrapper_bodies(text)
    out = []
    pos = 0
    for match in CALL.finditer(text):
        before = text[text.rfind('\n', 0, match.start()) + 1:match.start()]
        if before.rstrip().endswith('auto') or before.lstrip().startswith('//'):
            continue
        if match.start() < pos or any(a <= match.start() <= b for a, b in bodies):
            continue
        end = closing(text, match.end() - 1, '(', ')') + 1
        site = len(sites) + 1
        out.append(text[pos:match.start()])
        out.append(f'(::boost_mut::set({site}), {text[match.start():end]})')
        functions = list(FUNCTION.finditer(text[:match.start()]))
        function = (functions[-1].group(1) or functions[-1].group(2)) if functions else '?'
        line = text.count('\n', 0, match.start()) + 1
        sites.append((site, f'{name}:{line}', function, match.group(1).split('::')[-1]))
        pos = end
    out.append(text[pos:])
    open(path, 'w').write(''.join(out))


def replace(root, name, old, new):
    path = os.path.join(root, name)
    text = open(path).read()
    if old not in text:
        sys.exit(f'{name} has no text "{old}"')
    open(path, 'w').write(text.replace(old, new, 1))


def main():
    if len(sys.argv) != 3:
        sys.exit('usage: bound_mutate.py <include directory> <new directory>')
    include, output = sys.argv[1:]
    shutil.copytree(include, output)
    root = os.path.join(output, 'boost/decimal/detail/cmath')
    sites = []
    for folder, _, files in sorted(os.walk(root)):
        for file in sorted(files):
            path = os.path.join(folder, file)
            mark_file(path, os.path.relpath(path, root), sites)

    for name, digits in (('atan2.hpp', 'digits'), ('pow.hpp', 'sig'), ('beta.hpp', 'digits'), ('impl/gamma_impl.hpp', 'digits'),
                         ('impl/fixed_point.hpp', 'sig')):
        replace(root, name, 'return r < w || r > m - w;', f'return ::boost_mut::near({digits}, r, w, m);')
    replace(root, 'impl/fixed_point.hpp', 'if (r < w || r > 5U * unit - w || e + k + 1 < etiny_v<T>)',
            'if (::boost_mut::near(sig, r, w, 5U * unit) || e + k + 1 < etiny_v<T>)')
    replace(root, 'impl/trig_reduce.hpp', 'return r < w || r > m - w;', 'return ::boost_mut::near(v, r, w, m);')
    replace(root, 'impl/exp_impl.hpp', 'return rest < w || rest > unit - w;', 'return ::boost_mut::near(d, rest, w, unit);')
    replace(root, 'impl/trig_reduce.hpp',
            'return (r == 0U && dm.remainder < w) || (r == 499U && dm.remainder >= unit - w);',
            'return ::boost_mut::near_trig(v, r, dm.remainder, w, unit);')
    site = len(sites) + 1
    replace(root, 'impl/polynomial_impl.hpp', 'amb = units >= h || r <= units || r + units + 1U >= h;',
            f'::boost_mut::set({site}); amb = ::boost_mut::poly_ambiguous(out.d, units, h, r);')
    sites.append((site, 'impl/polynomial_impl.hpp', 'poly_digits', 'ambiguous'))
    site = len(sites) + 1
    replace(root, 'impl/polynomial_impl.hpp', 'if (((w.w[0] > eps && w.w[0] < ~eps) || e10 > hi10) && e10 + 1 >= lo10)',
            f'::boost_mut::set({site}); if ((!::boost_mut::poly_near(dd, w.w[0], eps) || e10 > hi10) && e10 + 1 >= lo10)')
    sites.append((site, 'impl/polynomial_impl.hpp', 'poly_try', 'near'))
    # The marked calls of ellint_near_boundary go to a new function that cuts the window of the check
    replace(root, 'impl/ellint_impl.hpp', 'constexpr auto ellint_near_boundary(', 'constexpr auto ellint_near_boundary_base(')
    replace(root, 'impl/ellint_impl.hpp', '\n// F or E for 0 < k <= 1 and phi > 0 in N words',
            '\nconstexpr auto ellint_near_boundary(const boost::int128::uint128_t& v, const int p, const unsigned bits) '
            'noexcept -> bool\n{\n    return ::boost_mut::near_bits(v, bits, [&v, p](const unsigned b) '
            '{ return ellint_near_boundary_base(v, p, b); });\n}\n\n// F or E for 0 < k <= 1 and phi > 0 in N words')
    replace(root, 'impl/fixed_point.hpp', '#include', '#include <bound_runtime.hpp>\n#include')

    for check in CHECKS:
        if not any(row[3] == check for row in sites):
            sys.exit(f'{check} has no site')
    with open(os.path.join(output, 'sites.txt'), 'w') as f:
        for row in sites:
            f.write(' '.join(str(cell) for cell in row) + '\n')
    print(len(sites), 'sites')


if __name__ == '__main__':
    main()
