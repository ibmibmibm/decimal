# Accuracy programs for the cmath functions

These programs check the results and the error bounds of the fixed-point code in `include/boost/decimal/detail/cmath`.
The tests and the CI do not build them.

A first try computes a function in few words and rounds the result only when an error bound keeps it away from a rounding boundary.
Else a more precise path computes the result again.
Thus a bound that is too small gives a wrong result, and a bound that is too large only costs time.
These programs show the measured error and the rounded results, so that you can compare them with each bound after a change.

## Build

The C++ programs need GCC or Clang with C++17, and `zeta_error`, `ellint_error` and `rounding_grid` also need MPFR and GMP.
`zeta_fit.py` needs Python 3 and mpmath.
Run these commands from the root of the repository:

```sh
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/rounding_grid.cpp -lmpfr -lgmp -o rounding_grid
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/zeta_error.cpp -lmpfr -lgmp -o zeta_error
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/ellint_error.cpp -lmpfr -lgmp -o ellint_error
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/erf_error.cpp -o erf_error
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/exp_error.cpp -o exp_error
```

The checks of the error bounds use a marked copy of the headers, in a directory that does not exist yet.
Build them with these commands:

```sh
python3 tools/cmath_accuracy/bound_mutate.py include marked
g++ -std=c++17 -O2 -Imarked -Itools/cmath_accuracy tools/cmath_accuracy/bound_mutants.cpp -o bound_mutants
g++ -std=c++17 -O2 -DBOUND_MARGIN -Imarked -Itools/cmath_accuracy tools/cmath_accuracy/rounding_grid.cpp -lmpfr -lgmp -o bound_margin
```

## Programs

`rounding_grid <type> <part>` compares the results with MPFR in all five rounding modes.
The type is `d32`, `d64`, `d128`, `f32`, `f64` or `f128`.
Each line gives a function, the number of checks and the number of wrong results.
The program returns 0 only if it finds no wrong result.
The parts are:

- `main`: a grid of short and long significands over a range of exponents.
- `exact`: arguments with exact results and ties, for `pow`, `hypot`, `exp2`, `log2`, `log10`, `tgamma`, `beta` and `atan2`.
- `near`: arguments next to the hard points of each function, such as multiples of pi/2, poles, zeros and overflow thresholds.
- `ext`: arguments and results next to the smallest and the largest values of the type.
- `ell`: `ellint_1`, `ellint_2` and `comp_ellint_2`.
- `assoc`: `assoc_laguerre` and `assoc_legendre`.
- `list <file>`: the inputs in a file, with a line `function|a|b|mode` for each input. The second argument is `0e0` for a function of one argument.

A result of zero must have the sign of the MPFR reference, except at a root of a polynomial, which has no sign.
The checks of `riemann_zeta(1)` always count as wrong in the parts `main` and `near`, because MPFR gives infinity and the library gives NaN.
The part `list` prints each input with its result, and the mark `wrong` if MPFR does not agree.
The part `list` gives MPFR the precision to hold a large argument exactly. Next to a rounding boundary, it doubles the precision up to 60000 bits, or takes the sign of the MPFR rounding error.
If `GRID_REC` names a file, the program writes each check to that file as `name|arguments|mode|result`.
For the fast types, the program does not check subnormal and overflow results.

`zeta_fit.py` prints the tables of `impl/riemann_zeta_impl.hpp` and the error of each interval.

`zeta_error <n>` prints the largest error of the first try of `riemann_zeta` for x > 0, against MPFR.

`ellint_error <n>` prints the largest error of the one-word first try of `comp_ellint_1` and `comp_ellint_2`, against an AGM (arithmetic-geometric mean) in MPFR.
It also prints the largest ratio of that error to the bound in the library.

`erf_error <count> <seed> <lo> <hi>` compares the two-word and three-word digits of `erf` and `erfc` for decimal128_t.

`exp_error <n>` compares the one-word and three-word digits of `exp` and `exp2`.

`bound_mutate.py <include> <new directory>` copies the headers and marks each check of an error bound as a site.
The file `sites.txt` in the new directory gives the number, the file and line, and the function of each site.

`bound_mutants <type> replay <record>...` reads the files that `GRID_REC` writes.
It runs each record again with each bound divided by 2, 4 and 16, and with the check removed, only where the site was next to a rounding boundary.
For each site, it prints the number of results that change.
The first line must show that no result differs without a cut, else the records do not match the headers.
MPFR checked the records, so a changed result is a wrong result.

`bound_mutants <type> hunt <range> <seconds> <seed>` runs random inputs in place of the records.
The ranges are in `bound_mutants.cpp`, and `beta_large` is the range for one argument above 10^15.
For each input and rounding mode that a cut changes, the program prints `site|cut|range|a|b|mode|new result|result`.
It prints the number of inputs for each site to the standard error stream.
Fields 3 to 6 of these lines are a file for the part `list` of `rounding_grid`, which checks the results against MPFR.

`bound_margin <type> <part>` is `rounding_grid` with a measure of each bound.
If `MARGIN_OUT` names a file, the program writes `site checks skipped ratio function arguments` to that file for each site.
The ratio is the largest error of the digits, divided by the bound of the site.
A ratio above 1 means that the bound is too small at those arguments.

## Limits

These programs measure errors at sample points, so they do not prove a bound.
A ratio below 1 from `bound_margin` is true only for the arguments of the grid, and a site that `bound_mutants` does not change can still need its bound.
The sweep of the polynomials for all n below 128 is not here.
