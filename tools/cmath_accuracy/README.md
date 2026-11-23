# Accuracy programs for the cmath functions

These programs check the results and the error bounds of the fixed-point code in `include/boost/decimal/detail/cmath`.
The tests and the CI do not build them.

A first try computes a function in few words and rounds the result only when an error bound keeps it away from a rounding boundary.
Else a more precise path computes the result again.
Thus a bound that is too small gives a wrong result, and a bound that is too large only costs time.
These programs show the measured error and the rounded results, so that you can compare them with each bound after a change.

## Build

The C++ programs need GCC 9 or Clang 9, or a later version, with C++17, and `rounding_grid` also needs MPFR and GMP.
`near_ties.py` needs Python 3 and mpmath.
Run these commands from the root of the repository:

```sh
g++ -std=c++17 -O2 -Iinclude tools/cmath_accuracy/rounding_grid.cpp -lmpfr -lgmp -o rounding_grid
```

## Programs

`rounding_grid <type> <part>` compares the results with MPFR in all five rounding modes.
The type is `d32`, `d64`, `d128`, `f32`, `f64` or `f128`.
Each line gives a function, the number of checks and the number of wrong results.
The program returns 0 only if it finds no wrong result.
If `GRID_STRICT` is set, an inconclusive result also gives the status 1.
The parts are:

- `main`: a grid of short and long significands over a range of exponents.
- `short`: each 16th argument of the grids of `main`, as a fast check of the comparison and of the records.
- `exact`: arguments with exact results and ties, for `pow`, `hypot`, `exp2`, `log2`, `log10`, `tgamma`, `beta` and `atan2`.
- `near`: arguments next to the hard points of each function, such as multiples of pi/2, poles, zeros and overflow thresholds.
- `ext`: arguments and results next to the smallest and the largest values of the type.
- `ell`: `ellint_1`, `ellint_2` and `comp_ellint_2`.
- `assoc`: `assoc_laguerre` and `assoc_legendre`.
- `self`: wrong and right results for the checker itself, next to 0, the overflow threshold and a power of ten. For a fast type, it has only results that are out of the range.
- `list <file>`: the inputs in a file, with a line `function|a|b|mode` for each input. The second argument is `0e0` for a function of one argument.
A line counts as wrong if the type cannot hold an argument exactly, for example an argument with too many digits.

A result of zero must have the sign of the MPFR reference, except at a root of a polynomial, which has no sign.
At the pole `riemann_zeta(1)`, MPFR gives infinity, and a result counts as right only if it is NaN, which the library documents.
The part `list` prints each input with its result, and the mark `wrong` if MPFR does not agree.
The part `list` gives MPFR the precision to hold a large argument exactly. Next to a rounding boundary, it doubles the precision up to 60000 bits, or takes the sign of the MPFR rounding error.
For an argument that binary does not hold exactly, that sign comes from the ends of its binary interval, in the grids and in the part `list`. This is true only for a function that is monotonic there, so `sin`, `cos`, `tan`, `tgamma`, `lgamma` and `riemann_zeta` for x up to -2.7 take more bits.
An argument in a record and in the part `list` is a number, or `nan`, `inf` or `-inf`.
If `GRID_REC` names a file, the program writes each result that the comparison accepts to that file as `name|arguments|mode|result`.
If the file does not open or does not get all records, the exit status is not 0.
A wrong result, a result of a fast type that has no comparison, and the results of the part `self` get no record.
A result is inconclusive if the reference stays on a rounding boundary at the largest precision and the sign of its error is not known.
The result is wrong only if it is wrong for a value on the boundary and for values on each side of it, and the result gets no record.
The program prints the first 8 of these results, and the part `list` marks each of them `inconclusive`.
Integer arithmetic finds an exact value on a boundary for `log10`, `sqrt`, `cbrt` and `hypot`.
It also finds one for `pow` with `x` equal to 1 or an exponent `P/Q`, where `|P|` is at most 10000 and `Q` is at most 1000, and for `beta` with an integer argument from 1 to 64.
Then the result is not inconclusive.
A reference above two times the largest value of the type is not on a boundary, and only the overflow rules apply.
For an infinite result, the overflow threshold is a boundary of the largest value, so the precision of the reference goes up next to it as next to other boundaries.
For the fast types, an overflow result has the same check as for the other types.
A fast type has no subnormal results, so a result below the smallest normal value must be at most that value, and no rounding check follows.
The last line gives the number of inconclusive results, and the number of results of a fast type that have no rounding check.

`near_ties.py <digits> <count> <seed>` writes a file `near_<case>.txt` for each case, for the part `list` of `rounding_grid`.
The cases are the functions that are a constant plus a small term, for example `acos(x)` next to pi/2 for a small x.
There the arguments that are next to each other give results that are very close, so a result can be much nearer to a rounding boundary than a random argument gives.
Each line has an argument that the program solved for a rounding boundary and then rounded to the digits of the type.
The ranges of the arguments are those of 34 digits, and fewer digits take a part of each range in proportion.

## Limits

These programs measure errors at sample points, so they do not prove a bound.
The sweep of the polynomials for all n below 128 is not here.
