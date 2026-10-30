# Implementation and development notes

For installation and the Python API, see the [user guide](../../README.md).
This document describes the C interfaces, numerical algorithms, bindings,
tests, and benchmarks.

## Design and shared fields

The implementation favors simple, explicit expressions and fast algebra
operations. Repetition is acceptable when it keeps the structure uniform and
bindings easy to generate mechanically. Internal numerical helpers should stay
focused on the needs of the algebra functions.

Shared complex arithmetic lives in `src/algpy/fields/complex.h`. The header at
`algebras/cayleydickson/complex.h` re-exports it for compatibility. Python scalar
types and NumPy wrappers live under `algebras`.

The fields `double` and `complex` provide ordinary derivatives of their
elementary unary functions, reciprocal (`inv`), and `square` through order four.
The C names are `{field}_{function}_deriv1` through
`{field}_{function}_deriv4`, for example `double_sin_deriv3` and
`complex_log_deriv4`. They return derivatives, without division by the
factorial used in Taylor coefficients. Evaluate them where the function is
differentiable; complex derivatives follow the field functions' principal
branches, with signed zeros selecting the side of a branch cut.

The source is organized into [fields](../../src/algpy/fields),
[algebras](../../src/algpy/algebras),
[vectorspaces](../../src/algpy/vectorspaces), and
[polynomial helpers](../../src/algpy/polynomials).
Algebra arithmetic headers are independent of Python and NumPy. The adjacent
`_*.c` files define the Python scalar types and NumPy integration.

## C API and Python bindings

Every algebra header exposes the following C functions with its own prefix,
for example `quaternion_exp(a)` or `clifford3_exp(k, a)`. For general algebras,
operations that depend on the parameters take `vecN k` first. Operations such
as addition, projections, and conjugation do not take parameters. Concrete
specializations omit `k`.

| C suffix                                          | Python interface                                                                                                      |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `make`, `from_re`, `to_re`                        | Type constructor, `Type.from_re(x)`, `x.to_re()`                                                                      |
| `zero`, `one`, `basis`, `coeffs`                  | Class methods with the same names                                                                                     |
| `is_real`, `re`, `im`                             | Scalar methods; `np.real(x)` and `np.imag(x)` also work on scalars                                                    |
| `pos`, `neg`                                      | Unary `+`/`-`, `np.positive`, `np.negative`                                                                           |
| `add`, `sub`, `mul`, `div`, `pow`                 | Arithmetic operators and `np.add`, `np.subtract`, `np.multiply`, `np.divide`, `np.power`                              |
| `lpow`, `rpow`                                    | `a.lpow(b) = exp(b * log(a))`, `a.rpow(b) = exp(log(a) * b)`                                                          |
| `ipow`                                            | `x.ipow(n)` for a signed 64-bit integer exponent                                                                      |
| `equal`, `not_equal`                              | `==`, `!=`, `np.equal`, `np.not_equal`                                                                                |
| `inv`, `square`, `conjugate`                      | Scalar methods, `np.reciprocal`, `np.square`, `np.conjugate`                                                          |
| `exp`, `expm1`, `log`, `log1p`, `sqrt`            | Scalar methods and NumPy ufuncs with the same names                                                                   |
| `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`       | Scalar methods and NumPy ufuncs with the same names                                                                   |
| `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh` | Scalar methods; NumPy uses `arcsin`, `arccos`, `arctan`, `arcsinh`, `arccosh`, `arctanh` (also scalar method aliases) |
| `norm`, `minpoly`, `eivals`                       | Scalar methods with the same names                                                                                    |
| `lmat`, `rmat`                                    | Instance methods returning `(N, N)` NumPy `float64` arrays                                                            |
| `mix`                                             | `x.mix(y, t)`, defined as `x * (inv(x) * y)**t` with those parentheses                                                |

Throughout the project, `minpoly` means the **generic minimal polynomial**
of the algebra, evaluated at the element. It can have repeated roots and can
have higher degree than the minimal polynomial of that particular element.
It can also have lower degree than the characteristic polynomial of a
multiplication matrix. The API keeps the concise `minpoly` name;
`generic_minpoly` would be the explicit spelling of the same concept.

Mixed real/algebra operations also have explicit C entry points: `radd`,
`ladd`, `rsub`, `lsub`, `rmul`, `lmul`, `rdiv`, `ldiv`, `rpows`, and `lpows`.
The `r` forms take `(double, algebra)` and the `l` forms take
`(algebra, double)`, preceded by `k` where needed. Python operators and
ufuncs accept real operands through the registered casts.

Every algebra exposes `{algebra}_lpow([k,] a, b)` for `exp(b * log(a))` and
`{algebra}_rpow([k,] a, b)` for `exp(log(a) * b)`, taking two algebra elements.
Python provides `a.lpow(b)` and `a.rpow(b)`, accepting an exponent from the
same algebra or a real number. Right-power remains the default for `a ** b`
and `np.power(a, b)`; the C function `{algebra}_pow` forwards to `rpow`.
In commutative algebras, `lpow` reuses `rpow` and its domain handling.
The internal mixed-scalar C functions are `lpows([k,] algebra, double)` and
`rpows([k,] double, algebra)`.

Every algebra provides `{algebra}_ipow([k,] a, int64_t n)`, using repeated
squaring with the algebra's multiplication. Negative exponents invert the base
once; exponent zero returns the identity, including `zero.ipow(0)`. Positive
powers therefore work for singular and nilpotent elements without a logarithm.
Python's `x.ipow(n)` accepts integer objects through `__index__`, including
NumPy integer scalars, and raises `OverflowError` outside the signed 64-bit
range. It preserves every exponent bit without converting to floating point.
Call `ipow` explicitly to select this operation. The general power functions,
`**` and `np.power`, use their algebra-specific formulas. The fields expose
`double_ipow` and `complex_ipow` as well.

Registered ufuncs support algebra scalars and arrays, broadcasting, strided
inputs, and `out=`. General types read their mutable `K` parameters at each
call. Component conversions are native class methods, `Type.to_float(a,
dtype=float)` and `Type.from_float(a)`, registered by the C wrappers along with
`Type.N` and `Type.dtype`. `to_float` appends a coefficient axis of length `N`;
`from_float` consumes that final axis, returning a zero-dimensional algebra
array for a single coefficient vector. Contiguous native-double buffers share
storage where possible; strided inputs and dtype conversions copy as needed.

Multiplication matrices use the same coefficient order as `Type.to_float`.
Left matrices act on coefficient columns and right matrices on coefficient rows:

```python
a.lmat() @ Type.to_float(b)  # Type.to_float(a * b)
Type.to_float(a) @ b.rmat()  # Type.to_float(a * b)
```

Both matrices are computed with `N` calls to the algebra's multiplication:
column `j` of `lmat(a)` contains the coefficients of `a*e_j`, and row `i` of
`rmat(a)` contains the coefficients of `e_i*a`. The same algorithm is repeated
in every implementation, using one temporary basis vector and one product
without constructing the structure tensor.
Each call returns an independent, writable array. Parameterized algebras read
their current `K` on each call. These operators describe a single multiplication,
including in nonassociative algebras. The C entry points are
`{algebra}_lmat([k,] a, result)` and `{algebra}_rmat([k,] a, result)`, filling
`N*N` doubles in row-major order.

The legacy NumPy dtype API does not dispatch array `.real`/`.imag` properties
to custom scalar methods. Use `Type.to_float(array)` for array component
access; the first component is the real coefficient. Scalar `re`/`im`
projections return values in the same algebra, with complementary zeroed
components.

All algebra families and their specializations provide the common elementary
functions. Algebra-specific helpers (such as grade
projections and `norm2`/`norm4`/`norm8`/`norm16` reductions) remain additional
functions with family-dependent result types; `norm` is the common entry point.

The return type of `norm` is fixed by the C algebra specialization:

| Algebras                                                                             | C return type | Python return type |
| ------------------------------------------------------------------------------------ | ------------- | ------------------ |
| `real`, `complex1`–`complex4`, `quaternion`, `octonion`, `sedenion`, `dual1`–`dual4` | `double`      | Built-in `float`   |
| `splitcomplex1`–`splitcomplex4` and all parameterized general families               | `complex`     | Built-in `complex` |

Complex-valued norms always return Python `complex`, including `0j` and `1+0j`.
Changing a general algebra's `K` parameters does not change this return type.
Wrappers call the algebra's C `norm` function and convert its declared result
directly; no runtime value inspection or custom scalar construction is needed.

## Functional calculus and specialization dependencies

### Multicomplex and recursive specializations

The fixed specializations are `complex1`–`complex4`, `dual1`–`dual4`, and
`splitcomplex1`–`splitcomplex4`.
The complex specializations build on the shared
complex field and the preceding complex specialization, with their own
arithmetic and recursive function evaluation. The split-complex specializations
likewise build on each other, using the real components `x+y` and `x-y` at the
base level. Their elementary functions and powers follow the domains of the
corresponding real functions in every component; values outside those domains
produce NaNs.

The dual specializations also build on each other, with commuting generators
whose squares are zero. For `a = u + v*e`, function evaluation recurses as
`f(a) = f(u) + v*f'(u)*e`. The C helpers provide derivatives through order three
in `dual1`, two in `dual2`, and one in `dual3`; `dual4` provides the value
functions. All four specializations have their own arithmetic and metadata,
without general multicomplex dependencies. Callbacks must be differentiable at
the real coefficient. General powers use `exp(log(a)*b)` and require a positive
real coefficient in the base; reciprocals require a nonzero real coefficient.

The general `multicomplex1`–`multicomplex4` implementations support mixed and
unequal real parameters. Positive generators split into real components, zero
generators carry derivatives, and negative generators supply complex
continuations. A negative outer generator is exchanged with an inner generator
when necessary to establish a complex direction. The standard all-negative,
all-positive, and all-zero parameter sets call the corresponding concrete
specializations directly. Dependencies run from the general algebras to these
specializations.

Elementary functions use the principal complex branches when a complex
direction is present, and the real function domains otherwise. Nilpotent
components require differentiability at each component's value. Powers split
over nonzero generators; zero generators use `exp(log(a)*b)`, with the same
domain requirements as the dual specializations. The C `multicomplex1_apply`
helper takes real value and derivative callbacks followed by a complex value
callback. Higher `apply` helpers take value and first-derivative callbacks in
the preceding general algebra. Derivative helpers follow the same order
limits as `dual1`–`dual3`.

### Cayley-Dickson

The quaternion, octonion, and sedenion implementations have fixed-parameter
arithmetic and no dependencies on the general Cayley-Dickson headers. Octonion
multiplication uses quaternion pairs; sedenion multiplication uses octonion
pairs. Their `{algebra}_apply(a, complex_callback)` helpers evaluate the callback
at `re(a) + i*|im(a)|` and lift the result along the input's imaginary direction.
Real inputs select the first imaginary basis vector for the principal branch.
General algebras with all parameters equal to `-1` use these specializations.

Every Cayley-Dickson element satisfies
`a² - 2*re(a)*a + norm2(a) = 0`. Writing `a = w + v` gives the real scalar
`q = v² = -norm2(v)`. General `cayleydickson2`–`cayleydickson4` provide
`{algebra}_apply(k, a, real_callback, real_deriv1, complex_callback)`: negative
`q` uses complex evaluation, positive `q` uses the real values at
`w ± sqrt(q)`, and zero `q` uses `f(w) + f'(w)*v`. The order-one API uses the
equivalent multicomplex implementation. These evaluations take place in the
associative subalgebra generated by the input, including for octonions and
sedenions.

For nonreal inputs, positive `q` requires both real eigenvalues to lie in the
real function's domain; zero `q` requires differentiability at `w`. Invalid
evaluations return NaN components. Real inputs can continue into the complex
plane using the first negative-square basis vector, when one exists. Binary
powers default to `rpow`, using `exp(log(a)*b)`; `lpow` uses
`exp(b*log(a))`. Orders two through four evaluate real exponents by quadratic
reduction, including `0**0 == 1`.
The `mix` formula is
`a * ((inv(a) * b)**t)`, preserving the displayed parentheses in nonassociative
algebras.

### Clifford and Hermite interpolation

`clifford3` evaluates elementary functions using a quadratic over the centre
`span{1, I}`, where `I = e0*e1*e2` and `I*I = -K[0]*K[1]*K[2]`.
The implementation treats complex, split, and dual centres explicitly, without
solving a quartic or constructing multiplication matrices. Its C helper is
`clifford3_apply(k, a, f, deriv1, deriv2, deriv3, cf, cderiv1)`: the real
callbacks provide the function and its first three ordinary derivatives, and
the complex callbacks provide the same analytic function and its first
derivative, with matching branch choices. Complex repeated roots use
`cderiv1`; fourth-order nilpotents use the terminating cubic Taylor expansion.
Nearly coincident roots and centres use derivative formulas to avoid
cancellation.
The C implementation keeps centre selection in `clifford3_apply`, with separate
helpers for the scalar/quadratic shortcuts, each centre type, and the numerical
corrections. Its quadratic interpolation and coalescing-root helpers live in
`polynomials/hermite_quadratic.h` and are shared with `clifford4`.

`clifford4_apply` takes the same callbacks. For `b = a - a[0]`, it computes
`b²` and `b³` and forms the centred generic quartic from their scalar parts;
the scalar part of `b⁴` is a coefficient contraction. The result is a real
cubic polynomial in `b`. Scalar, three-generator subalgebra, quadratic and
cubic element relations have separate paths, including nilpotents and
semisimple zero eigenvalues that do not require singular derivatives.

The algebra-independent `polynomials/hermite.h` provides
`poly4_apply(p, w, f, deriv1, deriv2, deriv3, cf, cderiv1, result)` for
ascending real quartic coefficients. It writes the four ascending coefficients
of the interpolant to `f(w+x)` and returns nonzero on success, or zero if the
function cannot be evaluated in the requested real polynomial calculus.
It uses the existing Ferrari solver,
with power-of-two variable scaling, followed by Newton/Hermite interpolation.
`hermite_interpolate` also accepts one to four supplied real/conjugate nodes.
Repeated nodes impose derivative conditions. Near native-function nodes use
derivative quadrature and interpolation without merging their separation;
custom callbacks use derivative interpolation when two quadrature or stencil
estimates agree, with ordinary differences as the fallback.

Real scalar inputs continue into the first negative-square basis blade when
needed, following the Cayley-Dickson convention. Nonreal inputs use the real
quadratic channels or the principal complex field functions; functions without
a finite value and the needed derivatives return NaN components. These branch
choices need not define a real primary matrix function on its branch cuts.
Use `ipow` for direct integer powers, including singular and nilpotent elements.

## Developer setup

Install an editable checkout with test dependencies:

```sh
python -m pip install -e ".[test]"
```

Rebuild the extensions after changing C sources or headers. For example,
with GCC:

```sh
CC=gcc CFLAGS='-O2 -Wall -Werror=return-type -Werror=incompatible-pointer-types' python setup.py build_ext --inplace
```

## Tests

Install test dependencies with `pip install -e '.[test]'`, then run `pytest`.
Each algebra has one `tests/test_{algebra}.py` file. It selects the common API,
conversion, multiplication-matrix, integer-power and SciPy checks and contains
the algebra's specific arithmetic and functional-calculus tests. For example,
`pytest tests/test_clifford4.py` runs the complete Clifford4 coverage.
Shared assertions, numerical references and fixtures live in support modules
without the `test_` prefix. `test_bindings.py` checks the complete C API and
extension loading; polynomial, Hermite, field and benchmark tests remain
separate library tests. No algebra test module imports another test module.

The SciPy comparisons check elementary functions through both multiplication
matrices for associative algebras, including split and nilpotent parameters.
For `cayleydickson4`, they compare the coefficient vector of `f(a)` with the
first column of `f(a.lmat())` and the first row of `f(a.rmat())`. Power
associativity makes the span of `1` and `a` invariant under these operators.
Inputs stay within the shared analytic domains. Separate polynomial and
Clifford tests cover repeated roots, nilpotents, nearly coincident spectra,
and all negative/positive/zero generator signatures. Matrix comparisons are
skipped if SciPy is unavailable.

## Benchmarks

Install the benchmark dependencies and compare the 16-component algebras with
SciPy matrix functions:

```sh
python -m pip install -e ".[benchmark]"
python scripts/benchmark_multicomplex4.py
python scripts/benchmark_cayleydickson4.py
python scripts/benchmark_clifford4.py
```

The scripts share a runner and measure nine functions across unequal negative,
unit negative, mixed, unit positive, and zero parameters. `multicomplex4` keeps
the original `complex`, `split`, and `dual` case names; the other scripts use
`unit-negative`, `unit-positive`, and `nilpotent`. Reports separate scalar
calls, batches, and calls including conversion to NumPy output. SciPy receives
precomputed input matrices. All outputs are checked for numerical agreement
before timing. BLAS uses one thread; timings include output allocations and
exclude input construction.

`multicomplex4` and `clifford4` compare complete 16×16 matrices and include a
timing for algebra calls followed by `lmat()`. `cayleydickson4` compares the
coefficient vector with the first column of SciPy's 16×16 result and includes
a timing for algebra calls followed by `Type.to_float()`. Its SciPy timings
include extracting that column. The recorded reference kind and matrix size
identify each comparison.

The Cayley-Dickson script also supports the smaller associative quadratic
slice `span{1,v}`, with `a=w+v` and `v²=q`. SciPy evaluates `[[w,q],[1,w]]`;
its first column gives the constant and `v` coefficients of `f(a)`. These
timings include reconstructing the 16 coefficients:

```sh
python scripts/benchmark_cayleydickson4.py --reference slice
```

Results are written to `benchmark-results/<algebra>/results.json` and
`results.csv`, with the Cayley-Dickson slice results in a `slice/` subdirectory.
JSON includes repeat samples, quartiles, environment metadata, and accuracy
checks; CSV summarizes timings and speedups. Older SciPy versions use a
recorded loop fallback when a matrix function lacks batch support.

Use `--cases`, `--functions`, `--samples`, `--repeats`, `--target`, and `--output`
to select parameter sets, functions, input count, repeat count, target seconds
per repeat, and output directory. For example:

```sh
python scripts/benchmark_multicomplex4.py --cases generic-negative --repeats 15 --target 0.05
python scripts/benchmark_clifford4.py --functions exp log sqrt --samples 64
```
