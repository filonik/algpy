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
The real field additionally supplies orders five through seven for the dual
derivative helpers, degenerate five-generator Clifford algebras and `jet5`.
The complex field supplies order five for `complex_jet5`; its fifth
derivatives of inverse functions are a polynomial times a power of the
principal first derivative, so they share its branch cuts.
The C names are `{field}_{function}_deriv{N}`, for example `double_sin_deriv3` and
`complex_log_deriv4`. They return derivatives, without division by the
factorial used in Taylor coefficients. Evaluate them where the function is
differentiable; complex derivatives follow the field functions' principal
branches, with signed zeros selecting the side of a branch cut.

The source is organized into [fields](../../src/algpy/fields),
[algebras](../../src/algpy/algebras),
[lattices](../../src/algpy/lattices),
[vectorspaces](../../src/algpy/vectorspaces), and
[polynomial helpers](../../src/algpy/polynomials).
Algebra arithmetic headers are independent of Python and NumPy. The adjacent
`_*.c` files define the Python scalar types and NumPy integration.

The fixed specializations, such as `quaternion`, `dual3`, `pga2`, or `cga3`,
are independent, optimized implementations: their headers and wrappers never
include or import the general parameterized algebras. They may use other
fixed specializations, for example `dual1` as the centre of `pga2`, and the
shared fields, vector spaces, polynomial helpers and lattices. The general
algebras, in turn, may build on fixed specializations where that makes sense.
They also share implementations among themselves: `clifford1` aliases the
general `cayleydickson1`, and `clifford3` and `clifford5` use `clifford1` for
their centre. The standalone-header tests reject fixed headers that include
general ones.

## C API and Python bindings

Every algebra header exposes the following C functions with its own prefix,
for example `quaternion_exp(a)` or `clifford3_exp(k, a)`. For general algebras,
operations that depend on the parameters take `vecN k` first. Operations such
as addition, projections, and conjugation do not take parameters. Concrete
specializations omit `k`.

The C constructor `{algebra}_make(const double coeffs[ALGPY_..._N])` copies
one full coefficient array in basis order into the returned value. The input
may be const or a slice of a larger array and must contain at least `N`
doubles. The result owns its coefficients; the input can be changed or go out
of scope after the call. For example:

```c
const double coeffs[ALGPY_QUATERNION_N] = {1, 2, 3, 4};
quaternion a = quaternion_make(coeffs);
quaternion b = quaternion_make((const double[]){1, 2, 3, 4});
```

Python constructors still take optional positional or keyword coefficients,
so Python code can use `quaternion(*coeffs)`. The wrappers pass their parsed
coefficient arrays directly to the C constructors. `from_re` remains the C
entry point for embedding a real scalar.

| C suffix                                          | Python interface                                                                                                      |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `make`, `from_re`, `to_re`                        | Type constructor, `Type.from_re(x)`, `x.to_re()`                                                                      |
| `zero`, `one`, `basis`, `gens`, `coeffs`          | Class methods with the same names                                                                                     |
| `is_real`, `re`, `im`                             | Scalar methods; `np.real(x)` and `np.imag(x)` also work on scalars                                                    |
| `pos`, `neg`                                      | Unary `+`/`-`, `np.positive`, `np.negative`                                                                           |
| `add`, `sub`, `mul`, `div`, `pow`                 | Arithmetic operators and `np.add`, `np.subtract`, `np.multiply`, `np.divide`, `np.power`                              |
| `lpow`, `rpow`                                    | `a.lpow(b) = exp(b * log(a))`, `a.rpow(b) = exp(log(a) * b)`                                                          |
| `ipow`                                            | `x.ipow(n)` for a signed 64-bit integer exponent                                                                      |
| `equal`, `not_equal`                              | `==`, `!=`, `np.equal`, `np.not_equal`                                                                                |
| `inv`, `square`, `conjugate`                      | Scalar methods, `np.reciprocal`, `np.square`, `np.conjugate`                                                          |
| `grade0`–`gradeN`                                 | Read-only scalar view `x.grade[k]`, with `len(x.grade) == N.bit_length()`                                             |
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

The Clifford algebras with an odd number of generators have the
two-dimensional centre `span{1, w}`, with `w` the pseudoscalar and `w*w`
given by `clifford1_pseudoscalar_square(k) = k0`,
`clifford3_pseudoscalar_square(k) = -k0*k1*k2` and
`clifford5_pseudoscalar_square(k) = k0*k1*k2*k3*k4`. The centre is the
`clifford1` algebra with that parameter. `clifford1_center_minpoly(k, a,
clifford1 p[2])`, `clifford3_center_minpoly(k, a, clifford1 p[3])` and
`clifford5_center_minpoly(k, a, clifford1 p[5])` write the minimal
polynomial over the centre in ascending order, with `p[i]` standing for
`p[i][0] + p[i][1]*w` and a leading coefficient of one. For `clifford1`,
which is commutative, it is `t - a`; for `clifford3` it is the central
quadratic `a*a - (a + conj(a))*a + a*conj(a)`. The real `minpoly` is its
product with the conjugate polynomial under `w -> -w`.

The fixed specializations with an odd number of generators implement the
same polynomials with their own arithmetic, since they do not include the
general headers, and type the coefficients by their centre:
`pga2_center_minpoly(a, dual1 p[3])`, `cga1_center_minpoly(a, splitcomplex1
p[3])`, `pga4_center_minpoly(a, dual1 p[5])` and `cga3_center_minpoly(a,
complex1 p[5])`, for `w*w` equal to 0, 1, 0 and -1. `cga3` reuses its
complex quartic, and `pga4` builds the quartic over the dual numbers from the
central traces of `b`, `b^2` and `b^3`, as `clifford5` does.

Python exposes `a.center_minpoly()`, returning an array of the type
`Type.center`. For `clifford1` that is the type itself. For `clifford3` and
`clifford5` it is a `clifford1` type created on first use, whose `K` the
parent type keeps equal to `[w*w]`. For the fixed types it is `dual1`,
`splitcomplex1` or `complex1`, imported on first use; a small metatype
provides the class attribute, as for the general types, so the extensions
still load without the package.

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

`{algebra}_gens(result)` writes the `n` generators of a `2ⁿ`-dimensional
algebra in parameter order, at coefficient positions `1 << i`. Like `basis`,
it does not take `k`. Each implementation constructs only those elements.
`real_gens(real *result)` writes nothing and accepts `NULL`; Python's
`real.gens()` returns an empty array. Other Python `gens()` methods return
arrays of shape `(n,)` with their algebra's dtype.

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
functions. Algebra-specific helpers (such as the
`norm2`/`norm4`/`norm8`/`norm16`/`norm32` reductions) remain additional
functions with family-dependent result types; `norm` is the common entry point.

Every algebra and lattice header provides grade projections
`{algebra}_grade0`–`{algebra}_gradeN` for its `N = 2ⁿ` coefficients:
`gradeK` keeps the coefficients whose basis index has `K` generator bits and
zeroes the rest, so the projections of an element sum to it. The jets, with
one generator `e`, have `jetN_grade0`–`jetN_gradeN` instead, where `gradeK`
keeps the coefficient of `e^K`; the complex jets likewise keep the complex
coefficient of `e^K`, its two real components. They do not take
`k`. Headers whose type aliases another algebra, such as `clifford1` and
`cayleydickson1`, forward with macros. The wrappers expose them as the
read-only sequence view `x.grade`, shared through `src/algpy/_grades.h`:
`x.grade[k]` is the grade-`k` projection, negative indices count from the
top grade, and out-of-range indices raise `IndexError`. Algebra scalars are
immutable and hashable, so `x.grade[k] = y` raises `TypeError`.

The return type of `norm` is fixed by the C algebra specialization:

| Algebras                                                                             | C return type | Python return type |
| ------------------------------------------------------------------------------------ | ------------- | ------------------ |
| `real`, `complex1`–`complex5`, `quaternion`, `octonion`, `sedenion`, `trigintaduonion`, `dual1`–`dual5`, `jet1`–`jet5`, `complex_jet1`–`complex_jet5`, `pga2`–`pga4`, `cga3` | `double`      | Built-in `float`   |
| `splitcomplex1`–`splitcomplex5`, `pga1`, `cga1`, `cga2`, and all parameterized general families | `complex`     | Built-in `complex` |

Complex-valued norms always return Python `complex`, including `0j` and `1+0j`.
Changing a general algebra's `K` parameters does not change this return type.
Wrappers call the algebra's C `norm` function and convert its declared result
directly; no runtime value inspection or custom scalar construction is needed.

## Functional calculus and specialization dependencies

### Multicomplex and recursive specializations

The fixed specializations are `complex1`–`complex5`, `dual1`–`dual5`, and
`splitcomplex1`–`splitcomplex5`.
The complex specializations build on the shared
complex field and the preceding complex specialization, with their own
arithmetic and recursive function evaluation. The split-complex specializations
likewise build on each other, using the real components `x+y` and `x-y` at the
base level. Their elementary functions and powers follow the domains of the
corresponding real functions in every component; values outside those domains
produce NaNs.

`complex5` and `splitcomplex5` have 32 real coefficients and recurse through
their fourth-level specializations. The complex `apply` helper evaluates
the two lower-dimensional values `f(u+i0*v)` and `f(u-i0*v)`, then reconstructs
the outer complex coefficient. The split-complex helper evaluates `f(u+v)`
and `f(u-v)` and reconstructs their half-sum and half-difference. Both take
one callback over the preceding specialization, with no extra derivatives
or polynomial solver. Their multiplication, inverses, and norm reductions
follow the same doubling as the earlier levels. `norm32` is real in both
specializations; `norm` returns `double` for `complex5` and `complex` for
`splitcomplex5`, whose determinant can be negative. Polynomial coefficients
are the convolution of the two degree-16 polynomials, and eigenvalues retain
the existing generator sign-bit ordering. Neither specialization depends on
the general multicomplex implementation.

The dual specializations also build on each other, with commuting generators
whose squares are zero. For `a = u + v*e`, function evaluation recurses as
`f(a) = f(u) + v*f'(u)*e`. The C helpers provide derivatives through order five
in `dual1`, four in `dual2`, three in `dual3`, two in `dual4`, and one in
`dual5`. All five specializations have their own arithmetic and metadata,
without general multicomplex dependencies. Callbacks must be differentiable at
the real coefficient. General powers use `exp(log(a)*b)` and require a positive
real coefficient in the base; reciprocals require a nonzero real coefficient.

`dual5` has 32 real coefficients and uses the existing fifth real field
derivatives. Its `apply` callbacks take and return `dual4`, following the
same structure as the earlier levels. Multiplication is the square-zero
doubling `(u+v*e)*(s+t*e) = u*s + (u*t+v*s)*e`; inversion and the norm
reductions likewise recurse through the lower specializations. The norm is
`fabs(a[0])`, with a static `double`/Python `float` return type. Polynomial
metadata retains the multicomplex-family convention `(x-a[0])^32` and 32
repeated eigenvalues, while elementary functions need derivatives only
through order five.

The general `multicomplex1`–`multicomplex5` implementations support mixed and
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
limits as `dual1`–`dual4`: four, three, two, and one respectively. These
additional derivatives support arbitrary mixtures of positive, negative, and
zero parameters in `multicomplex5`. Its 32 coefficients use the same bit-mask
basis order, recursive norm reductions through `norm32`, and degree-32 generic
polynomial metadata as the fifth-level fixed specializations. Its norm has
the general family's static `complex` C return type and Python `complex`
return type for every parameter set.

`multicomplex5_minpoly` multiplies factors from the explicit spectrum,
combining conjugate roots into real quadratics when a negative parameter is
present. This avoids the cancellation of the higher-degree trace recurrence
and retains all repeated roots when parameters are zero.

### Cayley-Dickson

The quaternion, octonion, sedenion, and trigintaduonion implementations have
fixed-parameter arithmetic and no dependencies on the general Cayley-Dickson
headers. Octonion multiplication uses quaternion pairs; sedenion multiplication
uses octonion pairs; trigintaduonion multiplication uses sedenion pairs.
Their `{algebra}_apply(a, complex_callback)` helpers evaluate the callback
at `re(a) + i*|im(a)|` and lift the result along the input's imaginary direction.
Real inputs select the first imaginary basis vector for the principal branch.
General `cayleydickson2`–`cayleydickson5` call their corresponding specializations
when all parameters equal `-1`.

Every Cayley-Dickson element satisfies
`a² - 2*re(a)*a + norm2(a) = 0`. Writing `a = w + v` gives the real scalar
`q = v² = -norm2(v)`. General `cayleydickson2`–`cayleydickson5` provide
`{algebra}_apply(k, a, real_callback, real_deriv1, complex_callback)`: negative
`q` uses complex evaluation, positive `q` uses the real values at
`w ± sqrt(q)`, and zero `q` uses `f(w) + f'(w)*v`. The order-one API uses the
equivalent multicomplex implementation. These evaluations take place in the
associative subalgebra generated by the input, including for octonions,
sedenions, and trigintaduonions.

`trigintaduonion` has 32 real coefficients and five generators at basis indices
1, 2, 4, 8, and 16. Its generic polynomial remains quadratic, and `eivals`
returns the two roots of that quadratic, as with the lower specializations.
`norm2` is the sum of the 32 coefficient squares, while `norm` has static
`double`/Python `float` type. Elementary functions use a balanced `hypot`
reduction of the imaginary direction to preserve very large and very small
components. These formulas do not require alternativity or norm multiplicativity;
the full algebra has zero divisors.

`cayleydickson5` doubles two `cayleydickson4` elements with an arbitrary fifth
generator square. The 32-term `norm2` expression weights each coefficient
square by the product of `-k[i]` over its basis bits. The generic polynomial
remains quadratic, so its `apply` helper needs only the real value and first
derivative callbacks and the complex value callback used at earlier levels.
All-negative unit parameters delegate multiplication, inverses, squares, and
function evaluation to `trigintaduonion`. The general type's norm always has
static `complex` C and Python return types, including for those parameters.

For nonreal inputs, positive `q` requires both real eigenvalues to lie in the
real function's domain; zero `q` requires differentiability at `w`. Invalid
evaluations return NaN components. Real inputs can continue into the complex
plane using the first negative-square basis vector, when one exists. Binary
powers default to `rpow`, using `exp(log(a)*b)`; `lpow` uses
`exp(b*log(a))`. Orders two through five evaluate real exponents by quadratic
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

`clifford5` has 32 coefficients and five arbitrary real parameters. Its
centre is `span{1,I}`, with `I = e0*e1*e2*e3*e4` and `I² = product(K)`.
Central traces determine a quartic over this centre. The public generic real
minimal polynomial remains degree eight, including multiplicities, and its
constant coefficient gives `norm8`. Like every parameterized type, `norm`
always has the static complex return type.

`clifford5_apply` dispatches explicitly by centre type. Negative `I²` uses
the complex quartic/Hermite engine. Positive `I²` projects `I` to its two real
values and evaluates in two four-generator algebras. Close split channels
use real Hermite interpolation to retain small coefficients. Zero `I²` uses
the roots of a real quartic, each twice, with `polynomials/hermite8.h`.
Scalar inputs, unused generators, and quadratic/cubic element relations have
separate paths. A dual-centre element annihilated by its quotient quartic
uses that quartic directly, preserving semisimple zero eigenvalues.

The degenerate case can contain an eighth-order nilpotent: `Cl(3,1,0)` is
`Mat(4,R)` over dual numbers. Its terminating Taylor series needs real
derivatives through order seven. `clifford5_apply` therefore takes `f`,
`deriv1` through `deriv7`, followed by `cf` and `cderiv1` through `cderiv3`.
The larger real Hermite helper uses these derivatives and the existing
complex derivative interpolation for nearly repeated nonreal roots.
For close native-function nodes, it samples the third complex derivative
on an analytic Cauchy contour to obtain stable higher divided differences.
These coefficients are shared across the interpolation table; sampling the
derivative preserves small changes when the function value has saturated.
It does not solve a degree-eight polynomial. All three general algebra families
support up to five generators.

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

### Fixed projective and conformal geometric algebras

The headers in `algebras/clifford/` provide `pga1`–`pga4`, `cga1`,
`cga2`, and `cga3`, independently of the general Clifford implementations. The suffix
counts geometry dimensions: PGA has one additional generator and CGA has two.
PGA's metric is `(0, 1, ...)`; CGA's is `(-1, 1, 1, ...)`. The first generator
indices are reserved for the null generator or the conformal pair. Basis
coefficients use the same bitmask order as the general Clifford algebras.

Multiplication and squaring use explicit expressions with fixed signs and
zero terms removed. `pga1` has a real quadratic generic minimal polynomial.
`pga2` evaluates a quadratic over its dual centre, while `cga1` uses its split
centre. Both have generic degree four, including multiplicities. Native
elementary functions use quadratic derivative formulas; custom callbacks use
the shared checked Hermite interpolation near repeated roots.

`pga3` and `cga2` use the quartic/Hermite engine and their own multiplication.
Scalar and lower-degree element relations have separate paths. When a regular
generator is unused, they can evaluate in the smaller fixed `pga2` or `cga1`
subalgebra. Dependencies never run from a fixed specialization to the general
Clifford implementation. Through `pga3` and `cga2`, the `apply` functions take
the same real and complex callbacks as `clifford3_apply`, without a metric
argument.

The null radical does not affect PGA's spectrum or norm. `pga2` and `pga3`
have nonnegative generic-polynomial constant terms and return real norms;
`pga1` can have a negative quadratic constant term. Its norm, and the `cga1`
and `cga2` norms, therefore have the static complex return type. Python wrappers use
these declared types directly.

`pga4` has 32 real coefficients and retains a generic degree-eight
polynomial, including repeated roots. Its Euclidean quotient has a real
quartic, whose square annihilates the full algebra. Null-generator
coefficients affect function values while leaving the polynomial,
eigenvalues, and norm unchanged. The norm has a static real return type.

Its inverse first inverts the Euclidean quotient using its quartic, then
restores the square-zero radical with `u_inv - u_inv*radical*u_inv`.
Its `apply` function uses the existing real Hermite engine through degree
eight, with real derivatives through order seven and complex derivatives
through order three. Scalars, unused regular generators, and quadratic/cubic
element relations have separate paths; a quotient quartic is used directly
when it also annihilates the element. All arithmetic remains independent of
the general Clifford family.

`cga3` is the fixed 32-component specialization of `clifford5(-1,1,1,1,1)`.
Its pseudoscalar `I = e0*e1*e2*e3*e4` is central and squares
to `-1`, identifying the centre with the complex field. The algebra is
isomorphic to 4x4 complex matrices. Scalar and pseudoscalar parts give
normalized complex traces, so its generic polynomial over the centre is
quartic. Explicit multiplication and central coefficient contractions compute
these traces without constructing a matrix. Quadratic/cubic element relations
and unused-generator subalgebras have separate paths.

`polynomials/poly4_complex.h` supplies closed-form complex quadratic, cubic,
and quartic helpers, with power-of-two scaling for the quartic. The separate
`polynomials/hermite_complex.h` evaluates cubic interpolants with arbitrary
complex nodes and multiplicities up to four. It mirrors the real engine's
checked derivative interpolation near coalescing nodes. Existing real
polynomial and Hermite interfaces are unchanged.

Unlike the smaller Clifford algebras, CGA3 can have a fourth-order Jordan
block at a nonreal eigenvalue. Its `cga3_apply` therefore takes the existing
real callbacks through `deriv3`, followed by `cf`, `cderiv1`, `cderiv2`, and
`cderiv3`. The complex field already provides all these derivatives. No
derivatives beyond order three or general-purpose degree-eight solver are
needed.

The public `cga3_minpoly` keeps the library's real, generic convention: it
returns the degree-eight product of the complex quartic with its coefficient
conjugate. `cga3_eivals` includes the four complex-quartic roots and their
conjugates. `norm8` is the squared absolute value of the complex determinant;
the norm is real and nonnegative, with a static `double`/Python `float` type.
The NumPy dtype has 32 coefficients and five generators, while shared vector
operations and coefficient conversions retain their existing interfaces.

### Jets

`src/algpy/algebras/jet/jet1.h`–`jet5.h` implement `R[e]/(e^(N+1))` on the
shared `vec2`–`vec6` types (`vec6` was added for `jet5`), coefficient `k`
holding `e^k`. They are fixed specializations, independent of the dual and
multicomplex algebras, and provide the common C API. Their one generator
sits at coefficient position 1, so `jetN_gens` writes a single element.

Multiplication is the truncated Cauchy product, expanded explicitly.
`jetN_inv` uses the series recurrence `y0 = 1/x0`,
`yk = -(x1*y(k-1) + ... + xk*y0)/x0`. Elementary functions use
`jetN_apply(a, f, derivs)`: with `a = x + n`, `n` nilpotent, it evaluates
`f(x) + sum f^(k)(x)/k! * n^k` by Horner's rule, using the real field's
`double_{function}_derivK` helpers for `K = 1, ..., N`. A function applied to
`x + e` therefore yields the normalized derivatives `f^(k)(x)/k!`. Integer
powers use repeated squaring; `lpow` and `rpow` agree, as jets commute. The
conjugate is `e -> -e`, the generic minimal polynomial is `(t - x)^(N+1)`
with the real part as its only eigenvalue, and `norm` is `|x|`.

The tests compare every elementary function with a truncated-Taylor reference
built from contour-integral derivatives, check that the `k`-th coefficient
of `f(x + e)` times `k!` equals the top mixed coefficient of `f` in `dualk`,
and check the chain rule on composed functions.

`complex_jet1.h`–`complex_jet5.h` implement `C[e]/(e^(N+1))` on `vec4`,
`vec6`, `vec8`, `vec10` and `vec12` (`vec10` and `vec12` were added for
them). Complex coefficient `k` is stored as `data[2k] + i*data[2k+1]`, so
over the reals the basis is `1, i, e, i*e, e^2, ...` and the generators `i`
and `e` sit at positions 1 and 2. They follow the real jets: an explicit
complex truncated product, the reciprocal recurrence in complex arithmetic,
and `complex_jetN_apply` with the complex field's derivatives on the
principal branches. The conjugate is complex conjugation of every
coefficient, as for `complex1`. Over the reals the generic minimal
polynomial is `((t - z)(t - conj(z)))^(N+1)` for the leading coefficient
`z`, with eigenvalues `z` and `conj(z)` each `N + 1` times, and `norm` is
`|z|`. The tests use a complex contour-integral reference and check that
complex jets with real arguments reproduce the real jets.

## Exterior lattices

The headers under `src/algpy/lattices/clifford` define `exterior1`–`exterior5`
on the shared `vecN` types. Their coefficient order is the bitmask order of
the Clifford algebras. The C API has `make`, `zero`, `bottom`, `top`,
`from_re`, `to_re`, `is_real`, `grade0`–`gradeN`, `pos`, `neg`, `add`, `sub`,
`rmul(double, a)`, `lmul(a, double)`, `ldiv(a, double)`, `equal`, `not_equal`,
`join`, `meet`, `ldual`, `rdual`, `norm2`, `basis`, and `gens`. The odd
dimensions add `dual`, which equals both `ldual` and `rdual`.

`join` is the Clifford product with all generator squares zero. On basis
blades, `rdual(e_a)` and `ldual(e_a)` are the signed complements satisfying
`join(rdual(e_a), e_a) = top` and `join(e_a, ldual(e_a)) = top`. `meet` is
`rdual(join(ldual(a), ldual(b)))`, which also equals
`ldual(join(rdual(a), rdual(b)))`. All products are expanded explicitly.

The wrappers map `join` to `nb_or` and `np.bitwise_or`, `meet` to `nb_and`
and `np.bitwise_and`, and, in odd dimensions, `dual` to `nb_invert` and
`np.invert`. In even dimensions `nb_invert` raises `TypeError` instead of
falling back to NumPy. The ufuncs `ldual` and `rdual` live in the small
extension `algpy._ufuncs`, alongside the other shared helpers, which creates
them without loops. Each extension that supports them imports it and
registers its loops, so one ufunc dispatches on every such dtype;
`algpy.lattices` re-exports them. Future ufuncs shared across families
belong in the same module.

Multiplication and division register only loops with one `float64` operand,
so lattice values can be scaled but not multiplied. The tests compare `join`
with the general `clifford(0, ..., 0)` product and all operations with an
independent bitmask reference in `tests/lattice_reference.py`.

### Lattice API of the fixed geometric algebras

The fixed PGA and CGA headers, `pga1.h`–`pga4.h` and `cga1.h`–`cga3.h`, also
implement this lattice API: `bottom`, `top`, `join`, `meet`, `ldual`, `rdual`,
and `dual` where the generator count is odd (`pga2`, `pga4`, `cga1`, `cga3`).
The exterior and regressive products and the Poincaré duals do not depend on
the metric, so these functions forward to the exterior lattice with the same
number of generators (`exterior2`–`exterior5`), which shares the coefficient
order and the `vecN` storage. `bottom` equals `one`. The wrappers register the
same operators, NumPy loops and `ldual`/`rdual` ufunc loops as the exterior
types, while `*` remains the geometric product. PGA and CGA keep only their
algebraic norms (`norm`, and the C `norm2`/`norm4`/`norm8` reductions); the
exterior types' coefficient `norm2` is not added. All of them share the
`x.grade[k]` projections.

For CGA's nondegenerate metric, the regressive product agrees with the
metric-dual meet `((a*I⁻¹) | (b*I⁻¹))*I` up to a sign that depends only on
the grades of `a` and `b`; the tests check this for all grade combinations.

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
