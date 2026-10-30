"""Quartic Clifford calculus, including defective and degenerate spectra."""
from functools import lru_cache
import itertools
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from algpy.algebras.general import clifford4
from field_reference import FUNCTIONS, reference
from field_reference import center, operation
from polynomial_reference import assert_roots
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


SIGNATURES = list(itertools.product((-1, 0, 1), repeat=4))


@lru_cache(None)
def scalar_type():
    # Factories allocate dtypes. Reuse the type and change its parameters.
    return clifford4(-1, -1, -1, -1)


@pytest.fixture(scope="module")
def matrix_functions():
    pytest.importorskip("scipy.linalg")
    from scipy_reference import FUNCTIONS as references
    return dict(references)


def multiply(k, a, b):
    """Independent blade products, with inversion-count signs."""
    result = np.zeros(16)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            factor = (-1)**sum(bool(i & (1 << p)) and bool(j & (1 << q))
                               for p in range(4) for q in range(p))
            for p in range(4):
                if i & j & (1 << p):
                    factor *= k[p]
            result[i ^ j] += factor*x*y
    return result


@pytest.mark.parametrize("parameters", SIGNATURES, ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_all_signatures(matrix_functions, parameters, name, function):
    cls = scalar_type()
    for k in dict.fromkeys([parameters, tuple(x*s for x, s in zip(parameters, (.25, 4, 2.25, .5)))]):
        cls.K = k
        rows = np.random.default_rng(94).uniform(-.025, .025, (2, 16))
        rows[:, 0] = center(name)
        values = cls.from_float(rows)
        output = np.empty(4, dtype=cls)[::2]
        assert function(values[::-1], out=output) is output
        for value, actual in zip(values[::-1], output):
            expected = matrix_functions[operation(name)](value.lmat())
            np.testing.assert_allclose(actual.lmat(), expected, rtol=3e-10, atol=3e-12)
            np.testing.assert_allclose(cls.to_float(getattr(value, name)()),
                                       cls.to_float(actual), rtol=3e-13, atol=3e-14)


@pytest.mark.parametrize("generator", range(4))
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_three_generator_subalgebras(matrix_functions, generator, name, function):
    cls = scalar_type()
    cls.K = (1, -2, 0, 3)
    row = np.random.default_rng(95).uniform(-.025, .025, 16)
    row[[bool(i & (1 << generator)) for i in range(16)]] = 0
    row[0] = center(name)
    value = cls(*row)
    np.testing.assert_allclose(function(value).lmat(),
                               matrix_functions[operation(name)](value.lmat()),
                               rtol=3e-10, atol=3e-12)


def conjugate(cls, value):
    # Similarity activates the fourth generator, exercising the quartic path
    # with a known Jordan structure rather than the subalgebra shortcut.
    s = cls(1, 0, 0, 0, 0, 0, 0, 0, .25)
    return s*value*s.ipow(-1)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_fourth_order_nilpotent(name, function):
    cls = scalar_type()
    cls.K = (1, 1, 0, 1)
    n = conjugate(cls, cls(0, .25, 0, .25, .25, 0, 0, .125))
    row = cls.to_float(n)
    w = center(name)
    expected = np.zeros(16)
    expected[0] = function(w)
    power = row.copy()
    factorial = 1
    for order in range(1, 4):
        factorial *= order
        assert np.any(power)
        expected += reference(function, complex(w), order).real/factorial*power
        power = multiply(cls.K, power, row)
    np.testing.assert_allclose(power, 0, atol=1e-16)
    np.testing.assert_allclose(cls.to_float(function(n + cls(w))), expected,
                               rtol=3e-11, atol=3e-12)


@pytest.mark.parametrize("delta", [0, -1e-2, -1e-4, -1e-6, -1e-12, 1e-12, 1e-6, 1e-4, 1e-2])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_complex_jordan_and_nearby_spectra(matrix_functions, delta, name, function):
    cls = scalar_type()
    cls.K = (1, 1-delta, 1, 1)
    value = conjugate(cls, cls(center(name), .25, 0, .25, 0, 0, 0, .4))
    np.testing.assert_allclose(function(value).lmat(),
                               matrix_functions[operation(name)](value.lmat()),
                               rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("delta", [0, -1e-2, -1e-4, -1e-6, -1e-12, 1e-12, 1e-6, 1e-4, 1e-2])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_four_nearly_coincident_roots(matrix_functions, delta, name, function):
    cls = scalar_type()
    cls.K = (1, 1-delta, 0, 1)
    value = conjugate(cls, cls(center(name), 1, 0, 1, 1, 0, 0, .2))
    np.testing.assert_allclose(function(value).lmat(),
                               matrix_functions[operation(name)](value.lmat()),
                               rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("three_eigenvalues", [False, True])
def test_square_root_of_semisimple_zero(three_eigenvalues):
    cls = scalar_type()
    cls.K = (1, -1, 1, 1)
    e, f = cls.basis()[1], cls.basis()[6]
    p = (cls(1) + e)*(cls(1) + f)/4
    q = (cls(1) - e)*(cls(1) + f)/4
    value = conjugate(cls, p + 2*q if three_eigenvalues else p)
    expected = conjugate(cls, p + np.sqrt(2)*q if three_eigenvalues else p)
    actual = np.sqrt(value)
    np.testing.assert_allclose(cls.to_float(actual), cls.to_float(expected), rtol=3e-11, atol=3e-12)
    np.testing.assert_allclose(cls.to_float(actual*actual), cls.to_float(value), rtol=3e-11, atol=3e-12)


@pytest.mark.parametrize("parameters", [(-1,)*4, (1,)*4, (0,)*4])
@pytest.mark.parametrize("scale", [1e-12, 1e-60, 1e-160])
@pytest.mark.parametrize("name", ["expm1", "log1p", "sqrt"])
def test_tiny_components(parameters, scale, name):
    cls = scalar_type()
    cls.K = parameters
    row = np.random.default_rng(96).uniform(-scale, scale, 16)
    row[0] = 1 if name == "sqrt" else scale
    actual = cls.to_float(getattr(np, name)(cls(*row)))
    delta = row.copy()
    delta[0] -= 1 if name == "sqrt" else 0
    square = multiply(parameters, delta, delta)
    cube = multiply(parameters, square, delta)
    if name == "sqrt":
        expected = delta/2 - square/8 + cube/16
        expected[0] += 1
    else:
        expected = row + square/2 + cube/6 if name == "expm1" else row - square/2 + cube/3
    np.testing.assert_allclose(actual, expected, rtol=5e-12, atol=scale*1e-14)


def test_real_continuation_and_invalid_domains():
    cls = scalar_type()
    cls.K = (1, 1, 0, 1)
    expected = np.zeros(16)
    expected[3] = 2
    np.testing.assert_allclose(cls.to_float(np.sqrt(cls(-4))), expected, atol=2e-14)
    cls.K = (1, 0, 0, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        assert np.isnan(cls.to_float(np.sqrt(cls(-4)))).all()
        assert np.isnan(cls.to_float(np.log(cls(-4)))).all()
        n = cls(0, 0, 1)
        assert np.isnan(cls.to_float(np.sqrt(n))).all()
        assert np.isnan(cls.to_float(np.log(n))).all()
    np.testing.assert_array_equal(cls.to_float(np.sqrt(cls(0))), np.zeros(16))


@pytest.mark.parametrize("parameters", [(-1,)*4, (1,)*4, (0,)*4, (1, -2, 0, .5)])
def test_powers_and_mix(matrix_functions, parameters):
    from scipy.linalg import expm, logm, fractional_matrix_power, solve
    cls = scalar_type()
    cls.K = parameters
    row = np.random.default_rng(97).uniform(-.025, .025, 16)
    row[0] = 1.5
    a = cls(*row)
    b = cls(*np.roll(row, 2))
    for exponent in [-3, -1, 0, .3, .5, 1, 2, 3, 6]:
        np.testing.assert_allclose((a**exponent).lmat(),
            fractional_matrix_power(a.lmat(), exponent), rtol=3e-10, atol=3e-12)
    np.testing.assert_allclose((a**b).lmat(), expm(logm(a.lmat())@b.lmat()),
                               rtol=3e-10, atol=3e-12)
    b = a + cls(.3)
    expected = a.lmat()@fractional_matrix_power(solve(a.lmat(), b.lmat()), .3)
    np.testing.assert_allclose(a.mix(b, .3).lmat(), expected, rtol=3e-10, atol=3e-12)


def test_custom_apply_callbacks(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    source = tmp_path/"apply.c"
    source.write_text('''
#include "src/algpy/algebras/general/clifford/clifford4.h"
#include <assert.h>
static double f(double x) { return x*x*x*x*x + 2*x + 4; }
static double d1(double x) { return 5*x*x*x*x + 2; }
static double d2(double x) { return 20*x*x*x; }
static double d3(double x) { return 60*x*x; }
static complex cf(complex z) {
    return complex_ladd(complex_add(complex_mul(complex_square(complex_square(z)), z),
                                   complex_lmul(z, 2)), 4);
}
static complex cd1(complex z) {
    return complex_ladd(complex_lmul(complex_square(complex_square(z)), 5), 2);
}
static double fast(double x) { return sin(1000*x); }
static double fast1(double x) { return 1000*cos(1000*x); }
static double fast2(double x) { return -1e6*sin(1000*x); }
static double fast3(double x) { return -1e9*cos(1000*x); }
static complex cfast(complex z) { return complex_sin(vec2_scale(z, 1000)); }
static complex cfast1(complex z) { return vec2_scale(complex_cos(vec2_scale(z, 1000)), 1000); }
int main(void) {
    clifford4 a = {{.3, .1, -.2, .15, .25, -.3, .125, .2,
                   -.1, .05, .2, -.125, .15, -.2, .1, .125}};
    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++)
    for (int l = -1; l <= 1; l++) for (int m = -1; m <= 1; m++) {
        vec4 k = {{i, j, l, m}};
        clifford4 actual = clifford4_apply(k, a, f, d1, d2, d3, cf, cd1);
        clifford4 square = clifford4_mul(k, a, a);
        clifford4 fourth = clifford4_mul(k, square, square);
        clifford4 want = clifford4_ladd(vec16_add(clifford4_mul(k, fourth, a), vec16_scale(a, 2)), 4);
        for (int n = 0; n < 16; n++) assert(fabs(actual.data[n]-want.data[n]) < 3e-11);
    }
    /* Two commuting involutions give four nearby real roots. A custom
     * high-frequency callback must not use a native function's length scale. */
    vec4 k = {{1, -1, 1, 1}};
    clifford4 e = {{0}}, h = {{0}}, expected = {{0}};
    e.data[1] = 1; h.data[6] = 1;
    a = clifford4_ladd(vec16_add(vec16_scale(e, .001), vec16_scale(h, .002)), .3);
    for (int i = -1; i <= 1; i += 2) for (int j = -1; j <= 1; j += 2) {
        clifford4 projection = clifford4_mul(k, clifford4_ladd(vec16_scale(e, i), 1),
                                              clifford4_ladd(vec16_scale(h, j), 1));
        expected = vec16_add(expected, vec16_scale(projection, fast(.3+i*.001+j*.002)/4));
    }
    clifford4 s = clifford4_one(); s.data[8] = .25;
    clifford4 invs = clifford4_inv(k, s);
    a = clifford4_mul(k, clifford4_mul(k, s, a), invs);
    expected = clifford4_mul(k, clifford4_mul(k, s, expected), invs);
    clifford4 actual = clifford4_apply(k, a, fast, fast1, fast2, fast3, cfast, cfast1);
    for (int n = 0; n < 16; n++) assert(fabs(actual.data[n]-expected.data[n]) < 3e-11);
    /* A similarity transform gives a fourth-order nilpotent with small
     * multiplication roundoff in its nominally zero quartic coefficients. */
    k = (vec4){{1, 1, 0, 1}};
    a = (clifford4){{.3, .25, 0, .25, .25, 0, 0, .125}};
    invs = clifford4_inv(k, s);
    a = clifford4_mul(k, clifford4_mul(k, s, a), invs);
    clifford4 square = clifford4_mul(k, a, a);
    clifford4 fourth = clifford4_mul(k, square, square);
    expected = clifford4_ladd(vec16_add(clifford4_mul(k, fourth, a), vec16_scale(a, 2)), 4);
    actual = clifford4_apply(k, a, f, d1, d2, d3, cf, cd1);
    for (int n = 0; n < 16; n++) assert(fabs(actual.data[n]-expected.data[n]) < 3e-11);
    return 0;
}
''')
    binary = tmp_path/"apply"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-I", str(root), str(source), "-lm", "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)


@pytest.mark.parametrize("parameters", [
    (-1, -2, -3, -4), (1, 1, 1, 1), (-1, 0, 1, 0), (0, 0, 0, 0),
])
def test_clifford4_eigenvalues(parameters):
    from algpy.algebras.general import clifford4

    cls = clifford4(*parameters)
    # Scalar and generator inputs give repeated eigenvalues.
    assert_roots(cls.one().eivals(), [1.0]*4)
    generator_root = np.sqrt(complex(parameters[0]))
    assert_roots(cls.basis()[1].eivals(), [generator_root, -generator_root]*2)
    rng = np.random.default_rng(44)
    for _ in range(25):
        value = cls(*rng.normal(size=16))
        np.testing.assert_allclose(
            np.polynomial.polynomial.polyfromroots(value.eivals()),
            value.minpoly(), rtol=2e-10, atol=2e-10,
        )

ALGEBRA_NAME = 'clifford4'


def test_basic_one(common_algebra):
    basic.check_one(common_algebra)

def test_basic_zero(common_algebra):
    basic.check_zero(common_algebra)

def test_basic_add(common_algebra):
    basic.check_add(common_algebra)

def test_basic_sub(common_algebra):
    basic.check_sub(common_algebra)

def test_basic_mul(common_algebra):
    basic.check_mul(common_algebra)

def test_basic_inv(common_algebra):
    basic.check_inv(common_algebra)

def test_basic_div(common_algebra):
    basic.check_div(common_algebra)

def test_basic_norm(common_algebra):
    basic.check_norm(common_algebra)

def test_basic_minpoly_eivals(common_algebra):
    basic.check_minpoly_eivals(common_algebra)

def test_basic_basis(common_algebra):
    basic.check_basis(common_algebra)

def test_basic_coeffs(common_algebra):
    basic.check_coeffs(common_algebra)

def test_basic_float_conversion(common_algebra):
    basic.check_float_conversion(common_algebra)

def test_basic_k_mutable(common_algebra):
    basic.check_k_mutable(common_algebra)

@pytest.mark.parametrize("name,method", UNARY)
def test_api_unary_scalar_array_and_method(common_algebra, name, method):
    api.check_unary_scalar_array_and_method(common_algebra, name, method)

@pytest.mark.parametrize("name", ("add", "subtract", "multiply", "divide", "power"))
def test_api_binary_cast_broadcast_and_operator(common_algebra, name):
    api.check_binary_cast_broadcast_and_operator(common_algebra, name)

def test_api_predicates_parts_and_inverse(common_algebra):
    api.check_predicates_parts_and_inverse(common_algebra)

def test_api_norm_static_return_type(common_algebra):
    api.check_norm_static_return_type(common_algebra)

def test_conversion_native_methods_and_scalar_shape(common_algebra):
    conversions.check_native_methods_and_scalar_shape(common_algebra)

def test_conversion_shapes_and_contiguous_sharing(common_algebra):
    conversions.check_shapes_and_contiguous_sharing(common_algebra)

def test_conversion_strided_fortran_and_reversed_inputs(common_algebra):
    conversions.check_strided_fortran_and_reversed_inputs(common_algebra)

@pytest.mark.parametrize("dtype", [None, float, np.float32, np.int64, np.complex128, object,
                                  np.dtype(float).newbyteorder("S"), "U"])
def test_conversion_optional_output_dtype(common_algebra, dtype):
    conversions.check_optional_output_dtype(common_algebra, dtype)

def test_conversion_readonly_and_byte_swapped_inputs(common_algebra):
    conversions.check_readonly_and_byte_swapped_inputs(common_algebra)

def test_conversion_empty_inputs_and_invalid_shapes(common_algebra):
    conversions.check_empty_inputs_and_invalid_shapes(common_algebra)

def test_matrix_shape_dtype_and_action(common_algebra):
    matrices.check_matrix_shape_dtype_and_action(common_algebra)

def test_linearity_and_independent_results(common_algebra):
    matrices.check_linearity_and_independent_results(common_algebra)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)


@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power_scipy(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power_scipy(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)

def test_norm_parameter_type():
    api.check_general_norm_type_does_not_follow_parameters('clifford', 4)

@pytest.mark.parametrize("parameters", [(-2, -2, -2, -2), (-2, -2, -2, 0), (-2, -2, -2, 0.5), (-2, -2, 0, -2), (-2, -2, 0, 0), (-2, -2, 0, 0.5), (-2, -2, 0.5, -2), (-2, -2, 0.5, 0), (-2, -2, 0.5, 0.5), (-2, 0, -2, -2), (-2, 0, -2, 0), (-2, 0, -2, 0.5), (-2, 0, 0, -2), (-2, 0, 0, 0), (-2, 0, 0, 0.5), (-2, 0, 0.5, -2), (-2, 0, 0.5, 0), (-2, 0, 0.5, 0.5), (-2, 0.5, -2, -2), (-2, 0.5, -2, 0), (-2, 0.5, -2, 0.5), (-2, 0.5, 0, -2), (-2, 0.5, 0, 0), (-2, 0.5, 0, 0.5), (-2, 0.5, 0.5, -2), (-2, 0.5, 0.5, 0), (-2, 0.5, 0.5, 0.5), (0, -2, -2, -2), (0, -2, -2, 0), (0, -2, -2, 0.5), (0, -2, 0, -2), (0, -2, 0, 0), (0, -2, 0, 0.5), (0, -2, 0.5, -2), (0, -2, 0.5, 0), (0, -2, 0.5, 0.5), (0, 0, -2, -2), (0, 0, -2, 0), (0, 0, -2, 0.5), (0, 0, 0, -2), (0, 0, 0, 0), (0, 0, 0, 0.5), (0, 0, 0.5, -2), (0, 0, 0.5, 0), (0, 0, 0.5, 0.5), (0, 0.5, -2, -2), (0, 0.5, -2, 0), (0, 0.5, -2, 0.5), (0, 0.5, 0, -2), (0, 0.5, 0, 0), (0, 0.5, 0, 0.5), (0, 0.5, 0.5, -2), (0, 0.5, 0.5, 0), (0, 0.5, 0.5, 0.5), (0.5, -2, -2, -2), (0.5, -2, -2, 0), (0.5, -2, -2, 0.5), (0.5, -2, 0, -2), (0.5, -2, 0, 0), (0.5, -2, 0, 0.5), (0.5, -2, 0.5, -2), (0.5, -2, 0.5, 0), (0.5, -2, 0.5, 0.5), (0.5, 0, -2, -2), (0.5, 0, -2, 0), (0.5, 0, -2, 0.5), (0.5, 0, 0, -2), (0.5, 0, 0, 0), (0.5, 0, 0, 0.5), (0.5, 0, 0.5, -2), (0.5, 0, 0.5, 0), (0.5, 0, 0.5, 0.5), (0.5, 0.5, -2, -2), (0.5, 0.5, -2, 0), (0.5, 0.5, -2, 0.5), (0.5, 0.5, 0, -2), (0.5, 0.5, 0, 0), (0.5, 0.5, 0, 0.5), (0.5, 0.5, 0.5, -2), (0.5, 0.5, 0.5, 0), (0.5, 0.5, 0.5, 0.5)], ids=str)
def test_matrix_parameters(parameters):
    matrices.check_general_parameters('clifford', parameters)

def test_matrix_current_parameters():
    matrices.check_matrices_read_current_parameters('clifford', 4)

def test_negative_norm():
    api.check_clifford_negative_norm(4, (1, 1, 1, -1), 11)

@pytest.mark.parametrize("parameters", [(-2, -0.5, -1.5, -0.75), (2, 0.5, 1.5, 0.75), (0, 0, 0, 0), (0, -2, 0.5, -1.5), (-1, -2, -3, -4), (-1, -1, -1, -1), (-2, 0, 0.5, -1.5), (1, 1, 1, 1)], ids=str)
@pytest.mark.parametrize("operation", ['exp', 'expm1', 'log', 'log1p', 'reciprocal', 'square', 'sqrt', 'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh'])
@pytest.mark.parametrize("operator", ["lmat", "rmat"])
def test_scipy_elementary_functions(common_algebra, scipy_functions, parameters, operation, operator):
    from scipy_reference import sample_values
    if parameters is not None:
        common_algebra.K = parameters
    for value in sample_values(common_algebra, operation):
        result = getattr(np, operation)(value)
        expected = scipy_functions[operation](getattr(value, operator)())
        actual = getattr(result, operator)()
        assert np.isfinite(actual).all() and np.isfinite(expected).all()
        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-12, equal_nan=False)




def test_binding_nan_hash(common_algebra):
    binding_checks.check_nan_hash_stable(common_algebra)


def test_binding_subclass_arithmetic(common_algebra):
    binding_checks.check_subclass_mixed_arithmetic(common_algebra)


def test_binding_byteorder(common_algebra):
    binding_checks.check_dtype_byteorder(common_algebra)


def test_binding_real_input_dtypes(common_algebra):
    binding_checks.check_real_input_dtypes(common_algebra)


@pytest.mark.parametrize("reference", ['matrix'])
def test_benchmark_entry_point(tmp_path, reference):
    from benchmark_checks import check_benchmark_entry_point
    check_benchmark_entry_point(tmp_path, ALGEBRA_NAME, reference)
