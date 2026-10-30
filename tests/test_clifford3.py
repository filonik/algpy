"""Clifford functional calculus, including confluent and degenerate spectra."""
from functools import lru_cache
import itertools
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from algpy.algebras.general import clifford3
from field_reference import FUNCTIONS, reference
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


SIGNATURES = list(itertools.product((-1, 0, 1), repeat=3))


@lru_cache(None)
def scalar_type():
    # Each factory call allocates a NumPy dtype. Reuse one and mutate its K.
    return clifford3(-1, -1, -1)


@pytest.fixture(scope="module")
def matrix_functions():
    pytest.importorskip("scipy.linalg")
    from scipy_reference import FUNCTIONS as references
    return dict(references)


def operation(name):
    return {"inv": "reciprocal", "asin": "arcsin", "acos": "arccos",
            "atan": "arctan", "asinh": "arcsinh", "acosh": "arccosh",
            "atanh": "arctanh"}.get(name, name)


def center(name):
    return 2 if name == "acosh" else 1.5 if name in {"log", "sqrt", "inv"} else 0.3


def multiply(k, a, b):
    """Independent basis-blade multiplication with inversion-count signs."""
    result = np.zeros(8)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            factor = (-1)**sum(bool(i & (1 << p)) and bool(j & (1 << q))
                               for p in range(3) for q in range(p))
            for p in range(3):
                if i & j & (1 << p):
                    factor *= k[p]
            result[i ^ j] += factor*x*y
    return result


@pytest.mark.parametrize("parameters", SIGNATURES, ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_all_signatures(matrix_functions, parameters, name, function):
    cls = scalar_type()
    for k in dict.fromkeys([parameters, tuple(x*s for x, s in zip(parameters, (.25, 4, 2.25)))]):
        cls.K = k
        rows = np.random.default_rng(85).uniform(-0.035, 0.035, (2, 8))
        rows[:, 0] = center(name)
        values = cls.from_float(rows)
        output = np.empty(4, dtype=cls)[::2]
        assert function(values[::-1], out=output) is output
        for value, actual in zip(values[::-1], output):
            expected = matrix_functions[operation(name)](value.lmat())
            np.testing.assert_allclose(actual.lmat(), expected, rtol=3e-10, atol=3e-12)
            np.testing.assert_allclose(cls.to_float(getattr(value, name)()),
                                       cls.to_float(actual), rtol=3e-13, atol=3e-14)


NILPOTENTS = [
    ((0, 0, 0), [0, 1, 0, 0, 0, 0, 0, 0], 2),
    ((0, 0, 0), [0, 1, 0, 0, 0, 0, 1, 0], 3),
    ((1, 1, 0), [0, 1, 0, 1, 1, 0, 0, 0], 4),
]


def permute(parameters, coefficients, permutation):
    k = np.empty(3)
    k[list(permutation)] = parameters
    row = np.zeros(8)
    for mask, value in enumerate(coefficients):
        bits = [permutation[p] for p in range(3) if mask & (1 << p)]
        sign = (-1)**sum(x > y for i, x in enumerate(bits) for y in bits[i+1:])
        row[sum(1 << p for p in bits)] = sign*value
    return tuple(k), row


@pytest.mark.parametrize("parameters,coefficients,order", NILPOTENTS, ids=["square", "cube", "fourth"])
@pytest.mark.parametrize("permutation", list(itertools.permutations(range(3))))
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nilpotent_taylor(parameters, coefficients, order, permutation, name, function):
    parameters, row = permute(parameters, coefficients, permutation)
    cls = scalar_type()
    cls.K = parameters
    n = row/4
    w = center(name)
    expected = np.zeros(8)
    expected[0] = function(w)
    power = n.copy()
    factorial = 1
    for j in range(1, order):
        factorial *= j
        expected += reference(function, complex(w), j).real/factorial*power
        assert np.any(power)
        power = multiply(parameters, power, n)
    np.testing.assert_array_equal(power, 0)
    row = n.copy()
    row[0] = w
    np.testing.assert_allclose(cls.to_float(function(cls(*row))), expected, rtol=3e-12, atol=3e-13)


@pytest.mark.parametrize("imaginary", [-0.4, 0.4])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_complex_jordan_block(imaginary, name, function):
    cls = scalar_type()
    cls.K = (1, 1, 1)
    n = np.array([0, .25, 0, .25, 0, 0, 0, 0])
    z = complex(center(name), imaginary)
    value = function(z)
    derivative = reference(function, z, 1)
    I = np.eye(8)[7]
    expected = derivative.real*n + derivative.imag*multiply(cls.K, I, n)
    expected[0], expected[7] = value.real, value.imag
    row = n.copy()
    row[0], row[7] = z.real, z.imag
    np.testing.assert_allclose(cls.to_float(function(cls(*row))), expected, rtol=3e-12, atol=3e-13)


@pytest.mark.parametrize("q", [-1e-3, -1e-4, -1e-5, -1e-6, -1e-8, -1e-12, 0,
                               1e-12, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nearly_nilpotent_dual_center(matrix_functions, q, name, function):
    cls = scalar_type()
    cls.K = (1, 1-q, 0)
    row = [center(name), 1, 0, 1, 1, 0, 0, .2]
    value = cls(*row)
    expected = matrix_functions[operation(name)](value.lmat())
    np.testing.assert_allclose(function(value).lmat(), expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("last", [0] + [sign*10.**exponent for sign in [-1, 1]
                         for exponent in [-24, -16, -12, -11, -10, -9, -8, -6, -4, -2]])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nearly_dual_center(matrix_functions, last, name, function):
    cls = scalar_type()
    cls.K = (1, 1, last)
    value = cls(center(name), .25, 0, .25, .25, 0, 0, .2)
    expected = matrix_functions[operation(name)](value.lmat())
    np.testing.assert_allclose(function(value).lmat(), expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("delta", [-1e-4, -1e-8, 1e-8, 1e-4])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nearly_complex_jordan_block(matrix_functions, delta, name, function):
    cls = scalar_type()
    cls.K = (1, 1-delta, 1)
    value = cls(center(name), .25, 0, .25, 0, 0, 0, .4)
    expected = matrix_functions[operation(name)](value.lmat())
    np.testing.assert_allclose(function(value).lmat(), expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("parameters,row", [((-1, -1, -1), [1, 0, 0, 0, 0, 0, 0, 1]),
                                           ((1, -1, 1), [.75, .25, 0, 0, 0, 0, .25, .75])])
def test_square_root_with_zero_scalar_channel(parameters, row):
    cls = scalar_type()
    cls.K = parameters
    value = cls(*row)
    result = np.sqrt(value)
    # The nonzero spectra are {2} and {1,2}, respectively. Interpolation
    # also fixes sqrt(0)=0 without SciPy's ill-conditioned Schur result at 0.
    expected = np.array(row)/np.sqrt(2) if parameters[0] == -1 else (
        (2-1/np.sqrt(2))*np.array(row) + (1/np.sqrt(2)-1)*multiply(parameters, row, row))
    np.testing.assert_allclose(cls.to_float(result), expected, atol=2e-13)
    np.testing.assert_allclose(cls.to_float(result*result), row, atol=2e-13)


@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0)])
@pytest.mark.parametrize("scale", [1e-12, 1e-60, 1e-160])
@pytest.mark.parametrize("name", ["expm1", "log1p", "sqrt"])
def test_tiny_components(parameters, scale, name):
    cls = scalar_type()
    cls.K = parameters
    row = np.array([scale, scale, 0, 0, 0, 0, 0, scale])
    row[0] = 1 if name == "sqrt" else scale
    result = cls.to_float(getattr(np, name)(cls(*row)))
    if name == "sqrt":
        delta = row.copy()
        delta[0] -= 1
        square = multiply(parameters, delta, delta)
        cube = multiply(parameters, square, delta)
        expected = delta/2 - square/8 + cube/16
        expected[0] += 1
    else:
        square = multiply(parameters, row, row)
        cube = multiply(parameters, square, row)
        expected = row + square/2 + cube/6 if name == "expm1" else row - square/2 + cube/3
    np.testing.assert_allclose(result, expected, rtol=5e-12, atol=scale*1e-14)


def test_real_continuation_and_invalid_domains():
    cls = scalar_type()
    cls.K = (1, 1, 0)
    expected = np.zeros(8)
    expected[3] = 2
    np.testing.assert_allclose(cls.to_float(np.sqrt(cls(-4))), expected, atol=2e-14)
    cls.K = (1, 0, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        assert np.isnan(cls.to_float(np.sqrt(cls(-4)))).all()
        assert np.isnan(cls.to_float(np.log(cls(-4)))).all()
        assert np.isnan(cls.to_float(np.sqrt(cls(0, 0, 1)))).all()
        assert np.isnan(cls.to_float(np.log(cls(0, 0, 1)))).all()
    np.testing.assert_array_equal(cls.to_float(np.sqrt(cls(0))), np.zeros(8))


def test_current_parameters_and_power_broadcasting():
    cls = scalar_type()
    values = cls.from_float(np.array([[1.2, .1, .2, -.05, .1, .03, -.02, .2]]))
    cls.K = (-1, -1, -1)
    before = cls.to_float(np.exp(values)).copy()
    cls.K = (1, 1, 0)
    after = np.exp(values)
    assert not np.allclose(before, cls.to_float(after))
    np.testing.assert_allclose(cls.to_float(after[0]), cls.to_float(values[0].exp()))
    exponents = np.array([[0], [1], [2], [3]])
    output = np.empty((4, 2), dtype=cls)[:, ::2]
    assert np.power(values, exponents, out=output) is output
    for row, exponent in zip(output, exponents[:, 0]):
        np.testing.assert_allclose(cls.to_float(row[0]), cls.to_float(values[0]**int(exponent)))


@pytest.mark.parametrize("parameters", SIGNATURES, ids=str)
def test_powers_and_mix(matrix_functions, parameters):
    cls = scalar_type()
    cls.K = parameters
    row = np.random.default_rng(86).uniform(-.025, .025, 8)
    row[0] = 1.5
    a = cls(*row)
    b = cls(*np.roll(row, 2))
    for exponent in [-3, -1, 0, .3, .5, 1, 2, 3, 6]:
        from scipy.linalg import fractional_matrix_power
        expected = fractional_matrix_power(a.lmat(), exponent)
        np.testing.assert_allclose((a**exponent).lmat(), expected, rtol=3e-10, atol=3e-12)
    from scipy.linalg import expm, logm, solve
    expected = expm(logm(a.lmat())@b.lmat())
    np.testing.assert_allclose((a**b).lmat(), expected, rtol=3e-10, atol=3e-12)
    row = np.random.default_rng(87).uniform(-.025, .025, 8)
    row[0] = 1.2
    b = cls(*row)
    expected = a.lmat()@fractional_matrix_power(solve(a.lmat(), b.lmat()), .3)
    np.testing.assert_allclose(a.mix(b, .3).lmat(), expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("parameters,coefficients,order", NILPOTENTS)
def test_singular_integer_powers(parameters, coefficients, order):
    cls = scalar_type()
    cls.K = parameters
    a = cls(*coefficients)
    expected = np.eye(8)[0]
    for exponent in range(7):
        np.testing.assert_allclose(cls.to_float(a.ipow(exponent)), expected, atol=1e-14)
        expected = multiply(parameters, expected, coefficients)
    for exponent in [-3, 0, 1, 2, 3, 6]:
        assert cls.to_float(cls(-2).ipow(exponent))[0] == pytest.approx((-2.)**exponent)
    assert cls.to_float(cls(0).ipow(0))[0] == 1


def test_custom_apply_callbacks(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    source = tmp_path/"apply.c"
    source.write_text('''
#include "src/algpy/algebras/general/clifford/clifford3.h"
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
static void check(vec3 k, clifford3 a) {
    clifford3 actual = clifford3_apply(k, a, f, d1, d2, d3, cf, cd1);
    clifford3 square = clifford3_mul(k, a, a);
    clifford3 fourth = clifford3_mul(k, square, square);
    clifford3 want = clifford3_ladd(vec8_add(clifford3_mul(k, fourth, a), vec8_scale(a, 2)), 4);
    for (int i = 0; i < 8; i++) assert(fabs(actual.data[i]-want.data[i]) < 2e-12);
}
int main(void) {
    for (int i = -1; i <= 1; i++) for (int j = -1; j <= 1; j++) for (int l = -1; l <= 1; l++) {
        vec3 k = {{i, j, l}};
        check(k, clifford3_make(.3, .1, -.2, .15, .25, -.3, .125, .2));
    }
    vec3 k = {{1, 1, 1}};
    check(k, clifford3_make(.3, .25, 0, .25, 0, 0, 0, .2));
    k.data[2] = 0;
    check(k, clifford3_make(.3, 1, 0, 1, 1, 0, 0, .2));
    for (int i = -1; i <= 1; i += 2) {
        k.data[1] = 1 + i*1e-2;
        check(k, clifford3_make(.3, 1, 0, 1, 1, 0, 0, .2));
    }
    k = (vec3){{1e-6, 0, 0}};
    clifford3 actual = clifford3_apply(k, clifford3_make(.3, 1, 0, 0, 0, 0, 0, 0),
        fast, fast1, fast2, fast3, cfast, cfast1);
    assert(fabs(actual.data[0] - sin(300)*cos(1)) < 1e-12);
    assert(fabs(actual.data[1] - 1000*cos(300)*sin(1)) < 1e-9);
    return 0;
}
''')
    binary = tmp_path/"apply"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-I", str(root), str(source), "-lm", "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)

ALGEBRA_NAME = 'clifford3'


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

@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)


@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power_scipy(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power_scipy(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", [(-1, -1, -1), (1, 1, 1), (0, 0, 0), (0, -2, 0.5)], ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)

def test_norm_parameter_type():
    api.check_general_norm_type_does_not_follow_parameters('clifford', 3)

@pytest.mark.parametrize("parameters", [(-2, -2, -2), (-2, -2, 0), (-2, -2, 0.5), (-2, 0, -2), (-2, 0, 0), (-2, 0, 0.5), (-2, 0.5, -2), (-2, 0.5, 0), (-2, 0.5, 0.5), (0, -2, -2), (0, -2, 0), (0, -2, 0.5), (0, 0, -2), (0, 0, 0), (0, 0, 0.5), (0, 0.5, -2), (0, 0.5, 0), (0, 0.5, 0.5), (0.5, -2, -2), (0.5, -2, 0), (0.5, -2, 0.5), (0.5, 0, -2), (0.5, 0, 0), (0.5, 0, 0.5), (0.5, 0.5, -2), (0.5, 0.5, 0), (0.5, 0.5, 0.5)], ids=str)
def test_matrix_parameters(parameters):
    matrices.check_general_parameters('clifford', parameters)

def test_matrix_current_parameters():
    matrices.check_matrices_read_current_parameters('clifford', 3)

def test_negative_norm():
    api.check_clifford_negative_norm(3, (1, 1, -1), 7)

def test_ipow_current_parameters():
    powers.check_current_parameters()

@pytest.mark.parametrize("parameters", [(-2, -0.5, -1.5), (2, 0.5, 1.5), (0, 0, 0), (0, -2, 0.5)], ids=str)
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
