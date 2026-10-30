"""Tests for the 16-component complex specialization."""
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from algpy.algebras import complex4
from complex_reference import FUNCTIONS, character_matrix, reconstruct
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.fixture(scope="module")
def algebra():
    return complex4, character_matrix(complex4)


@pytest.mark.parametrize("operation,method", FUNCTIONS)
def test_elementary_functions(algebra, operation, method):
    cls, matrix = algebra
    coefficients = list(np.random.default_rng(34).uniform(-0.2, 0.2, (12, cls.N)))
    for real, first_imag, outer_imag in [(2, 0, 0), (-2, 0, 0), (0, 2, 0),
                                         (0, 0, 2), (0, 0, 0)]:
        row = np.zeros(cls.N)
        row[0], row[1], row[cls.N//2] = real, first_imag, outer_imag
        coefficients.append(row)
    # A nonzero zero divisor, whose characters are either zero or two.
    row = np.zeros(cls.N)
    row[0] = row[cls.N//2+1] = 1
    coefficients.append(row)
    coefficients = np.array(coefficients)
    components = coefficients @ matrix.T
    if operation == "log":
        # log is undefined on a zero complex component.
        valid = np.all(components != 0, axis=-1)
        coefficients, components = coefficients[valid], components[valid]
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(components))
    values = cls.from_float(coefficients)
    storage = np.empty(2*len(values), dtype=cls)
    out = storage[::2]
    assert function(values[::-1], out=out) is out
    np.testing.assert_allclose(cls.to_float(out), expected[::-1], rtol=2e-12, atol=2e-14)
    for value, reference in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(getattr(value, method)()), reference,
                                   rtol=2e-12, atol=2e-14)
        np.testing.assert_allclose(cls.to_float(function(value)), reference,
                                   rtol=2e-12, atol=2e-14)


@pytest.mark.parametrize("forward,inverse", [
    ("sin", "arcsin"), ("cos", "arccos"), ("tan", "arctan"),
    ("sinh", "arcsinh"), ("cosh", "arccosh"), ("tanh", "arctanh"),
    ("exp", "log"), ("square", "sqrt"),
])
def test_inverse_functions(algebra, forward, inverse):
    cls, _ = algebra
    value = cls(0.3, *np.linspace(-0.02, 0.02, cls.N-1))
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value),
                               rtol=2e-12, atol=2e-14)


def test_trigonometric_identities(algebra):
    cls, _ = algebra
    value = cls(0.3, *np.linspace(-0.02, 0.02, cls.N-1))
    circular = np.square(np.sin(value)) + np.square(np.cos(value))
    hyperbolic = np.square(np.cosh(value)) - np.square(np.sinh(value))
    np.testing.assert_allclose(cls.to_float(circular), cls.to_float(cls.one()), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(hyperbolic), cls.to_float(cls.one()), atol=2e-14)


def test_power_and_mix(algebra):
    cls, matrix = algebra
    x = cls(1.5, *np.linspace(-0.02, 0.02, cls.N-1))
    y = cls(1.3, *np.linspace(0.01, -0.01, cls.N-1))
    xp, yp = cls.to_float(x) @ matrix.T, cls.to_float(y) @ matrix.T
    expected = reconstruct(matrix, np.power(xp, yp))
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=2e-12, atol=2e-14)
    values = np.array([x, y], dtype=cls)
    np.testing.assert_allclose(cls.to_float(np.power(values, y))[0], expected,
                               rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(x**2), cls.to_float(np.square(x)), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), reconstruct(matrix, np.power(2, xp)),
                               rtol=2e-12, atol=2e-14)
    for t in [0, 0.4, 1]:
        expected = reconstruct(matrix, xp*np.power(yp/xp, t))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected,
                                   rtol=2e-12, atol=2e-14)


def test_arithmetic_and_metadata(algebra):
    cls, matrix = algebra
    rng = np.random.default_rng(35)
    x = cls(1.5, *rng.uniform(-0.02, 0.02, cls.N-1))
    y = cls(1.2, *rng.uniform(-0.01, 0.01, cls.N-1))
    xp, yp = cls.to_float(x) @ matrix.T, cls.to_float(y) @ matrix.T
    for actual, expected in ((x*y, xp*yp), (x/y, xp/yp), (x.inv(), 1/xp),
                             (np.square(x), xp*xp)):
        np.testing.assert_allclose(cls.to_float(actual), reconstruct(matrix, expected),
                                   rtol=2e-12, atol=2e-14)
    assert x.norm() == pytest.approx(np.prod(np.abs(xp))**(1/len(xp)))
    # Preserve the original sign-bit ordering of all generator characters.
    eigenvalues = []
    order = cls.N.bit_length()-1
    for index in range(cls.N):
        positive = index & 1
        row = 0
        for generator in range(1, order):
            bit = (index >> generator) & 1
            row = 2*row + (bit if positive else 1-bit)
        eigenvalues.append(xp[row] if positive else xp[row].conjugate())
    np.testing.assert_allclose(x.eivals(), eigenvalues, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(x.minpoly(), np.polynomial.polynomial.polyfromroots(eigenvalues),
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))


def test_standalone_header(tmp_path):
    order = 4
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix, lower = f"complex{order}", f"complex{order-1}"
    dim, half = 2**order, 2**(order-1)
    source = tmp_path/"standalone.c"
    source.write_text(f'''
#include "src/algpy/algebras/multicomplex/{prefix}.h"
#include <assert.h>
#if defined(ALGPY_MULTICOMPLEX1_H) || defined(ALGPY_MULTICOMPLEX2_H) || \\
    defined(ALGPY_MULTICOMPLEX3_H) || defined(ALGPY_MULTICOMPLEX4_H)
#error Concrete complex specializations must not depend on general multicomplex headers
#endif
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[1] = 0.25;
    a.data[{half}] = 0.125;
    {lower} n2 = {prefix}_norm2(a);
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    for (int i = 0; i < {half}; i++) {{
        assert(fabs(product.data[i] - n2.data[i]) < 1e-12);
        assert(fabs(product.data[i+{half}]) < 1e-12);
    }}
    assert(fabs({prefix}_norm{dim}(a) - pow({prefix}_norm(a), {dim})) < 1e-9);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)

ALGEBRA_NAME = 'complex4'


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

@pytest.mark.parametrize("parameters", [None], ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [None], ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [None], ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [None], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", [None], ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", [None], ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)

def test_api_specialization_structure():
    api.check_specialization_structure('complex4', 16)

@pytest.mark.parametrize("parameters", [None], ids=str)
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
