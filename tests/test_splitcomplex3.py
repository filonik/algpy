"""Tests for the 8-component split-complex specialization."""
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from algpy.algebras import splitcomplex3
from splitcomplex_reference import FUNCTIONS, character_matrix, reconstruct
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.fixture(scope="module")
def algebra():
    return splitcomplex3, character_matrix(splitcomplex3)


@pytest.mark.parametrize("operation,method,low,high", FUNCTIONS)
def test_elementary_functions(algebra, operation, method, low, high):
    cls, matrix = algebra
    components = np.random.default_rng(41).uniform(low, high, (12, cls.N))
    coefficients = reconstruct(matrix, components)
    # Compute the reference from the actual input coefficients after rounding.
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(coefficients @ matrix.T))
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


@pytest.mark.parametrize("operation", ["exp", "expm1", "sqrt", "sin", "cos", "tanh"])
def test_zero_divisors(algebra, operation):
    cls, matrix = algebra
    # A nonzero idempotent, with real characters zero and one.
    x = cls(0.5, 0.5)
    components = cls.to_float(x) @ matrix.T
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(components))
    np.testing.assert_allclose(cls.to_float(function(x)), expected, atol=2e-14)


@pytest.mark.parametrize("operation,value", [
    ("sqrt", -1), ("log", -1), ("log1p", -2),
    ("arcsin", 2), ("arccos", 2), ("arccosh", 0), ("arctanh", 2),
])
def test_real_function_domains(algebra, operation, value):
    cls, _ = algebra
    with np.errstate(all="ignore"):
        result = getattr(np, operation)(cls(value))
    assert np.isnan(cls.to_float(result)).all()


@pytest.mark.parametrize("forward,inverse,components", [
    ("sin", "arcsin", [-1, -0.5, 0.5, 1]),
    ("cos", "arccos", [-1, -0.5, 0.5, 1]),
    ("tan", "arctan", [-2, -0.5, 0.5, 2]),
    ("sinh", "arcsinh", [-2, -0.5, 0.5, 2]),
    ("cosh", "arccosh", [1, 1.5, 2, 3]),
    ("tanh", "arctanh", [-0.75, -0.5, 0.5, 0.75]),
    ("exp", "log", [0.25, 0.5, 1, 2]),
    ("square", "sqrt", [0, 0.5, 1, 2]),
])
def test_inverse_functions(algebra, forward, inverse, components):
    cls, matrix = algebra
    value = cls.from_float(reconstruct(matrix, np.resize(components, cls.N)))
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value), atol=2e-14)


def test_power_and_mix(algebra):
    cls, matrix = algebra
    rng = np.random.default_rng(42)
    xp = rng.uniform(0.5, 2, cls.N)
    yp = rng.uniform(0.5, 2, cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    expected = reconstruct(matrix, xp**yp)
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], expected,
                               rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), reconstruct(matrix, 2**xp), atol=2e-14)
    for t in [0, 0.4, 1]:
        expected = reconstruct(matrix, xp*(yp/xp)**t)
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected,
                                   rtol=2e-12, atol=2e-14)
    # Real pow handles negative components with integer exponents and zero**0.
    xp = np.resize([-2, -1, 0, 2], cls.N)
    yp = np.resize([2, 3, 0, 4], cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    np.testing.assert_allclose(cls.to_float(x**y), reconstruct(matrix, xp**yp), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(x**2), cls.to_float(np.square(x)), atol=2e-14)
    np.testing.assert_allclose(cls.to_float((-2)**y), reconstruct(matrix, (-2)**yp),
                               atol=2e-14)


def test_arithmetic_and_metadata(algebra):
    cls, matrix = algebra
    xp = np.random.default_rng(43).uniform(0.5, 2, cls.N)
    xp[0] = -xp[0]
    yp = np.random.default_rng(44).uniform(0.5, 2, cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    for actual, expected in ((x*y, xp*yp), (x/y, xp/yp), (x.inv(), 1/xp),
                             (np.square(x), xp*xp)):
        np.testing.assert_allclose(cls.to_float(actual), reconstruct(matrix, expected),
                                   rtol=2e-12, atol=2e-14)
    assert type(x.norm()) is complex
    assert x.norm() == pytest.approx(complex(np.prod(xp))**(1/cls.N))
    # The original four-component API uses its own character ordering.
    eigenvalues = xp[[3, 1, 2, 0]] if cls.N == 4 else xp
    np.testing.assert_allclose(x.eivals(), eigenvalues, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(x.minpoly(), np.polynomial.polynomial.polyfromroots(xp),
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))
    # Products of basis elements are indexed by the XOR of their bit masks.
    expected = np.eye(cls.N)[np.arange(cls.N)[:, None] ^ np.arange(cls.N)]
    np.testing.assert_array_equal(cls.to_float(cls.coeffs()), expected)


def test_standalone_header_and_apply(tmp_path):
    order = 3
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix, lower = f"splitcomplex{order}", f"splitcomplex{order-1}"
    half = 2**(order-1)
    callback = "square_real" if order == 1 else f"{lower}_square"
    norm_check = (f"assert(fabs(product.data[0] - {prefix}_norm2(a)) < 1e-12);"
                  if order == 1 else f'''
    {lower} n2 = {prefix}_norm2(a);
    for (int i = 0; i < {half}; i++) {{
        assert(fabs(product.data[i] - n2.data[i]) < 1e-12);
        assert(fabs(product.data[i+{half}]) < 1e-12);
    }}''')
    source = tmp_path/"standalone.c"
    source.write_text(f'''
#include "src/algpy/algebras/multicomplex/{prefix}.h"
#include <assert.h>
#if defined(ALGPY_MULTICOMPLEX1_H) || defined(ALGPY_MULTICOMPLEX2_H) || \\
    defined(ALGPY_MULTICOMPLEX3_H) || defined(ALGPY_MULTICOMPLEX4_H)
#error Concrete split-complex specializations must not depend on general multicomplex headers
#endif
static inline double square_real(double a) {{ return a*a; }}
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[{half}] = 0.25;
    {prefix} applied = {prefix}_apply(a, {callback});
    {prefix} squared = {prefix}_square(a);
    for (int i = 0; i < {2**order}; i++) {{
        assert(fabs(applied.data[i] - squared.data[i]) < 1e-12);
    }}
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    {norm_check}
    {prefix} two = {prefix}_from_re(2);
    {prefix} three = {prefix}_from_re(3);
    assert({prefix}_rpows(2, three).data[0] == 8);
    assert({prefix}_lpows(two, 3).data[0] == 8);
    assert(fabs({prefix}_mix(two, {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)

ALGEBRA_NAME = 'splitcomplex3'


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
    api.check_specialization_structure('splitcomplex3', 8)

def test_negative_norm():
    api.check_splitcomplex_negative_norm(3)

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
