"""Jet algebra jet5: R[e]/(e^6), normalized derivatives by forward AD."""
import math
import numpy as np
import pytest
from algpy import algebras
from algpy.algebras import jet5
from dual_reference import FUNCTIONS, derivative, real_points
from jet_reference import product, taylor
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import numpy_checks
import binding_checks

ORDER = 5


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_elementary_functions(name, function):
    cls = jet5
    coefficients = np.random.default_rng(61 + ORDER).uniform(-0.3, 0.3, (3, cls.N))
    coefficients[:, 0] = real_points(name)
    expected = np.array([taylor(function, row) for row in coefficients])
    values = cls.from_float(coefficients)
    storage = np.empty(2*len(values), dtype=cls)
    out = storage[::2]
    assert function(values[::-1], out=out) is out
    np.testing.assert_allclose(cls.to_float(out), expected[::-1], rtol=3e-10, atol=3e-10)
    for value, result in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(function(value)), result, rtol=3e-10, atol=3e-10)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_normalized_derivatives(name, function):
    """f(x + e) holds f^(k)(x)/k!; the dual algebras give f^(k)(x) unnormalized."""
    cls = jet5
    e = cls.gens()[0]
    for x in real_points(name):
        jet = cls.to_float(function(x + e))
        expected = [derivative(function, x, k)/math.factorial(k) for k in range(cls.N)]
        np.testing.assert_allclose(jet, expected, rtol=3e-10, atol=3e-10)
        for k in range(1, cls.N):
            dual = getattr(algebras, f"dual{k}")
            mixed = dual.to_float(function(x + sum(dual.gens())))[-1]
            np.testing.assert_allclose(jet[k]*math.factorial(k), mixed, rtol=1e-11, atol=1e-11)


def test_chain_rule():
    """Composition of jets is composition of truncated Taylor series."""
    cls = jet5
    e = cls.gens()[0]
    for x in (-0.4, 0.1, 0.7):
        composed = np.exp(np.sin(x + e)*np.cos(x + e)) / (2 + np.tanh(x + e))
        reference = lambda z: np.exp(np.sin(z)*np.cos(z))/(2 + np.tanh(z))
        expected = [derivative(reference, x, k)/math.factorial(k) for k in range(cls.N)]
        np.testing.assert_allclose(cls.to_float(composed), expected, rtol=1e-9, atol=1e-9)


def test_structure_and_conjugate():
    cls = jet5
    e = cls.gens()[0]
    assert e.ipow(cls.N) == cls.zero() and (cls.N == 1 or e.ipow(cls.N - 1) != cls.zero())
    rng = np.random.default_rng(7 + ORDER)
    x, y = rng.normal(size=(2, cls.N))
    x[0] = 1.5 + abs(x[0])  # well-conditioned reciprocal series
    a, b = cls(*x), cls(*y)
    np.testing.assert_allclose(cls.to_float(a*b), product(x, y), rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(cls.to_float(a*b), cls.to_float(b*a), rtol=1e-14, atol=1e-14)
    signs = np.array([(-1)**k for k in range(cls.N)])
    np.testing.assert_array_equal(cls.to_float(a.conjugate()), signs*x)
    np.testing.assert_allclose(cls.to_float((a*b).conjugate()),
                               cls.to_float(a.conjugate()*b.conjugate()), rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(cls.to_float(a*a.inv()), np.eye(cls.N)[0], atol=1e-12)
    np.testing.assert_allclose(a.minpoly(), np.polynomial.Polynomial.fromroots([x[0]]*cls.N).coef)
    np.testing.assert_array_equal(a.eivals(), [x[0]]*cls.N)
    assert a.norm() == abs(x[0])
    keywords = cls(x=1, **{f"x{k}": k + 1 for k in range(1, cls.N)})
    np.testing.assert_array_equal(cls.to_float(keywords), np.arange(1, cls.N + 1))


ALGEBRA_NAME = 'jet5'


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

@pytest.mark.parametrize("parameters", [None], ids=str)
def test_basic_gens(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    basic.check_gens(common_algebra)

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

def test_api_grade_projection(common_algebra):
    api.check_grade_projection(common_algebra)

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
    api.check_specialization_structure('jet5', 6)


def test_numpy_exp_log_sqrt_power(common_algebra):
    numpy_checks.check_exp_log_sqrt_power(common_algebra)

def test_numpy_pow_operator(common_algebra):
    numpy_checks.check_pow_operator(common_algebra)

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
