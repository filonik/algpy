from algpy.algebras import complex2
import numpy as np
import pytest
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import numpy_checks
import binding_checks


FUNCTIONS = [
    ("cos", "cos"), ("sin", "sin"), ("tan", "tan"),
    ("arccos", "acos"), ("arcsin", "asin"), ("arctan", "atan"),
    ("cosh", "cosh"), ("sinh", "sinh"), ("tanh", "tanh"),
    ("arccosh", "acosh"), ("arcsinh", "asinh"), ("arctanh", "atanh"),
]


def component_reference(function, coefficients):
    w, x, y, z = np.asarray(coefficients).T
    p = function((w-z) + 1j*(x+y))
    q = function((w+z) + 1j*(x-y))
    return np.stack(((p.real+q.real)/2, (p.imag+q.imag)/2,
                     (p.imag-q.imag)/2, (q.real-p.real)/2), axis=-1)


@pytest.mark.parametrize("ufunc,method", FUNCTIONS)
def test_elementary_functions(ufunc, method):
    rng = np.random.default_rng(22)
    coefficients = np.vstack((
        rng.uniform(-2, 2, (16, 4)),
        # Real branch cuts and a nonzero zero divisor (p=0, q=2).
        [2, 0, 0, 0], [-2, 0, 0, 0], [0, 2, 0, 0], [0, -2, 0, 0],
        [1, 0, 0, 1], [0, 0, 0, 0],
    ))
    function = getattr(np, ufunc)
    expected = component_reference(function, coefficients)
    values = complex2.from_float(coefficients)
    # Exercise strided input and output in the registered NumPy loop.
    storage = np.empty(2*len(values), dtype=complex2)
    out = storage[::2]
    result = function(values[::-1], out=out)
    assert result is out
    np.testing.assert_allclose(complex2.to_float(out), expected[::-1],
                               rtol=2e-13, atol=2e-14)
    for value, reference in zip(values, expected):
        np.testing.assert_allclose(complex2.to_float(getattr(value, method)()),
                                   reference, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(complex2.to_float(function(value)),
                                   reference, rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize("forward,inverse", [
    ("sin", "arcsin"), ("cos", "arccos"), ("tan", "arctan"),
    ("sinh", "arcsinh"), ("cosh", "arccosh"), ("tanh", "arctanh"),
])
def test_inverse_functions(forward, inverse):
    value = complex2(0.3, 0.2, -0.1, 0.05)
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(complex2.to_float(result), complex2.to_float(value),
                               rtol=2e-13, atol=2e-14)


def test_trigonometric_identities():
    values = complex2.from_float([[0.3, 0.2, -0.1, 0.05], [1, 0, 0, 1]])
    expected = np.tile([1, 0, 0, 0], (2, 1))
    np.testing.assert_allclose(complex2.to_float(np.square(np.sin(values)) + np.square(np.cos(values))),
                               expected, atol=2e-14)
    np.testing.assert_allclose(complex2.to_float(np.square(np.cosh(values)) - np.square(np.sinh(values))),
                               expected, atol=2e-14)


ALGEBRA_NAME = 'complex2'


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
    api.check_specialization_structure('complex2', 4)


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
