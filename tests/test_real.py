"""The real specialization follows ordinary double arithmetic directly."""
import numpy as np
import pytest
from algpy.algebras import real
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


ALGEBRA_NAME = "real"


@pytest.mark.parametrize("name", ["expm1", "log1p"])
@pytest.mark.parametrize("value", [-1e-300, -1e-20, -1e-12, 0, 1e-300, 1e-20, 1e-12])
def test_small_elementary_values(name, value):
    function = getattr(np, name)
    expected = function(value)
    np.testing.assert_allclose(real.to_float(function(real(value))), [expected], rtol=2e-15, atol=0)
    np.testing.assert_allclose(real.to_float(getattr(real(value), name)()), [expected], rtol=2e-15, atol=0)


@pytest.mark.parametrize("base,exponent", [(-2, 3), (-2, 2), (-2, -3), (-2, 0),
    (0, 0), (0, 1), (0, 2), (.25, -.5), (4, .5), (-2, .5)])
def test_power_matches_double(base, exponent):
    with np.errstate(all="ignore"):
        expected = np.power(float(base), float(exponent))
        scalar = real(base)**real(exponent)
        left = real(base)**float(exponent)
        right = float(base)**real(exponent)
        array = np.power(np.array([real(base)]), exponent)
    for actual in [scalar, left, right, array]:
        np.testing.assert_allclose(real.to_float(actual), expected, rtol=2e-15, atol=0)


@pytest.mark.parametrize("value", [-np.inf, -1e300, -1e-300, -0., 0., 1e-300, 1e300, np.inf])
def test_reciprocal_matches_double(value):
    with np.errstate(all="ignore"):
        expected = np.reciprocal(value)
        scalar = real(value).inv()
        array = np.reciprocal(np.array([real(value)]))
    for actual in [scalar, array]:
        np.testing.assert_allclose(real.to_float(actual), expected, rtol=2e-15, atol=0)
        assert np.signbit(real.to_float(actual)).all() == np.signbit(expected)


@pytest.mark.parametrize("numerator,denominator", [(1e-310, 1e-310), (1e300, 1e300),
    (1e-300, 1e300), (1e300, 1e-300), (2, -0.), (-2, -0.), (np.inf, 1), (1, np.inf)])
def test_division_matches_double(numerator, denominator):
    with np.errstate(all="ignore"):
        expected = np.divide(numerator, denominator)
        scalar = real(numerator)/real(denominator)
        left = real(numerator)/denominator
        right = numerator/real(denominator)
        array = np.divide(np.array([real(numerator)]), denominator)
    for actual in [scalar, left, right, array]:
        np.testing.assert_allclose(real.to_float(actual), expected, rtol=2e-15, atol=0)


@pytest.mark.parametrize("value", [-np.inf, -1e300, -1e-300, -0., 0., 1e-300, 1e300, np.inf])
def test_norm_matches_absolute_value(value):
    assert real(value).norm() == abs(value)


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
    api.check_specialization_structure('real', 1)

def test_small_matrices(common_algebra):
    matrices.check_small_commutative_matrices(common_algebra, [[2]])

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
