from algpy.algebras import complex1
import numpy as np
import pytest
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.mark.parametrize("operation", [
    "arccos", "arcsin", "arctan", "arccosh", "arcsinh", "arctanh",
])
def test_complex_inverse_principal_branches(operation):
    function = getattr(np, operation)
    # Both sides of real and imaginary branch cuts, preserving signed zeros.
    for re, im in [(2, 0.0), (2, -0.0), (-2, 0.0), (-2, -0.0),
                   (0.0, 2), (-0.0, 2), (0.0, -2), (-0.0, -2),
                   (0.3, 0.2), (-0.4, -0.6), (0, 0)]:
        expected = function(complex(re, im))
        actual = complex1.to_float(function(complex1(re, im)))
        np.testing.assert_allclose(actual, [expected.real, expected.imag], atol=2e-15)


def test_complex_minpoly_and_eivals():
    value = complex1(2, 3)
    assert value.minpoly().tolist() == [13, -4, 1]
    assert value.eivals().tolist() == [2 - 3j, 2 + 3j]


@pytest.mark.parametrize("scale", [1e-308, 1e-200, 1e-160, 1e160, 1e200, 1e308])
@pytest.mark.parametrize("real,imag", [(1, 0), (0, 1), (1, 1), (-1, 1)])
def test_complex_reciprocal_extreme_scales(scale, real, imag):
    # Compute the exact scaled reciprocal without forming scale**2. NumPy's
    # complex division can itself overflow for the largest balanced input.
    norm = real*real + imag*imag
    expected = np.array([(real/norm)/scale, (-imag/norm)/scale])
    actual = complex1.to_float(np.reciprocal(complex1(real*scale, imag*scale)))
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0)


def test_complex_reciprocal_subnormal_components():
    # Computing 1/x before dividing by two would overflow even though the
    # reciprocal's components are representable.
    scale = 5e-309
    expected = np.array([.5/scale, -.5/scale])
    actual = complex1.to_float(np.reciprocal(complex1(scale, scale)))
    np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0)


@pytest.mark.parametrize("scale", [1e-320, 1e-308, 1e200, 1e308])
def test_complex_division_equal_extreme_operands(scale):
    value = complex1(scale, scale)
    np.testing.assert_array_equal(complex1.to_float(value/value), [1, 0])


def test_complex_division_subnormal_components():
    smallest = np.nextafter(0., 1.)
    result = complex1.to_float(complex1(smallest)/complex1(smallest, 2*smallest))
    np.testing.assert_allclose(result, [.2, -.4], rtol=3e-15, atol=0)


def test_complex_division_large_components():
    a, b = complex1(1e308, .5e308), complex1(.9e308, .9e308)
    expected = complex(1, .5)/complex(.9, .9)
    np.testing.assert_allclose(complex1.to_float(a/b), [expected.real, expected.imag],
                               rtol=3e-15, atol=0)


@pytest.mark.parametrize("operation,real,imag", [
    ("tan", .25, 400), ("tan", .25, -400),
    ("tan", 0, 1000), ("tan", 0, -1000),
    ("tanh", 400, .25), ("tanh", -400, .25),
    ("tanh", 1000, 0), ("tanh", -1000, 0),
])
def test_complex_bounded_tangents_at_large_arguments(operation, real, imag):
    function = getattr(np, operation)
    expected = function(complex(real, imag))
    actual = complex1.to_float(function(complex1(real, imag)))
    assert np.isfinite(actual).all()
    np.testing.assert_allclose(actual, [expected.real, expected.imag], rtol=3e-15, atol=0)


ALGEBRA_NAME = 'complex1'


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
    api.check_specialization_structure('complex1', 2)

def test_conversion_dimension_limit():
    conversions.check_numpy_dimension_limit()

def test_small_matrices(common_algebra):
    matrices.check_small_commutative_matrices(common_algebra, [[2, -3], [3, 2]])

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
