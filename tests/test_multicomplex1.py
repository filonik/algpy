from algpy.algebras.general import multicomplex1
import itertools
import math
import numpy as np
import pytest
from algpy import algebras
from field_reference import FUNCTIONS
from multicomplex_reference import (scalar_type, scaled, expected, coefficients,
                                    multiply, c_derivatives, check_standalone_header_and_apply)
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.fixture(params=list(itertools.product([-1, 0, 1], repeat=1)), ids=lambda signs: str(signs))
def algebra(request):
    signs = request.param
    cls = scalar_type(len(signs))
    return cls, signs


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_functions_for_every_signature(algebra, name, function):
    cls, signs = algebra
    for parameters in dict.fromkeys([signs, scaled(signs)]):
        cls.K = parameters
        rows = np.array([coefficients(parameters, name, seed) for seed in [61, 62]])
        wanted = np.array([expected(parameters, function, row) for row in rows])
        values = cls.from_float(rows)
        output = np.empty(2*len(values), dtype=cls)[::2]
        assert function(values[::-1], out=output) is output
        np.testing.assert_allclose(cls.to_float(output), wanted[::-1], rtol=3e-9, atol=3e-10)
        for value, result in zip(values, wanted):
            np.testing.assert_allclose(cls.to_float(getattr(value, name)()), result,
                                       rtol=3e-9, atol=3e-10)
            np.testing.assert_allclose(cls.to_float(function(value)), result, rtol=3e-9, atol=3e-10)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("family,parameter", [("complex", -1), ("splitcomplex", 1), ("dual", 0)])
def test_specialization_agreement(name, function, family, parameter):
    order = 1
    cls = scalar_type(order)
    cls.K = [parameter]*order
    concrete = getattr(algebras, f"{family}{order}")
    row = coefficients(tuple(cls.K), name)
    actual = cls.to_float(function(cls.from_float(row)))
    wanted = concrete.to_float(function(concrete.from_float(row)))
    np.testing.assert_allclose(actual, wanted, rtol=2e-13, atol=2e-14)


def test_power_and_mix(algebra):
    cls, signs = algebra
    parameters = scaled(signs)
    cls.K = parameters
    xp = coefficients(parameters, "log")
    yp = coefficients(parameters, "log", 62)
    yp[0] = 1.25
    x, y = cls.from_float(np.array([xp, yp]))
    logarithm = expected(parameters, np.log, xp)
    wanted = expected(parameters, np.exp, multiply(parameters, logarithm, yp))
    np.testing.assert_allclose(cls.to_float(x**y), wanted, rtol=3e-9, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], wanted,
                               rtol=3e-9, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(x**2), multiply(parameters, xp, xp), atol=3e-12)
    wanted = expected(parameters, np.exp, math.log(2)*xp)
    np.testing.assert_allclose(cls.to_float(2**x), wanted, rtol=3e-9, atol=3e-10)
    for t in [0, 0.4, 1]:
        delta = t*(expected(parameters, np.log, yp)-logarithm)
        wanted = multiply(parameters, xp, expected(parameters, np.exp, delta))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), wanted, rtol=3e-9, atol=3e-10)


def test_arithmetic_and_metadata(algebra):
    cls, signs = algebra
    parameters = scaled(signs)
    cls.K = parameters
    xp = coefficients(parameters, "log")
    yp = coefficients(parameters, "log", 62)
    x, y = cls.from_float(np.array([xp, yp]))
    for actual, wanted in [(x*y, multiply(parameters, xp, yp)),
                           (x/y, multiply(parameters, xp, expected(parameters, np.reciprocal, yp))),
                           (x.inv(), expected(parameters, np.reciprocal, xp))]:
        np.testing.assert_allclose(cls.to_float(actual), wanted, rtol=3e-9, atol=3e-10)
    eigenvalues = []
    for signs in itertools.product([-1, 1], repeat=len(parameters)):
        factors = [sign*np.sqrt(complex(k)) for sign, k in zip(signs, parameters)]
        eigenvalues.append(sum(xp[i]*np.prod([v for bit, v in enumerate(factors) if i & (1 << bit)])
                               for i in range(cls.N)))
    np.testing.assert_allclose(np.sort_complex(x.eivals()), np.sort_complex(eigenvalues),
                               rtol=2e-13, atol=2e-14)
    polynomial = np.polynomial.polynomial.polyfromroots(eigenvalues)
    np.testing.assert_allclose(x.minpoly(), polynomial, rtol=3e-10, atol=3e-10)
    assert type(x.norm()) is complex
    assert x.norm() == pytest.approx(complex(polynomial[0].real)**(1/cls.N))


def test_real_component_powers():
    order = 1
    cls = scalar_type(order)
    cls.K = [4, 0.25, 9, 2][:order]
    row = coefficients(tuple(cls.K), "log")
    row[0] = -2
    x = cls.from_float(row)
    np.testing.assert_allclose(cls.to_float(x**2), multiply(cls.K, row, row), atol=3e-13)
    divisor = cls(1, 0.5)
    np.testing.assert_array_equal(cls.to_float(divisor**0), cls.to_float(cls.one()))


@pytest.mark.parametrize("parameters", [(-4,)])
@pytest.mark.parametrize("name,function", [("log", np.log), ("sqrt", np.sqrt), ("asin", np.arcsin),
                                           ("acos", np.arccos), ("acosh", np.arccosh)])
def test_complex_continuation(parameters, name, function):
    cls = scalar_type(len(parameters))
    cls.K = parameters
    row = np.zeros(cls.N)
    row[0] = -2
    row[1 << next(i for i, k in enumerate(parameters) if k < 0)] = 0.4
    for bit, k in enumerate(parameters):
        if k == 0:
            row[1 << bit] = 0.1
    wanted = expected(parameters, function, row)
    np.testing.assert_allclose(cls.to_float(function(cls.from_float(row))), wanted,
                               rtol=3e-9, atol=3e-10)


@pytest.mark.parametrize("offset", [1, 2, 3])
@pytest.mark.parametrize("index,name,function", [(i, name, f) for i, (name, f) in enumerate(FUNCTIONS)],
                         ids=[name for name, _ in FUNCTIONS])
def test_derivative_helpers(c_derivatives, offset, index, name, function):
    order = 1
    for signs in itertools.product([-1, 0, 1], repeat=order):
        parameters = scaled(signs)
        row = coefficients(parameters, name)
        wanted = expected(parameters, function, row, offset)
        actual = c_derivatives(parameters, index, offset, row)
        np.testing.assert_allclose(actual, wanted, rtol=3e-9, atol=3e-10)


def test_standalone_header_and_apply(tmp_path):
    check_standalone_header_and_apply(1, tmp_path)

ALGEBRA_NAME = 'multicomplex1'


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

@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", [(-1,), (1,), (0,)], ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)

def test_norm_parameter_type():
    api.check_general_norm_type_does_not_follow_parameters('multicomplex', 1)

@pytest.mark.parametrize("parameters", [(-2,), (0,), (0.5,)], ids=str)
def test_matrix_parameters(parameters):
    matrices.check_general_parameters('multicomplex', parameters)

def test_matrix_current_parameters():
    matrices.check_matrices_read_current_parameters('multicomplex', 1)

@pytest.mark.parametrize("parameters", [(-2,), (2,), (0,)], ids=str)
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
