"""Complex jet algebra complex_jet5: C[e]/(e^6), normalized complex derivatives."""
import math
import numpy as np
import pytest
from algpy.algebras import complex_jet5, jet5
from dual_reference import FUNCTIONS
from jet_reference import (complex_derivative, complex_product, complex_taylor,
                           from_complex, to_complex)
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import numpy_checks
import binding_checks

ORDER = 5


def complex_points(name):
    # Away from the principal branch cuts and poles.
    if name in ("log", "sqrt", "inv", "acosh"):
        return [1.3 + 0.4j, 0.8 - 0.5j]
    return [0.3 + 0.2j, -0.2 - 0.35j]


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_elementary_functions(name, function):
    cls = complex_jet5
    rng = np.random.default_rng(71 + ORDER)
    coefficients = rng.uniform(-0.2, 0.2, (2, ORDER + 1)) + 1j*rng.uniform(-0.2, 0.2, (2, ORDER + 1))
    coefficients[:, 0] = complex_points(name)
    expected = np.array([complex_taylor(function, row) for row in coefficients])
    values = from_complex(cls, coefficients)
    out = np.empty_like(values)
    assert function(values, out=out) is out
    np.testing.assert_allclose(to_complex(cls, out), expected, rtol=2e-8, atol=2e-8)
    for value, result in zip(values, expected):
        np.testing.assert_allclose(to_complex(cls, function(value)), result, rtol=2e-8, atol=2e-8)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_normalized_derivatives_and_real_restriction(name, function):
    """f(z + e) holds f^(k)(z)/k!; on the real axis it is the real jet."""
    cls = complex_jet5
    i, e = cls.gens()
    for z in complex_points(name):
        jet = to_complex(cls, function(z.real + z.imag*i + e))
        expected = [complex_derivative(function, z, k)/math.factorial(k) for k in range(ORDER + 1)]
        np.testing.assert_allclose(jet, expected, rtol=2e-8, atol=2e-8)
        x = abs(z.real) + 1 if name in ("log", "sqrt", "inv", "acosh") else z.real
        real = jet5.to_float(function(x + jet5.gens()[0]))
        np.testing.assert_allclose(to_complex(cls, function(x + e)), real, rtol=1e-14, atol=1e-14)


def test_chain_rule():
    cls = complex_jet5
    i, e = cls.gens()
    reference = lambda z: np.exp(np.sin(z)*np.cos(z))/(2 + np.tanh(z))
    for z in (0.3 + 0.2j, -0.4 + 0.1j):
        composed = np.exp(np.sin(z.real + z.imag*i + e)*np.cos(z.real + z.imag*i + e)) / (2 + np.tanh(z.real + z.imag*i + e))
        expected = [complex_derivative(reference, z, k)/math.factorial(k) for k in range(ORDER + 1)]
        np.testing.assert_allclose(to_complex(cls, composed), expected, rtol=1e-8, atol=1e-8)


def test_structure_and_conjugate():
    cls = complex_jet5
    i, e = cls.gens()
    assert i*i == cls(-1) and e.ipow(ORDER + 1) == cls.zero()
    assert ORDER == 0 or e.ipow(ORDER) != cls.zero()
    rng = np.random.default_rng(9 + ORDER)
    x, y = rng.normal(size=(2, ORDER + 1)) + 1j*rng.normal(size=(2, ORDER + 1))
    x[0] = 1.5 + abs(x[0])  # well-conditioned reciprocal series
    a, b = from_complex(cls, x)[()], from_complex(cls, y)[()]
    np.testing.assert_allclose(to_complex(cls, a*b), complex_product(x, y), rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(cls.to_float(a*b), cls.to_float(b*a), rtol=1e-14, atol=1e-14)
    np.testing.assert_array_equal(to_complex(cls, a.conjugate()), np.conj(x))
    np.testing.assert_allclose(cls.to_float((a*b).conjugate()),
                               cls.to_float(a.conjugate()*b.conjugate()), rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(cls.to_float(a*a.inv()), np.eye(cls.N)[0], atol=1e-12)
    np.testing.assert_allclose(a.minpoly(), np.poly([x[0], np.conj(x[0])]*(ORDER + 1))[::-1].real,
                               rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(a.eivals(), [x[0], np.conj(x[0])]*(ORDER + 1))
    assert a.norm() == pytest.approx(abs(x[0]))
    keywords = cls(x=1, y=2, **{f"{p}{k}": 2*k + (p == "y") + 1 for k in range(1, ORDER + 1) for p in "xy"})
    np.testing.assert_array_equal(cls.to_float(keywords), np.arange(1, cls.N + 1))


ALGEBRA_NAME = 'complex_jet5'


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
    api.check_specialization_structure('complex_jet5', 12)


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
