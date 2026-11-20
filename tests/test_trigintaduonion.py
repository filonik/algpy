"""Tests for the 32-component Cayley-Dickson specialization."""
from algpy.algebras import sedenion, trigintaduonion
import numpy as np
import pytest
from field_reference import FUNCTIONS
from cayleydickson_reference import basis_coefficient, multiply, check_standalone_specialization_header
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_fixed_complex_slices(name, function):
    cls = trigintaduonion
    rng = np.random.default_rng(73)
    rows = rng.uniform(-0.3, 0.3, (3, cls.N))
    rows[:, 0] = [0.3, -0.4, 1.6]
    for real in [-2, -0.3, 0.3, 2]:
        row = np.zeros(cls.N)
        row[0] = real
        rows = np.vstack([rows, row])
    expected = []
    for row in rows:
        radius = np.linalg.norm(row[1:])
        z = function(complex(row[0], radius))
        result = np.zeros(cls.N)
        result[0] = z.real
        if radius:
            result[1:] = row[1:]/radius*z.imag
        else:
            result[1] = z.imag
        expected.append(result)
    values = cls.from_float(rows)
    np.testing.assert_allclose(cls.to_float(function(values)), expected, rtol=3e-13, atol=3e-14)
    for value, wanted in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(getattr(value, name)()), wanted,
                                   rtol=3e-13, atol=3e-14)


def test_specialized_arithmetic_is_independent():
    order = 5
    cls = trigintaduonion
    k = (-1,)*order
    a, b = [np.random.default_rng(seed).uniform(-0.3, 0.3, cls.N) for seed in [76, 77]]
    x, y = cls(*a), cls(*b)
    np.testing.assert_allclose(cls.to_float(x*y), multiply(k, a, b), atol=3e-14)
    np.testing.assert_allclose(cls.to_float(np.square(x)), multiply(k, a, a), atol=3e-14)
    np.testing.assert_allclose(cls.to_float(x*x.conjugate()), [a@a]+[0]*(cls.N-1), atol=3e-14)
    assert x.norm() == pytest.approx(np.linalg.norm(a))


@pytest.mark.parametrize("operation", ["sqrt", "expm1", "log1p"])
@pytest.mark.parametrize("axis", [1, 16, 31])
def test_small_components(operation, axis):
    cls = trigintaduonion
    for scale in [1e-12, 1e-60, 1e-160]:
        w = 1 if operation == "sqrt" else scale
        row = np.zeros(cls.N)
        row[0], row[axis] = w, scale
        value = cls(*row)
        z = complex(w, scale)
        # NumPy's complex log1p loses the real part for these tiny inputs.
        expected = z - z*z/2 + z*z*z/3 if operation == "log1p" else getattr(np, operation)(z)
        actual = cls.to_float(getattr(np, operation)(value))
        np.testing.assert_allclose(actual[[0, axis]], [expected.real, expected.imag], rtol=5e-14, atol=0)
        np.testing.assert_array_equal(np.delete(actual, [0, axis]), 0)


@pytest.mark.parametrize("operation", ["log", "sqrt"])
@pytest.mark.parametrize("axis", [1, 16, 31])
def test_large_directions(operation, axis):
    cls = trigintaduonion
    for scale in [1e150, 1e200, 1e300]:
        row = np.zeros(cls.N)
        row[0], row[axis] = 0.5, scale
        value = cls(*row)
        expected = getattr(np, operation)(complex(0.5, scale))
        actual = cls.to_float(getattr(np, operation)(value))
        np.testing.assert_allclose(actual[[0, axis]], [expected.real, expected.imag], rtol=5e-14, atol=0)
        np.testing.assert_array_equal(np.delete(actual, [0, axis]), 0)


def test_basis_products_against_doubling_reference():
    cls = trigintaduonion
    identity = np.eye(cls.N)
    expected = np.empty((cls.N, cls.N, cls.N))
    for i in range(cls.N):
        for j in range(cls.N):
            expected[i, j] = basis_coefficient((-1,)*5, i, j)*identity[i ^ j]
    np.testing.assert_array_equal(cls.to_float(cls.coeffs()), expected)


def test_embedded_sedenions():
    cls = trigintaduonion
    a, b = np.random.default_rng(79).uniform(-0.3, 0.3, (2, 16))
    for actual, expected in [(cls(*a)*cls(*b), sedenion(*a)*sedenion(*b)),
                             (np.exp(cls(*a)), np.exp(sedenion(*a))),
                             (np.log(cls(*a)), np.log(sedenion(*a)))]:
        coefficients = cls.to_float(actual)
        np.testing.assert_allclose(coefficients[:16], sedenion.to_float(expected), atol=3e-14)
        np.testing.assert_array_equal(coefficients[16:], 0)


def test_quadratic_polynomial_and_spectrum():
    cls = trigintaduonion
    row = np.random.default_rng(80).uniform(-0.3, 0.3, cls.N)
    a = cls(*row)
    polynomial = a.minpoly()
    assert polynomial.shape == (3,) and a.eivals().shape == (2,)
    np.testing.assert_allclose(polynomial, [row@row, -2*row[0], 1], atol=3e-14)
    np.testing.assert_allclose(cls.to_float(a*a + polynomial[1]*a + polynomial[0]), 0, atol=3e-14)
    np.testing.assert_allclose(a.eivals(), row[0] + np.array([-1j, 1j])*np.linalg.norm(row[1:]), atol=3e-14)
    np.testing.assert_allclose(cls.to_float(a*a.conjugate()), cls.to_float(a.conjugate()*a), atol=3e-14)


def test_zero_divisors_and_failure_of_alternativity():
    cls = trigintaduonion
    basis = cls.basis()
    a, b = basis[1] + basis[18], basis[4] - basis[23]
    assert a != cls.zero() and b != cls.zero()
    assert a*b == cls.zero()
    assert a.norm() == pytest.approx(np.sqrt(2))
    assert b.norm() == pytest.approx(np.sqrt(2))
    assert a.inv()*a == cls.one() and a*a.inv() == cls.one()
    assert (a*a)*b == -2*b
    assert a*(a*b) == cls.zero()
    np.testing.assert_array_equal(a.lmat() @ cls.to_float(b), 0)
    np.testing.assert_array_equal(cls.to_float(a) @ b.rmat(), 0)


def test_standalone_specialization_header(tmp_path):
    check_standalone_specialization_header("trigintaduonion", 32, tmp_path)

ALGEBRA_NAME = 'trigintaduonion'


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
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)


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
    api.check_specialization_structure('trigintaduonion', 32)

def test_matrix_nonassociative_products(common_algebra):
    matrices.check_nonassociative_products(common_algebra)




def test_binding_nan_hash(common_algebra):
    binding_checks.check_nan_hash_stable(common_algebra)


def test_binding_subclass_arithmetic(common_algebra):
    binding_checks.check_subclass_mixed_arithmetic(common_algebra)


def test_binding_byteorder(common_algebra):
    binding_checks.check_dtype_byteorder(common_algebra)


def test_binding_real_input_dtypes(common_algebra):
    binding_checks.check_real_input_dtypes(common_algebra)
