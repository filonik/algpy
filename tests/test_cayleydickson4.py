from algpy.algebras.general import cayleydickson4
import itertools
import numpy as np
import pytest
from field_reference import FUNCTIONS, reference
from cayleydickson_reference import (scalar_type, basis_coefficient, multiply,
                                    sample, matrix_functions, check_general_apply_custom_callbacks)
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import binding_checks


@pytest.mark.parametrize("parameters", list(itertools.product([-1, 0, 1], repeat=4)), ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_general_functions(matrix_functions, parameters, name, function):
    cls = scalar_type(len(parameters))
    operation = {"inv": "reciprocal", "asin": "arcsin", "acos": "arccos", "atan": "arctan",
                 "asinh": "arcsinh", "acosh": "arccosh", "atanh": "arctanh"}.get(name, name)
    # Unequal, nonunit parameters exercise the general path as well as dispatch.
    scaled = tuple(k*m for k, m in zip(parameters, [0.25, 4, 2.25, 9]))
    for k in dict.fromkeys([parameters, scaled]):
        cls.K = k
        rows = np.array([sample(k, name, seed) for seed in [71, 72]])
        expected = []
        for row in rows:
            q = sum(basis_coefficient(k, i, i)*row[i]**2 for i in range(1, cls.N))
            # In the associative span of 1,v, multiplication by w+v has this
            # 2x2 matrix. Its first column is (constant, coefficient of v).
            matrix = np.array([[row[0], q], [1, row[0]]])
            value = matrix_functions[operation](matrix)
            wanted = row*value[1, 0].real
            wanted[0] = value[0, 0].real
            np.testing.assert_allclose(value.imag, 0, atol=2e-12)
            expected.append(wanted)
        values = cls.from_float(rows)
        output = np.empty(4, dtype=cls)[::2]
        assert function(values[::-1], out=output) is output
        np.testing.assert_allclose(cls.to_float(output), expected[::-1], rtol=3e-10, atol=3e-12)
        for value, wanted in zip(values, expected):
            np.testing.assert_allclose(cls.to_float(getattr(value, name)()), wanted,
                                       rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("parameters", ["negative", "positive", "zero", "mixed"])
def test_arithmetic_and_minimal_polynomial(parameters):
    order = 4
    k = {"negative": (-1,)*order, "positive": (1,)*order,
         "zero": (0,)*order, "mixed": (-2, 0, 0.5, -1.5)[:order]}[parameters]
    cls = scalar_type(order)
    cls.K = k
    a, b = [np.random.default_rng(seed).uniform(-0.2, 0.2, cls.N) for seed in [74, 75]]
    a[0], b[0] = 1.5, 1.2
    x, y = cls(*a), cls(*b)
    product = multiply(k, a, b)
    square = multiply(k, a, a)
    np.testing.assert_allclose(cls.to_float(x*y), product, atol=3e-14)
    np.testing.assert_allclose(cls.to_float(np.square(x)), square, atol=3e-14)
    conjugate = a.copy()
    conjugate[1:] *= -1
    norm2 = multiply(k, a, conjugate)[0]
    np.testing.assert_allclose(cls.to_float(x.inv()), conjugate/norm2, atol=3e-14)
    polynomial = x.minpoly()
    np.testing.assert_allclose(polynomial, [norm2, -2*a[0], 1], atol=3e-14)
    np.testing.assert_allclose(cls.to_float(np.square(x) - 2*a[0]*x + cls(norm2)), 0, atol=3e-14)
    np.testing.assert_allclose(np.sort_complex(x.eivals()),
                               np.sort_complex(np.roots(polynomial[::-1])), atol=2e-12)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nonzero_nilpotent_slices(name, function):
    order = 4
    cls = scalar_type(order)
    cls.K = [1]*order
    w = 2 if name == "acosh" else 1.5 if name in {"log", "sqrt", "inv"} else 0.3
    row = np.zeros(cls.N)
    row[0], row[1], row[3] = w, 0.25, 0.25  # (e0 + e0*e1)^2 = 0
    expected = row*reference(function, complex(w), 1).real
    expected[0] = function(w)
    np.testing.assert_allclose(cls.to_float(function(cls(*row))), expected, rtol=3e-12, atol=3e-13)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nearly_nilpotent_slices(name, function):
    order = 4
    cls = scalar_type(order)
    cls.K = [1e-24]+[0]*(order-1)
    w = 2 if name == "acosh" else 1.5 if name in {"log", "sqrt", "inv"} else 0.3
    row = np.zeros(cls.N)
    row[0], row[1], row[2] = w, 1, 1
    # v^2 = 1e-24, so deviations from the derivative limit are below double
    # precision here, even though the imaginary coefficients are of unit size.
    expected = row*reference(function, complex(w), 1).real
    expected[0] = function(w)
    np.testing.assert_allclose(cls.to_float(function(cls(*row))), expected, rtol=3e-12, atol=3e-13)


def test_real_continuation_and_invalid_real_slices():
    order = 4
    cls = scalar_type(order)
    cls.K = [4, 9]+[0]*(order-2)
    # The first negative-square basis vector is e0*e1, whose square is -36.
    root = np.zeros(cls.N)
    root[3] = 1/3
    np.testing.assert_allclose(cls.to_float(np.sqrt(cls(-4))), root, atol=3e-15)
    logarithm = np.zeros(cls.N)
    logarithm[0], logarithm[3] = np.log(4), np.pi/6
    np.testing.assert_allclose(cls.to_float(np.log(cls(-4))), logarithm, atol=3e-15)
    cls.K = [0]*order
    assert np.isnan(cls.to_float(np.sqrt(cls(-4)))).all()
    assert np.isnan(cls.to_float(np.log(cls(-4)))).all()
    nil = cls.one() + 0.25*cls.basis()[1]
    np.testing.assert_allclose(cls.to_float(nil**3), cls.to_float(cls.one()+0.75*cls.basis()[1]), atol=3e-14)
    nil = -cls.one() + 0.25*cls.basis()[1]
    np.testing.assert_allclose(cls.to_float(nil**3), cls.to_float(-cls.one()+0.75*cls.basis()[1]), atol=3e-14)
    for operation in [np.log, np.sqrt]:
        with np.errstate(divide="ignore", invalid="ignore"):
            assert np.isnan(cls.to_float(operation(cls.basis()[1]))).all()


def test_powers_and_mix():
    order = 4
    cls = scalar_type(order)
    for k in [(-1,)*order, (0,)*order, (-2, 0, 0.5, -1.5)[:order]]:
        cls.K = k
        x = cls(*sample(k, "log", 78))
        y = cls(*sample(k, "log", 79))
        np.testing.assert_allclose(cls.to_float(x**y), cls.to_float(np.exp(np.log(x)*y)), atol=3e-13)
        for b in [-2, -0.5, 0, 0.5, 1, 2, 3]:
            np.testing.assert_allclose(cls.to_float(x**b), cls.to_float(np.exp(np.log(x)*b)), atol=3e-12)
        for t in [0, 0.4, 1]:
            np.testing.assert_allclose(cls.to_float(x.mix(y, t)),
                                       cls.to_float(x*((x.inv()*y)**t)), atol=3e-12)
        np.testing.assert_array_equal(cls.to_float(cls.zero()**0), cls.to_float(cls.one()))
        np.testing.assert_array_equal(cls.to_float(cls.zero()**2), cls.to_float(cls.zero()))


def test_general_apply_custom_callbacks(tmp_path):
    check_general_apply_custom_callbacks(4, tmp_path)

def test_scipy_requires_identity_action(common_algebra, scipy_functions):
    cls = common_algebra
    cls.K = (-1, -1, -1, -1)
    row = np.random.default_rng(106).uniform(-.025, .025, 16)
    row[0] = .35
    value = cls(*row)
    result = np.exp(value)
    for operator in ["lmat", "rmat"]:
        reference = scipy_functions["exp"](getattr(value, operator)())
        # This example really distinguishes the valid identity-action check
        # from an invalid full-matrix comparison in a nonassociative algebra.
        assert np.max(np.abs(getattr(result, operator)() - reference)) > 1e-4
        expected = reference[:, 0] if operator == "lmat" else reference[0, :]
        np.testing.assert_allclose(cls.to_float(result), expected, rtol=2e-10, atol=2e-12)

ALGEBRA_NAME = 'cayleydickson4'


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

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)


@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", [(-1, -1, -1, -1), (1, 1, 1, 1), (0, 0, 0, 0), (0, -2, 0.5, -1.5)], ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)

def test_norm_parameter_type():
    api.check_general_norm_type_does_not_follow_parameters('cayleydickson', 4)

@pytest.mark.parametrize("parameters", [(-2, -2, -2, -2), (-2, -2, -2, 0), (-2, -2, -2, 0.5), (-2, -2, 0, -2), (-2, -2, 0, 0), (-2, -2, 0, 0.5), (-2, -2, 0.5, -2), (-2, -2, 0.5, 0), (-2, -2, 0.5, 0.5), (-2, 0, -2, -2), (-2, 0, -2, 0), (-2, 0, -2, 0.5), (-2, 0, 0, -2), (-2, 0, 0, 0), (-2, 0, 0, 0.5), (-2, 0, 0.5, -2), (-2, 0, 0.5, 0), (-2, 0, 0.5, 0.5), (-2, 0.5, -2, -2), (-2, 0.5, -2, 0), (-2, 0.5, -2, 0.5), (-2, 0.5, 0, -2), (-2, 0.5, 0, 0), (-2, 0.5, 0, 0.5), (-2, 0.5, 0.5, -2), (-2, 0.5, 0.5, 0), (-2, 0.5, 0.5, 0.5), (0, -2, -2, -2), (0, -2, -2, 0), (0, -2, -2, 0.5), (0, -2, 0, -2), (0, -2, 0, 0), (0, -2, 0, 0.5), (0, -2, 0.5, -2), (0, -2, 0.5, 0), (0, -2, 0.5, 0.5), (0, 0, -2, -2), (0, 0, -2, 0), (0, 0, -2, 0.5), (0, 0, 0, -2), (0, 0, 0, 0), (0, 0, 0, 0.5), (0, 0, 0.5, -2), (0, 0, 0.5, 0), (0, 0, 0.5, 0.5), (0, 0.5, -2, -2), (0, 0.5, -2, 0), (0, 0.5, -2, 0.5), (0, 0.5, 0, -2), (0, 0.5, 0, 0), (0, 0.5, 0, 0.5), (0, 0.5, 0.5, -2), (0, 0.5, 0.5, 0), (0, 0.5, 0.5, 0.5), (0.5, -2, -2, -2), (0.5, -2, -2, 0), (0.5, -2, -2, 0.5), (0.5, -2, 0, -2), (0.5, -2, 0, 0), (0.5, -2, 0, 0.5), (0.5, -2, 0.5, -2), (0.5, -2, 0.5, 0), (0.5, -2, 0.5, 0.5), (0.5, 0, -2, -2), (0.5, 0, -2, 0), (0.5, 0, -2, 0.5), (0.5, 0, 0, -2), (0.5, 0, 0, 0), (0.5, 0, 0, 0.5), (0.5, 0, 0.5, -2), (0.5, 0, 0.5, 0), (0.5, 0, 0.5, 0.5), (0.5, 0.5, -2, -2), (0.5, 0.5, -2, 0), (0.5, 0.5, -2, 0.5), (0.5, 0.5, 0, -2), (0.5, 0.5, 0, 0), (0.5, 0.5, 0, 0.5), (0.5, 0.5, 0.5, -2), (0.5, 0.5, 0.5, 0), (0.5, 0.5, 0.5, 0.5)], ids=str)
def test_matrix_parameters(parameters):
    matrices.check_general_parameters('cayleydickson', parameters)

def test_matrix_current_parameters():
    matrices.check_matrices_read_current_parameters('cayleydickson', 4)

@pytest.mark.parametrize("parameters", [(-2, -0.5, -1.5, -0.75), (2, 0.5, 1.5, 0.75), (0, 0, 0, 0), (0, -2, 0.5, -1.5), (-1, -2, -3, -4), (-1, -1, -1, -1), (-2, 0, 0.5, -1.5), (1, 1, 1, 1)], ids=str)
@pytest.mark.parametrize("operation", ['exp', 'expm1', 'log', 'log1p', 'reciprocal', 'square', 'sqrt', 'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh'])
@pytest.mark.parametrize("operator", ["lmat", "rmat"])
def test_scipy_elementary_functions(common_algebra, scipy_functions, parameters, operation, operator):
    from scipy_reference import sample_values
    if parameters is not None:
        common_algebra.K = parameters
    for value in sample_values(common_algebra, operation):
        result = getattr(np, operation)(value)
        expected = scipy_functions[operation](getattr(value, operator)())
        actual = common_algebra.to_float(result)
        expected = expected[:, 0] if operator == "lmat" else expected[0, :]
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


@pytest.mark.parametrize("reference", ['matrix', 'slice'])
def test_benchmark_entry_point(tmp_path, reference):
    from benchmark_checks import check_benchmark_entry_point
    check_benchmark_entry_point(tmp_path, ALGEBRA_NAME, reference)
