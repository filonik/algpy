"""Tests for the 4-component dual specialization."""
import math
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from algpy.algebras import dual2
from field_reference import reference
from dual_reference import FUNCTIONS, product, real_points, taylor, c_derivatives
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
from algebra_reference import UNARY
import numpy_checks
import binding_checks


@pytest.fixture(scope="module")
def algebra():
    return dual2


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_elementary_functions(algebra, name, function):
    cls = algebra
    coefficients = np.random.default_rng(51).uniform(-0.3, 0.3, (3, cls.N))
    coefficients[:, 0] = real_points(name)
    expected = np.array([taylor(function, row) for row in coefficients])
    values = cls.from_float(coefficients)
    storage = np.empty(2*len(values), dtype=cls)
    out = storage[::2]
    assert function(values[::-1], out=out) is out
    np.testing.assert_allclose(cls.to_float(out), expected[::-1], rtol=3e-10, atol=3e-10)
    for value, result in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(getattr(value, name)()), result,
                                   rtol=3e-10, atol=3e-10)
        np.testing.assert_allclose(cls.to_float(function(value)), result, rtol=3e-10, atol=3e-10)


@pytest.mark.parametrize("offset", range(1, 3))
@pytest.mark.parametrize("index,name,function", [(i, name, function)
                         for i, (name, function) in enumerate(FUNCTIONS)],
                         ids=[name for name, _ in FUNCTIONS])
def test_derivative_helpers(c_derivatives, offset, index, name, function):
    order = 2
    coefficients = np.random.default_rng(52).uniform(-0.3, 0.3, 2**order)
    coefficients[0] = real_points(name)[0]
    expected = taylor(function, coefficients, offset)
    actual = c_derivatives(order, index, offset, coefficients)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-10)


def test_pure_nilpotents(algebra):
    cls = algebra
    coefficients = np.zeros(cls.N)
    for generator in range(cls.N.bit_length()-1):
        coefficients[1 << generator] = 1
    value = cls.from_float(coefficients)
    # Every subset occurs degree! times in the corresponding power.
    np.testing.assert_array_equal(cls.to_float(np.exp(value)), np.ones(cls.N))
    expected = [0]+[(-1)**(i.bit_count()-1)*math.factorial(i.bit_count()-1)
                    for i in range(1, cls.N)]
    np.testing.assert_array_equal(cls.to_float(np.log1p(value)), expected)


@pytest.mark.parametrize("function", [np.expm1, np.log1p])
def test_small_inputs(algebra, function):
    cls = algebra
    coefficients = np.linspace(-1, 1, cls.N)*1e-20
    coefficients[0] = 1e-20
    result = function(cls.from_float(coefficients))
    np.testing.assert_allclose(cls.to_float(result), coefficients, rtol=2e-15, atol=1e-35)


@pytest.mark.parametrize("forward,inverse", [
    (np.sin, np.arcsin), (np.cos, np.arccos), (np.tan, np.arctan),
    (np.sinh, np.arcsinh), (np.cosh, np.arccosh), (np.tanh, np.arctanh),
    (np.exp, np.log), (np.square, np.sqrt),
])
def test_inverse_functions(algebra, forward, inverse):
    cls = algebra
    value = cls(1.5 if inverse is np.arccosh else 0.4,
                *np.linspace(-0.1, 0.1, cls.N-1))
    result = forward(inverse(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value), rtol=2e-12, atol=2e-14)


def test_power_and_mix(algebra):
    cls = algebra
    xp = np.r_[1.5, np.linspace(-0.1, 0.1, cls.N-1)]
    yp = np.r_[1.3, np.linspace(0.1, -0.1, cls.N-1)]
    x, y = cls.from_float(np.array([xp, yp]))
    expected = taylor(np.exp, product(taylor(np.log, xp), yp))
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=3e-10, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], expected,
                               rtol=3e-10, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(x**2), product(xp, xp), rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), taylor(np.exp, math.log(2)*xp),
                               rtol=3e-10, atol=3e-10)
    for t in [0, 0.4, 1]:
        expected = product(xp, taylor(np.exp, t*(taylor(np.log, yp)-taylor(np.log, xp))))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected, rtol=3e-10, atol=3e-10)


def test_arithmetic_and_metadata(algebra):
    cls = algebra
    xp = np.r_[1.5, np.random.default_rng(53).uniform(-0.2, 0.2, cls.N-1)]
    yp = np.r_[1.25, np.random.default_rng(54).uniform(-0.2, 0.2, cls.N-1)]
    x, y = cls.from_float(np.array([xp, yp]))
    for actual, expected in [(x*y, product(xp, yp)),
                             (x/y, product(xp, taylor(np.reciprocal, yp))),
                             (x.inv(), taylor(np.reciprocal, xp))]:
        np.testing.assert_allclose(cls.to_float(actual), expected, rtol=3e-10, atol=3e-10)
    for real in [0, -0.75, 1.5]:
        value = cls(real, *xp[1:])
        assert type(value.norm()) is float
        assert value.norm() == abs(real)
        np.testing.assert_array_equal(value.eivals(), np.full(cls.N, real, dtype=complex))
        expected = np.polynomial.polynomial.polyfromroots([real]*cls.N)
        np.testing.assert_allclose(value.minpoly(), expected, rtol=2e-13, atol=1e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))
    expected = np.zeros((cls.N, cls.N, cls.N))
    for i in range(cls.N):
        for j in range(cls.N):
            if not i & j:
                expected[i, j, i | j] = 1
    np.testing.assert_array_equal(cls.to_float(cls.coeffs()), expected)


def test_standalone_header_and_apply(tmp_path):
    order = 2
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix = f"dual{order}"
    lower = f"dual{order-1}" if order > 1 else "double"
    half = 2**(order-1)
    cube = 'return a*a*a;' if order == 1 else f'return {lower}_mul({lower}_square(a), a);'
    cube_deriv = 'return 3*a*a;' if order == 1 else f'return {lower}_lmul({lower}_square(a), 3);'
    norm_check = (f'assert(fabs(product.data[0] - {prefix}_norm2(a)) < 1e-12);'
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
#error Concrete dual specializations must not depend on general multicomplex headers
#endif
static inline {lower} cube({lower} a) {{ {cube} }}
static inline {lower} cube_deriv1({lower} a) {{ {cube_deriv} }}
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[1] = 0.25;
    a.data[{half}] += 0.125;
    {prefix} applied = {prefix}_apply(a, cube, cube_deriv1);
    {prefix} expected = {prefix}_mul({prefix}_square(a), a);
    for (int i = 0; i < {2**order}; i++) {{
        assert(fabs(applied.data[i] - expected.data[i]) < 1e-12);
    }}
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    {norm_check}
    assert(fabs({prefix}_norm{2**order}(a) - pow({prefix}_norm(a), {2**order})) < 1e-9);
    {prefix} two = {prefix}_from_re(2), three = {prefix}_from_re(3);
    assert(fabs({prefix}_rpows(2, three).data[0] - 8) < 1e-12);
    assert(fabs({prefix}_lpows(two, 3).data[0] - 8) < 1e-12);
    assert(fabs({prefix}_mix(two, {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)


@pytest.mark.parametrize("name,function", FUNCTIONS[:4]+FUNCTIONS[6:],
                         ids=[name for name, _ in FUNCTIONS[:4]+FUNCTIONS[6:]])
def test_field_derivative_callers(name, function):
    x = 1.5 if name in ["log", "sqrt", "acosh"] else 0.4
    first = reference(function, complex(x), 1).real
    second = reference(function, complex(x), 2).real
    coefficients = [x, 0.3, -0.2, 0.1]
    expected = [function(x), 0.3*first, -0.2*first, 0.1*first-0.06*second]
    actual = getattr(dual2(*coefficients), name)()
    np.testing.assert_allclose(dual2.to_float(actual), expected, rtol=3e-10, atol=3e-10)

ALGEBRA_NAME = 'dual2'


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
    api.check_specialization_structure('dual2', 4)


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
