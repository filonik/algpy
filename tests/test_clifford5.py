"""Fixed-metric Clifford5 arithmetic, geometry and NumPy functions."""
import numpy as np
import pytest
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
import binding_checks
from geometric_reference import multiply
from algpy.algebras.general import clifford, clifford5
from functools import lru_cache
import itertools
from field_reference import center, operation
from algebra_reference import UNARY
from field_reference import FUNCTIONS

ALGEBRA_NAME = 'clifford5'
PARAMETERS = [(-1,)*5, (1,)*5, (0,)*5, (1, 1, 1, -1, 0), (-2, .5, 1.5, 0, -.75)]


@lru_cache(None)
def scalar_type():
    return clifford5(*([-1]*5))


@pytest.mark.parametrize("parameters", list(itertools.product((-1, 0, 1), repeat=5)), ids=str)
@pytest.mark.parametrize("name", ["exp", "log", "sqrt", "sin"])
def test_all_signatures(scipy_functions, parameters, name):
    algebra = scalar_type()
    algebra.K = tuple(p*s for p, s in zip(parameters, (.25, 4, 2.25, .5, 1.5)))
    row = np.random.default_rng(211).uniform(-.015, .015, 32)
    row[0] = center(name)
    value = algebra(*row)
    actual = getattr(np, name)(value)
    expected = scipy_functions[name](value.lmat())
    np.testing.assert_allclose(actual.lmat(), expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
def test_independent_arithmetic(parameters):
    algebra = scalar_type()
    algebra.K = parameters
    x, y = np.random.default_rng(212).uniform(-.2, .2, (2, 32))
    a, b = algebra(*x), algebra(*y)
    np.testing.assert_allclose(algebra.to_float(a*b), multiply(parameters, x, y), atol=2e-14)
    np.testing.assert_allclose(algebra.to_float(np.square(a)), multiply(parameters, x, x), atol=2e-14)
    assert len(a.minpoly()) == 9 and len(a.eivals()) == 8
    np.testing.assert_allclose(np.polynomial.polynomial.polyfromroots(a.eivals()), a.minpoly(),
                               rtol=3e-10, atol=3e-12)


@pytest.fixture(scope="module")
def dual_matrix_basis():
    """Cl(3,1,0) = Mat(4,R) over dual numbers, represented by 8x8 blocks."""
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    y = np.array([[0, 1], [-1, 0]])
    generators = [np.kron(x, np.eye(2)), np.kron(z, x), np.kron(z, z), np.kron(z, y)]
    pseudoscalar = generators[0]@generators[1]@generators[2]@generators[3]
    basis = []
    for mask in range(32):
        value = np.eye(4)
        for i, generator in enumerate(generators):
            if mask & (1 << i):
                value = value@generator
        if mask & 16:
            basis.append(np.block([[np.zeros((4, 4)), -value@pseudoscalar],
                                    [np.zeros((4, 4)), np.zeros((4, 4))]]))
        else:
            basis.append(np.kron(np.eye(2), value))
    return np.array(basis)


def dual_element(algebra, basis, matrix):
    coefficients = np.einsum("kij,ij->k", basis, matrix)/np.einsum("kij,kij->k", basis, basis)
    return algebra(*coefficients)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_eighth_order_nilpotent(dual_matrix_basis, name, function):
    algebra = scalar_type()
    algebra.K = (1, 1, 1, -1, 0)
    n = np.diag(np.ones(3), 1)
    link = np.zeros((4, 4)); link[3, 0] = 1
    nilpotent = np.block([[n, link], [np.zeros((4, 4)), n]])
    np.testing.assert_array_equal(np.linalg.matrix_power(nilpotent, 8), 0)
    assert np.linalg.matrix_power(nilpotent, 7).any()
    w = center(name)
    expected = function(w)*np.eye(8)
    power = np.eye(8)
    circle = np.exp(2j*np.pi*np.arange(512)/512)
    for order in range(1, 8):
        power = power@nilpotent
        coefficient = np.mean(function(w+.5*circle)*circle**(-order)).real/.5**order
        expected += coefficient*power
    value = dual_element(algebra, dual_matrix_basis, w*np.eye(8) + nilpotent)
    actual = np.einsum("k,kij->ij", algebra.to_float(function(value)), dual_matrix_basis)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_complex_fourth_order_jordan(scipy_functions, dual_matrix_basis, name, function):
    algebra = scalar_type()
    algebra.K = (1, 1, 1, -1, 0)
    rotation = np.array([[center(name), -.25], [.25, center(name)]])
    a = np.block([[rotation, .25*np.eye(2)], [np.zeros((2, 2)), rotation]])
    link = np.zeros((4, 4)); link[2:, :2] = .25*np.eye(2)
    matrix = np.block([[a, link], [np.zeros((4, 4)), a]])
    value = dual_element(algebra, dual_matrix_basis, matrix)
    actual = np.einsum("k,kij->ij", algebra.to_float(function(value)), dual_matrix_basis)
    expected = scipy_functions[operation(name)](matrix)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("scale", [1e-8, 1e-60, 1e-160])
@pytest.mark.parametrize("name", ["expm1", "log1p", "sqrt"])
def test_tiny_components(parameters, scale, name):
    algebra = scalar_type()
    algebra.K = parameters
    row = np.random.default_rng(214).uniform(-scale, scale, 32)
    row[0] = 1 if name == "sqrt" else scale
    delta = row.copy()
    if name == "sqrt": delta[0] = 0
    square = multiply(parameters, delta, delta)
    cube = multiply(parameters, square, delta)
    if name == "sqrt":
        expected = delta/2 - square/8 + cube/16
        expected[0] += 1
    elif name == "expm1": expected = row + square/2 + cube/6
    else: expected = row - square/2 + cube/3
    actual = algebra.to_float(getattr(np, name)(algebra(*row)))
    np.testing.assert_allclose(actual, expected, rtol=5e-12, atol=scale*1e-14)


def test_dispatch_mutability_and_pickle():
    import pickle
    algebra = clifford(-1, 1, 2, 0, .5)
    assert algebra.N == 32 and algebra.K == [-1, 1, 2, 0, .5]
    value = algebra(*np.linspace(-.1, .1, 32))
    restored = pickle.loads(pickle.dumps(value))
    assert type(restored).K == algebra.K
    np.testing.assert_array_equal(type(restored).to_float(restored), algebra.to_float(value))
    algebra.K = [1]*5
    np.testing.assert_array_equal(algebra.to_float(algebra.basis()[16]*algebra.basis()[16]),
                                  algebra.to_float(algebra.one()))
    with pytest.raises(ValueError): algebra.K = [1]*4
    with pytest.raises(TypeError): clifford(*([1]*6))


def test_negative_norm():
    algebra = scalar_type()
    algebra.K = (1, 1, 1, -1, -1)
    basis = algebra.basis()
    projector = (1 + basis[1])*(1 + basis[10])*(1 + basis[31])/8
    value = 1 - 2*projector
    assert value.minpoly()[0] == pytest.approx(-1)
    assert type(value.norm()) is complex
    assert value.norm() == pytest.approx(np.exp(1j*np.pi/8))


def test_sqrt_semisimple_zero_with_four_eigenvalues():
    algebra = scalar_type()
    algebra.K = (1, 1, 1, -1, 0)
    basis = algebra.basis()
    p, q = (1 + basis[1])/2, (1 + basis[10])/2
    # This exact unipotent similarity activates all five generators.
    s, inverse = 1 + .125*basis[22], 1 - .125*basis[22]
    value = s*(p + 2*q)*inverse
    expected = s*(p*(1-q) + np.sqrt(2)*(1-p)*q + np.sqrt(3)*p*q)*inverse
    with np.errstate(all="raise"):
        actual = np.sqrt(value)
    np.testing.assert_allclose(algebra.to_float(actual), algebra.to_float(expected), atol=3e-14)


def test_standalone_custom_callbacks(tmp_path):
    from pathlib import Path
    import shutil
    import subprocess
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    functions = ["static double f(double x) { return double_ipow(x, 9) + 2*x + 4; }"]
    factor = 1
    for order in range(1, 8):
        factor *= 10-order
        functions.append(f"static double d{order}(double x) {{ return {factor}*double_ipow(x, {9-order})" +
                         (" + 2" if order == 1 else "") + "; }")
    functions += ["static complex cf(complex x) { return complex_ladd(complex_add(complex_ipow(x, 9), complex_lmul(x, 2)), 4); }"]
    for order, factor in [(1, 9), (2, 72), (3, 504)]:
        functions.append(f"static complex cd{order}(complex x) {{ return complex_ladd(complex_lmul(complex_ipow(x, {9-order}), {factor}), " +
                         ("2" if order == 1 else "0") + "); }")
    source = tmp_path/"callbacks.c"
    source.write_text('#include "src/algpy/algebras/general/clifford/clifford5.h"\n#include <assert.h>\n' +
                     "\n".join(functions) + '''
int main(void) {
    vec5 parameters[] = {{{-1, 1, 1, 1, 1}}, {{1, 1, 1, 1, 1}}, {{1, 1, 1, -1, 0}}, {{0, 0, 0, 0, 0}}};
    clifford5 a = clifford5_from_re(.3);
    for (int i = 1; i < 32; i++) a.data[i] = .01*(i%5 - 2);
    for (int j = 0; j < 4; j++) {
        vec5 k = parameters[j];
        clifford5 actual = clifford5_apply(k, a, f, d1, d2, d3, d4, d5, d6, d7, cf, cd1, cd2, cd3);
        clifford5 expected = clifford5_ladd(vec32_add(clifford5_ipow(k, a, 9), vec32_scale(a, 2)), 4);
        for (int i = 0; i < 32; i++) assert(fabs(actual.data[i] - expected.data[i]) < 2e-10);
    }
    return 0;
}
''')
    executable = tmp_path/"callbacks"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)], check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)


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

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
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

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("exponent", [0, 1, 2, 3, 4, 5, 7, 13])
def test_ipow_repeated_multiplication(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_repeated_multiplication(common_algebra, exponent)

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("exponent", [-1, -2, -3, -5, -13])
def test_ipow_negative_exponents(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_negative_exponents(common_algebra, exponent)

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("exponent", [-(2**63), -(2**63)+1, -(2**53)-1, 2**53+1, 2**63-1])
def test_ipow_exact_integer_range(common_algebra, parameters, exponent):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_exact_integer_range(common_algebra, exponent)

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power(common_algebra, operation)


@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
def test_ipow_zero_nilpotents_and_idempotents(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_zero_nilpotents_and_idempotents(common_algebra)

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
def test_ipow_index_protocol_and_argument_errors(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_index_protocol_and_argument_errors(common_algebra)


@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
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


@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)

@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power_scipy(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power_scipy(common_algebra, operation)


@pytest.mark.parametrize("parameters", [[-1]*5, [1]*5, [.5, -2] + [1]*3, [0] + [1]*4, [2, 0, -1] + [-1]*2], ids=str)
def test_api_center_minpoly(common_algebra, parameters):
    api.check_center_minpoly(common_algebra, parameters)
