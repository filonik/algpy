"""Fixed-metric CGA3 arithmetic, geometry and NumPy functions."""
import numpy as np
import pytest
import basic_checks as basic
import api_checks as api
import conversion_checks as conversions
import matrix_checks as matrices
import power_checks as powers
import binding_checks
import lattice_checks
import geometric_reference as geometric
from algebra_reference import UNARY
from field_reference import FUNCTIONS

ALGEBRA_NAME = 'cga3'


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
    api.check_specialization_structure('cga3', 32)

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


@pytest.mark.parametrize("parameters", [None], ids=str)
def test_power_order(common_algebra, parameters):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_power_order(common_algebra)

@pytest.mark.parametrize("parameters", [None], ids=str)
@pytest.mark.parametrize("operation", ["lpow", "rpow"])
def test_algebra_power_scipy(common_algebra, parameters, operation):
    if parameters is not None:
        common_algebra.K = parameters
    powers.check_algebra_power_scipy(common_algebra, operation)


def test_geometry_dimension(common_algebra):
    assert common_algebra.N == 32
    assert len(common_algebra.gens()) == 5


@pytest.fixture(scope="module")
def complex_basis():
    """Independent faithful 4x4 complex representation of Cl(4,1)."""
    identity = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    positive = [np.kron(x, identity), np.kron(z, x), np.kron(z, y), np.kron(z, z)]
    negative = 1j*positive[0]@positive[1]@positive[2]@positive[3]
    generators = [negative, *positive]
    basis = []
    for mask in range(32):
        value = np.eye(4, dtype=complex)
        for i, generator in enumerate(generators):
            if mask & (1 << i):
                value = value@generator
        basis.append(value)
    np.testing.assert_array_equal(basis[31], 1j*np.eye(4))
    return np.array(basis)


def matrix_of(algebra, basis, value):
    return np.einsum("k,kij->ij", algebra.to_float(value), basis)


def element_of(algebra, basis, matrix):
    coefficients = np.einsum("kij,ij->k", basis.conj(), matrix).real/4
    return algebra(*coefficients)


def test_complex_representation_and_polynomial(common_algebra, complex_basis):
    linear_sum_assignment = pytest.importorskip("scipy.optimize").linear_sum_assignment
    algebra = common_algebra
    rng = np.random.default_rng(147)
    for _ in range(12):
        a, b = [algebra(*row) for row in rng.uniform(-.2, .2, (2, 32))]
        matrix = matrix_of(algebra, complex_basis, a)
        np.testing.assert_allclose(matrix_of(algebra, complex_basis, a*b),
                                   matrix@matrix_of(algebra, complex_basis, b), atol=2e-15)
        polynomial = np.poly(matrix)[::-1]
        expected = np.convolve(polynomial, polynomial.conj()).real
        assert len(a.minpoly()) == 9 and len(a.eivals()) == 8
        np.testing.assert_allclose(a.minpoly(), expected, rtol=3e-12, atol=3e-14)
        assert type(a.norm()) is float
        assert a.norm() == pytest.approx(abs(np.linalg.det(matrix))**.25, rel=3e-12)
        spectrum = np.r_[np.linalg.eigvals(matrix), np.linalg.eigvals(matrix).conj()]
        rows, columns = linear_sum_assignment(abs(a.eivals()[:, None] - spectrum))
        np.testing.assert_allclose(a.eivals()[rows], spectrum[columns], rtol=3e-11, atol=3e-13)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("radius", [.3, 1e-4, 0])
def test_complex_jordan_block(common_algebra, complex_basis, name, function, radius):
    from field_reference import center, operation
    from scipy_reference import FUNCTIONS as references
    unitary, _ = np.linalg.qr(np.random.default_rng(148).normal(size=(4, 4)) +
                              1j*np.random.default_rng(149).normal(size=(4, 4)))
    jordan = (center(name) + .25j)*np.eye(4) + radius*np.diag(np.ones(3), 1)
    matrix = unitary@jordan@unitary.conj().T
    value = element_of(common_algebra, complex_basis, matrix)
    actual = matrix_of(common_algebra, complex_basis, function(value))
    expected = dict(references)[operation(name)](matrix)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_complex_centre(common_algebra, complex_basis, name, function):
    from field_reference import center
    w = center(name) + .25j
    value = common_algebra(w.real, x01234=w.imag)
    expected = function(w)*np.eye(4)
    np.testing.assert_allclose(matrix_of(common_algebra, complex_basis, function(value)), expected,
                               rtol=3e-13, atol=3e-14)


def test_3d_translator(common_algebra):
    em, ep, x, y, z = common_algebra.gens()
    origin, infinity = (em - ep)/2, em + ep
    displacement = 2*x - y + 3*z
    translator = common_algebra.one() - .5*displacement*infinity
    point = translator*origin*translator.inv()
    expected = origin + displacement + 7*infinity
    np.testing.assert_allclose(common_algebra.to_float(point), common_algebra.to_float(expected), atol=2e-14)
    np.testing.assert_allclose(common_algebra.to_float(point*point), 0, atol=2e-14)


def test_independent_arithmetic(common_algebra):
    geometric.check_arithmetic(common_algebra)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nilpotent_taylor(common_algebra, name, function):
    geometric.check_nilpotent(common_algebra, name, function)


def test_standalone_header(common_algebra, tmp_path):
    geometric.check_standalone_header(common_algebra, tmp_path)


def test_projector_sqrt(common_algebra):
    geometric.check_projector_sqrt(common_algebra)


def test_conformal_null_points(common_algebra):
    geometric.check_cga_points(common_algebra)


def test_lattice_matches_exterior5(common_algebra):
    from algpy.lattices import exterior5
    lattice_checks.check_matches_exterior(common_algebra, exterior5)

def test_lattice_constants(common_algebra):
    lattice_checks.check_constants(common_algebra)

def test_lattice_join_and_meet_reference(common_algebra):
    lattice_checks.check_join_and_meet_reference(common_algebra)

def test_lattice_identities_and_absorbing_elements(common_algebra):
    lattice_checks.check_identities_and_absorbing_elements(common_algebra)

def test_lattice_associativity(common_algebra):
    lattice_checks.check_associativity(common_algebra)

def test_lattice_generators(common_algebra):
    lattice_checks.check_generators(common_algebra)

def test_lattice_duals_on_basis(common_algebra):
    lattice_checks.check_duals_on_basis(common_algebra)

def test_lattice_duals_reference_and_inverse(common_algebra):
    lattice_checks.check_duals_reference_and_inverse(common_algebra)

def test_lattice_poincare_duality(common_algebra):
    lattice_checks.check_poincare_duality(common_algebra)

def test_lattice_dual_ufuncs(common_algebra):
    lattice_checks.check_dual_ufuncs(common_algebra)

def test_lattice_odd_dual(common_algebra):
    lattice_checks.check_odd_dual(common_algebra)

def test_lattice_vector_space_operations(common_algebra):
    lattice_checks.check_vector_space_operations(common_algebra)

def test_lattice_numpy_ufuncs(common_algebra):
    lattice_checks.check_numpy_ufuncs(common_algebra)

def test_lattice_cga_incidence_and_metric_meet(common_algebra):
    lattice_checks.check_cga_incidence_and_metric_meet(common_algebra)


def test_api_center_minpoly(common_algebra):
    from algpy.algebras import complex1
    assert common_algebra.center is complex1
    api.check_center_minpoly(common_algebra)
