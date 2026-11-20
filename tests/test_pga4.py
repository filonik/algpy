"""Fixed-metric PGA4 arithmetic, geometry and NumPy functions."""
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

ALGEBRA_NAME = 'pga4'


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
    api.check_specialization_structure('pga4', 32)

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


def test_independent_arithmetic(common_algebra):
    geometric.check_arithmetic(common_algebra)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nilpotent_taylor(common_algebra, name, function):
    geometric.check_nilpotent(common_algebra, name, function)


def test_standalone_header(common_algebra, tmp_path):
    geometric.check_standalone_header(common_algebra, tmp_path, degree=9)


def test_projector_sqrt(common_algebra):
    geometric.check_projector_sqrt(common_algebra)


def test_null_generator_and_radical(common_algebra):
    geometric.check_pga_radical(common_algebra)


def test_projected_nilpotent_sqrt(common_algebra):
    geometric.check_pga_projected_nilpotent_sqrt(common_algebra)


@pytest.fixture(scope="module")
def complex_basis():
    return geometric.pga_matrix_basis(4)


def as_matrix(algebra, basis, value):
    return np.einsum("k,kij->ij", algebra.to_float(value), basis)


def test_faithful_representation_and_polynomial(common_algebra, complex_basis):
    assignment = pytest.importorskip("scipy.optimize").linear_sum_assignment
    cls = common_algebra
    rng = np.random.default_rng(241)
    size = complex_basis.shape[-1]//2
    for _ in range(10):
        a, b = [cls(*row) for row in rng.uniform(-.2, .2, (2, cls.N))]
        matrix = as_matrix(cls, complex_basis, a)
        np.testing.assert_allclose(as_matrix(cls, complex_basis, a*b),
                                   matrix@as_matrix(cls, complex_basis, b), atol=4e-15)
        quotient = matrix[:size, :size]
        polynomial = np.poly(quotient)[::-1].real
        if 4 == 4:
            polynomial = np.convolve(polynomial, polynomial)
        assert len(a.minpoly()) == 9 and len(a.eivals()) == 8
        np.testing.assert_allclose(a.minpoly(), polynomial, rtol=3e-11, atol=3e-13)
        roots = np.linalg.eigvals(quotient)
        if 4 == 4:
            roots = np.tile(roots, 2)
        rows, cols = assignment(abs(a.eivals()[:, None] - roots))
        np.testing.assert_allclose(a.eivals()[rows], roots[cols], rtol=3e-10, atol=3e-12)
        assert type(a.norm()) is float
        assert a.norm() == pytest.approx(abs(np.linalg.det(quotient))**(1/size), rel=3e-11)
        value = cls.one()
        for coefficient in a.minpoly()[-2::-1]:
            value = value*a + coefficient
        np.testing.assert_allclose(cls.to_float(value), 0, atol=3e-12)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_faithful_matrix_functions(common_algebra, complex_basis, scipy_functions, name, function):
    from field_reference import center, operation
    cls = common_algebra
    row = np.random.default_rng(242).uniform(-.05, .05, cls.N)
    row[0] = center(name)
    value = cls(*row)
    result = function(value)
    expected = scipy_functions[operation(name)](as_matrix(cls, complex_basis, value))
    np.testing.assert_allclose(as_matrix(cls, complex_basis, result), expected, rtol=3e-10, atol=3e-12)


def test_generic_multiplicities_are_retained(common_algebra):
    cls = common_algebra
    value = cls(2)
    expected = np.polynomial.polynomial.polyfromroots([2]*8)
    np.testing.assert_array_equal(value.minpoly(), expected)
    np.testing.assert_array_equal(value.eivals(), [2]*8)


def test_full_radical_nilpotent_and_inverse(common_algebra):
    cls = common_algebra
    coefficients = np.zeros(cls.N)
    coefficients[1::2] = np.linspace(-.3, .3, cls.N//2)
    radical = cls(*coefficients)
    np.testing.assert_array_equal(cls.to_float(radical.ipow(2)), 0)
    np.testing.assert_allclose(cls.to_float(np.exp(radical)), cls.to_float(1+radical), atol=2e-14)
    np.testing.assert_allclose(cls.to_float((1+radical).inv()), cls.to_float(1-radical), atol=2e-14)


def test_constructor_keywords_and_pickle(common_algebra):
    import pickle
    cls = common_algebra
    last = "x"+"".join(str(i) for i in range(5))
    value = cls(x=1, x0=2, **{last: 3})
    assert value[0] == 1 and value[1] == 2 and value[-1] == 3
    assert value == eval(repr(value), {cls.__name__: cls})
    assert value == pickle.loads(pickle.dumps(value))


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

def test_lattice_pga_incidence(common_algebra):
    lattice_checks.check_pga_incidence(common_algebra)


def test_api_center_minpoly(common_algebra):
    from algpy.algebras import dual1
    assert common_algebra.center is dual1
    api.check_center_minpoly(common_algebra)
