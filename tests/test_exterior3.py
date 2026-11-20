"""Exterior3 lattice: join, meet, Poincaré duals and NumPy bindings."""
import numpy as np
import pytest
import api_checks as api
import conversion_checks as conversions
import lattice_checks
from lattice_reference import lattice_type

LATTICE_NAME = 'exterior3'


@pytest.fixture
def lattice():
    return lattice_type(LATTICE_NAME)


def test_lattice_constants(lattice):
    lattice_checks.check_constants(lattice)

def test_lattice_join_and_meet_reference(lattice):
    lattice_checks.check_join_and_meet_reference(lattice)

def test_lattice_join_is_degenerate_clifford_product(lattice):
    lattice_checks.check_join_is_degenerate_clifford_product(lattice)

def test_lattice_identities_and_absorbing_elements(lattice):
    lattice_checks.check_identities_and_absorbing_elements(lattice)

def test_lattice_associativity(lattice):
    lattice_checks.check_associativity(lattice)

def test_lattice_generators(lattice):
    lattice_checks.check_generators(lattice)

def test_lattice_duals_on_basis(lattice):
    lattice_checks.check_duals_on_basis(lattice)

def test_lattice_duals_reference_and_inverse(lattice):
    lattice_checks.check_duals_reference_and_inverse(lattice)

def test_lattice_poincare_duality(lattice):
    lattice_checks.check_poincare_duality(lattice)

def test_lattice_dual_ufuncs(lattice):
    lattice_checks.check_dual_ufuncs(lattice)

def test_lattice_odd_dual(lattice):
    lattice_checks.check_odd_dual(lattice)

def test_lattice_vector_space_operations(lattice):
    lattice_checks.check_vector_space_operations(lattice)

def test_lattice_no_lattice_products(lattice):
    lattice_checks.check_no_lattice_products(lattice)

def test_lattice_numpy_ufuncs(lattice):
    lattice_checks.check_numpy_ufuncs(lattice)

def test_lattice_scalar_protocols(lattice):
    lattice_checks.check_scalar_protocols(lattice)

def test_conversion_native_methods_and_scalar_shape(lattice):
    conversions.check_native_methods_and_scalar_shape(lattice)

def test_conversion_shapes_and_contiguous_sharing(lattice):
    conversions.check_shapes_and_contiguous_sharing(lattice)

def test_conversion_strided_fortran_and_reversed_inputs(lattice):
    conversions.check_strided_fortran_and_reversed_inputs(lattice)

@pytest.mark.parametrize("dtype", [None, float, np.float32, np.int64, np.complex128, object,
                                  np.dtype(float).newbyteorder("S"), "U"])
def test_conversion_optional_output_dtype(lattice, dtype):
    conversions.check_optional_output_dtype(lattice, dtype)

def test_conversion_readonly_and_byte_swapped_inputs(lattice):
    conversions.check_readonly_and_byte_swapped_inputs(lattice)

def test_api_grade_projection(lattice):
    api.check_grade_projection(lattice)
