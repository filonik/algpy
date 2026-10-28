from algebra_suite import make_tests
from algpy.algebras import complex1
import numpy as np
import pytest


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


make_tests(globals(), complex1)
