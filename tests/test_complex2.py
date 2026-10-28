from algebra_suite import make_tests
from algpy.algebras import complex2
import numpy as np
import pytest


FUNCTIONS = [
    ("cos", "cos"), ("sin", "sin"), ("tan", "tan"),
    ("arccos", "acos"), ("arcsin", "asin"), ("arctan", "atan"),
    ("cosh", "cosh"), ("sinh", "sinh"), ("tanh", "tanh"),
    ("arccosh", "acosh"), ("arcsinh", "asinh"), ("arctanh", "atanh"),
]


def component_reference(function, coefficients):
    w, x, y, z = np.asarray(coefficients).T
    p = function((w-z) + 1j*(x+y))
    q = function((w+z) + 1j*(x-y))
    return np.stack(((p.real+q.real)/2, (p.imag+q.imag)/2,
                     (p.imag-q.imag)/2, (q.real-p.real)/2), axis=-1)


@pytest.mark.parametrize("ufunc,method", FUNCTIONS)
def test_elementary_functions(ufunc, method):
    rng = np.random.default_rng(22)
    coefficients = np.vstack((
        rng.uniform(-2, 2, (16, 4)),
        # Real branch cuts and a nonzero zero divisor (p=0, q=2).
        [2, 0, 0, 0], [-2, 0, 0, 0], [0, 2, 0, 0], [0, -2, 0, 0],
        [1, 0, 0, 1], [0, 0, 0, 0],
    ))
    function = getattr(np, ufunc)
    expected = component_reference(function, coefficients)
    values = complex2.from_float(coefficients)
    # Exercise strided input and output in the registered NumPy loop.
    storage = np.empty(2*len(values), dtype=complex2)
    out = storage[::2]
    result = function(values[::-1], out=out)
    assert result is out
    np.testing.assert_allclose(complex2.to_float(out), expected[::-1],
                               rtol=2e-13, atol=2e-14)
    for value, reference in zip(values, expected):
        np.testing.assert_allclose(complex2.to_float(getattr(value, method)()),
                                   reference, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(complex2.to_float(function(value)),
                                   reference, rtol=2e-13, atol=2e-14)


@pytest.mark.parametrize("forward,inverse", [
    ("sin", "arcsin"), ("cos", "arccos"), ("tan", "arctan"),
    ("sinh", "arcsinh"), ("cosh", "arccosh"), ("tanh", "arctanh"),
])
def test_inverse_functions(forward, inverse):
    value = complex2(0.3, 0.2, -0.1, 0.05)
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(complex2.to_float(result), complex2.to_float(value),
                               rtol=2e-13, atol=2e-14)


def test_trigonometric_identities():
    values = complex2.from_float([[0.3, 0.2, -0.1, 0.05], [1, 0, 0, 1]])
    expected = np.tile([1, 0, 0, 0], (2, 1))
    np.testing.assert_allclose(complex2.to_float(np.square(np.sin(values)) + np.square(np.cos(values))),
                               expected, atol=2e-14)
    np.testing.assert_allclose(complex2.to_float(np.square(np.cosh(values)) - np.square(np.sinh(values))),
                               expected, atol=2e-14)

make_tests(globals(), complex2)
