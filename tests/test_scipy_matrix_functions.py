"""Analytic functions commute with associative multiplication operators."""
from functools import lru_cache

import numpy as np
import pytest

from algpy import algebras
from algpy.algebras import general
from test_api_consistency import CONCRETE


linalg = pytest.importorskip("scipy.linalg")


# SciPy has no named inverse trigonometric matrix functions. Use its matrix
# logarithms and square roots: funm's eigenvalue interpolation loses derivatives
# on the nondiagonalizable matrices of dual and degenerate general algebras.
def arcsinm(a):
    identity = np.eye(len(a))
    return -1j*linalg.logm(1j*a + linalg.sqrtm(identity - a@a))


def arccosm(a):
    return np.pi/2*np.eye(len(a)) - arcsinm(a)


def arctanm(a):
    identity = np.eye(len(a))
    return (linalg.logm(identity + 1j*a) - linalg.logm(identity - 1j*a))/(2j)


def arcsinhm(a):
    return linalg.logm(a + linalg.sqrtm(a@a + np.eye(len(a))))


def arccoshm(a):
    identity = np.eye(len(a))
    return linalg.logm(a + linalg.sqrtm(a - identity)@linalg.sqrtm(a + identity))


def arctanhm(a):
    identity = np.eye(len(a))
    return (linalg.logm(identity + a) - linalg.logm(identity - a))/2


FUNCTIONS = [
    ("exp", linalg.expm),
    ("expm1", lambda a: linalg.expm(a) - np.eye(len(a))),
    ("log", linalg.logm),
    ("log1p", lambda a: linalg.logm(np.eye(len(a)) + a)),
    ("reciprocal", linalg.inv),
    ("square", lambda a: linalg.fractional_matrix_power(a, 2)),
    ("sqrt", linalg.sqrtm),
    ("sin", linalg.sinm), ("cos", linalg.cosm), ("tan", linalg.tanm),
    ("arcsin", arcsinm), ("arccos", arccosm), ("arctan", arctanm),
    ("sinh", linalg.sinhm), ("cosh", linalg.coshm), ("tanh", linalg.tanhm),
    ("arcsinh", arcsinhm), ("arccosh", arccoshm), ("arctanh", arctanhm),
]

INVERSE_FUNCTIONS = {"arcsin", "arccos", "arctan", "arcsinh", "arccosh", "arctanh"}
# Keep existing TODOs visible. Strict expected failures must be removed when
# these functions are implemented; a newly implemented function fails as XPASS.
UNIMPLEMENTED = {
    "quaternion": INVERSE_FUNCTIONS | {"sqrt", "sin", "cos", "tan", "sinh", "cosh", "tanh"},
    "cayleydickson2": INVERSE_FUNCTIONS,
    "clifford2": INVERSE_FUNCTIONS,
    "clifford3": {name for name, _ in FUNCTIONS} - {"reciprocal", "square"},
    "clifford4": {name for name, _ in FUNCTIONS} - {"reciprocal", "square"},
}


def parameter_sets(order):
    scales = (2, 0.5, 1.5, 0.75)[:order]
    yield tuple(-k for k in scales)
    yield scales
    yield (0,)*order
    if order > 1:
        yield (0, -2, 0.5, -1.5)[:order]


# Cayley-Dickson is associative only through order two in general. All Clifford
# and multicomplex orders here, and the remaining concrete types, are associative.
ALGEBRAS = [
    pytest.param((name, None), id=name)
    for name, _ in CONCRETE if name not in {"octonion", "sedenion"}
] + [
    pytest.param((f"{family}{order}", parameters), id=f"{family}{order}{parameters}")
    for family, orders in [("cayleydickson", (1, 2)), ("clifford", (1, 2, 3, 4)),
                           ("multicomplex", (1, 2, 3, 4))]
    for order in orders for parameters in parameter_sets(order)
]


@lru_cache(None)
def parameterized_type(name):
    # Factories allocate a new NumPy dtype on each call. Reuse one per order,
    # setting K for each case instead of allocating a type for every test.
    return getattr(general, name)(*([-1]*int(name[-1])))


@pytest.fixture(params=ALGEBRAS)
def algebra(request):
    name, parameters = request.param
    if parameters is None:
        return name, getattr(algebras, name)
    cls = parameterized_type(name)
    cls.K = parameters
    return name, cls


def sample_values(cls, operation):
    # The perturbation has operator norm at most 0.2 in both representations.
    # These centers keep spectra away from branch cuts, poles and singularities,
    # while exercising every coefficient (including higher nilpotent products).
    if operation == "arccosh":
        centers = [1.5]
    elif operation in {"log", "sqrt", "reciprocal"}:
        centers = [0.8]
    else:
        centers = [-0.35, 0, 0.35]
    coefficients = np.random.default_rng(81).uniform(-1, 1, cls.N)
    coefficients[0] = 0
    perturbation = cls(*coefficients)
    scale = max(1, np.linalg.norm(perturbation.lmat(), 1),
                np.linalg.norm(perturbation.rmat(), 1))
    coefficients *= 0.2/scale
    for center in centers:
        yield cls(center)
        coefficients[0] = center
        yield cls(*coefficients)


@pytest.mark.parametrize("operation,matrix_function", FUNCTIONS,
                         ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("operator", ["lmat", "rmat"])
def test_elementary_functions(algebra, operation, matrix_function, operator, request):
    name, cls = algebra
    if operation in UNIMPLEMENTED.get(name, set()):
        request.applymarker(pytest.mark.xfail(
            strict=True, raises=AssertionError, reason=f"TODO: {name}.{operation}"))
    for value in sample_values(cls, operation):
        result = getattr(np, operation)(value)
        actual = getattr(result, operator)()
        assert np.isfinite(actual).all()
        expected = matrix_function(getattr(value, operator)())
        assert np.isfinite(expected).all()
        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-12,
                                   equal_nan=False)
