"""SciPy reference functions and safe input construction for algebra tests."""

import numpy as np

import pytest





linalg = pytest.importorskip("scipy.linalg")

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
