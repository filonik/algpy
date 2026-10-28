import numpy as np
import pytest
from algebra_suite import make_tests
from algpy.algebras import splitcomplex2

make_tests(globals(), splitcomplex2)


def _pos(*c):
    return splitcomplex2(*c)


def test_exp_log_sqrt_numpy():
    x = _pos(5.0, 0.5, 0.25, 0.125)
    arr = np.array([x, x * 2])
    for f in (np.exp, np.log, np.sqrt, np.sin, np.cos, np.tan,
              np.sinh, np.cosh, np.tanh):
        out = f(arr)
        assert out.dtype == splitcomplex2.dtype
    assert np.allclose(splitcomplex2.to_float(np.log(np.exp(x))), splitcomplex2.to_float(x))
    s = np.sqrt(x)
    assert np.allclose(splitcomplex2.to_float(s * s), splitcomplex2.to_float(x))
    assert np.allclose(
        splitcomplex2.to_float(np.sin(x) * np.sin(x) + np.cos(x) * np.cos(x)),
        splitcomplex2.to_float(splitcomplex2.one()),
    )
    assert np.allclose(splitcomplex2.to_float(np.power(x, 2)), splitcomplex2.to_float(x * x))


def test_pow_operator():
    x = _pos(5.0, 0.5, 0.25, 0.125)
    assert np.allclose(splitcomplex2.to_float(x ** 2), splitcomplex2.to_float(x * x))
    assert np.allclose(splitcomplex2.to_float(x ** 0.5), splitcomplex2.to_float(np.sqrt(x)))
