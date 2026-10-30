"""Integer powers use multiplication, including at singular elements."""

from functools import lru_cache




import numpy as np

import pytest

from algpy import algebras

from algpy.algebras import general



@lru_cache(None)
def parameterized_type(name):
    return getattr(general, name)(*([-1]*int(name[-1])))

def sample(cls, center):
    row = np.random.default_rng(94).uniform(-.025, .025, cls.N)
    row[0] = center
    return cls(*row)

def check_repeated_multiplication(algebra, exponent):
    cls = algebra
    value = sample(cls, -.7)
    expected = cls.one()
    for _ in range(exponent):
        expected = expected*value
    result = value.ipow(exponent)
    assert type(result) is cls
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(expected), rtol=3e-12, atol=3e-14)

def check_negative_exponents(algebra, exponent):
    cls = algebra
    value = sample(cls, 1.5)
    inverse = value.inv()
    expected, positive = cls.one(), cls.one()
    for _ in range(-exponent):
        expected = expected*inverse
        positive = positive*value
    result = value.ipow(exponent)
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(expected), rtol=3e-12, atol=3e-14)
    np.testing.assert_allclose(cls.to_float(result*positive), cls.to_float(cls.one()), rtol=3e-12, atol=3e-14)

def check_exact_integer_range(algebra, exponent):
    cls = algebra
    expected = cls(-1 if exponent % 2 else 1)
    for argument in [exponent, np.int64(exponent)]:
        np.testing.assert_array_equal(cls.to_float(cls(-1).ipow(argument)), cls.to_float(expected))

def check_zero_nilpotents_and_idempotents(algebra):
    cls = algebra
    np.testing.assert_array_equal(cls.to_float(cls.zero().ipow(0)), cls.to_float(cls.one()))
    for exponent in [1, 2, 3, 13, 2**63-1]:
        np.testing.assert_array_equal(cls.to_float(cls.zero().ipow(exponent)), np.zeros(cls.N))
    for blade in cls.basis()[1:]:
        square = cls.to_float(blade*blade)
        if not square.any():
            np.testing.assert_array_equal(cls.to_float(blade.ipow(2)), np.zeros(cls.N))
            np.testing.assert_array_equal(cls.to_float(blade.ipow(13)), np.zeros(cls.N))
            for exponent in [-5, 0, 1, 3, 13]:
                np.testing.assert_allclose(cls.to_float((cls.one()+blade).ipow(exponent)),
                                           cls.to_float(cls.one()+exponent*blade), atol=3e-14)
            break
        if square[0] == 1 and not square[1:].any():
            projector = (cls.one()+blade)/2
            np.testing.assert_array_equal(cls.to_float(projector.ipow(13)), cls.to_float(projector))
            break

class Index:
    def __index__(self):
        return -3

def check_index_protocol_and_argument_errors(algebra):
    value = algebra(-2)
    assert value.ipow(Index()).to_re() == pytest.approx(-.125)
    assert value.ipow(np.uint64(3)).to_re() == -8
    assert value.ipow(True).to_re() == -2
    assert value.ipow(False).to_re() == 1
    for argument in [3., np.float64(3), 2.5, 3+0j, "3", None, [3]]:
        with pytest.raises(TypeError):
            value.ipow(argument)
    for argument in [-(2**63)-1, 2**63, np.uint64(2**64-1)]:
        with pytest.raises(OverflowError):
            value.ipow(argument)
    with pytest.raises(TypeError):
        value.ipow()
    with pytest.raises(TypeError):
        value.ipow(1, 2)

def check_current_parameters():
    cls = parameterized_type("clifford3")
    value = cls(0, 1)
    cls.K = (-1, -1, -1)
    before = cls.to_float(value.ipow(2))
    cls.K = (1, 0, 0)
    after = cls.to_float(value.ipow(2))
    assert before[0] == -1 and after[0] == 1


def check_algebra_power(algebra, operation):
    coefficients = np.linspace(.01, .04, algebra.N)
    coefficients[0] = 1.5
    base = algebra(*coefficients)
    exponent = algebra(*np.linspace(.2, .3, algebra.N))
    logarithm = np.log(base)
    product = exponent*logarithm if operation == "lpow" else logarithm*exponent
    expected = np.exp(product)
    actual = getattr(base, operation)(exponent)
    assert type(actual) is algebra
    assert np.isfinite(algebra.to_float(actual)).all()
    np.testing.assert_allclose(algebra.to_float(actual), algebra.to_float(expected),
                               rtol=3e-12, atol=3e-14)

    if operation == "rpow":
        for result in [base**exponent, np.power(base, exponent)]:
            np.testing.assert_allclose(algebra.to_float(result), algebra.to_float(actual),
                                       rtol=3e-12, atol=3e-14)
        values = np.array([base, base, base], dtype=algebra)[::2, None]
        exponents = np.array([exponent, algebra(.5)], dtype=algebra)
        out = np.empty((2, 4), dtype=algebra)[:, ::2]
        assert np.power(values, exponents, out=out) is out
        for row in out:
            for result, exponent in zip(row, exponents):
                np.testing.assert_allclose(algebra.to_float(result),
                                           algebra.to_float(base.rpow(exponent)),
                                           rtol=3e-12, atol=3e-14)

    for exponent in [.25, np.float32(.3), np.float64(.5), 2, True]:
        expected = np.exp(float(exponent)*np.log(base))
        np.testing.assert_allclose(algebra.to_float(getattr(base, operation)(exponent)),
                                   algebra.to_float(expected), rtol=3e-12, atol=3e-14)

    class Derived(algebra):
        def __new__(cls, *components):
            return super().__new__(cls, *components)

    derived = Derived(*coefficients)
    actual = getattr(derived, operation)(Derived(.3))
    assert type(actual) is Derived
    np.testing.assert_allclose(algebra.to_float(actual),
                               algebra.to_float(getattr(base, operation)(.3)), rtol=3e-12, atol=3e-14)
    for invalid in [None, object(), "2", [2], (2,), 1+0j, np.array([2])]:
        with pytest.raises(TypeError):
            getattr(base, operation)(invalid)
    with pytest.raises(OverflowError):
        getattr(base, operation)(1 << 20000)
    with pytest.raises(TypeError):
        getattr(base, operation)()
    with pytest.raises(TypeError):
        getattr(base, operation)(.2, .3)


def _generator_exp(square, coefficient):
    if square < 0:
        scale = np.sqrt(-square)
        return np.cos(coefficient*scale), np.sin(coefficient*scale)/scale
    if square > 0:
        scale = np.sqrt(square)
        return np.cosh(coefficient*scale), np.sinh(coefficient*scale)/scale
    return 1., coefficient


def check_power_order(algebra):
    # e1 and e2 anticommute in every noncommutative algebra here. With
    # a = exp(t*e1), b = e2, left/right powers are exp(-/+ t*e12).
    # Evaluate those two-dimensional exponentials directly from e_i**2.
    parameters = algebra.K if hasattr(algebra, "K") else (-1, -1)
    scalar, component = _generator_exp(parameters[0], .25)
    base = algebra(scalar, component)
    exponent = algebra(0, 0, 1)
    scalar, component = _generator_exp(-parameters[0]*parameters[1], .25)
    left, right = np.zeros((2, algebra.N))
    left[0] = right[0] = scalar
    left[3], right[3] = -component, component
    np.testing.assert_allclose(algebra.to_float(base.lpow(exponent)), left,
                               rtol=3e-12, atol=3e-14)
    for result in [base.rpow(exponent), base**exponent, np.power(base, exponent)]:
        np.testing.assert_allclose(algebra.to_float(result), right,
                                   rtol=3e-12, atol=3e-14)
    assert not np.allclose(left, right)


def check_algebra_power_scipy(algebra, operation):
    pytest.importorskip("scipy.linalg")
    from scipy.linalg import expm, logm
    coefficients = np.linspace(.01, .04, algebra.N)
    coefficients[0] = 1.5
    base = algebra(*coefficients)
    exponent = algebra(*np.linspace(.2, .3, algebra.N))
    logarithm = logm(base.lmat())
    product = exponent.lmat() @ logarithm if operation == "lpow" else logarithm @ exponent.lmat()
    expected = expm(product)
    actual = getattr(base, operation)(exponent).lmat()
    assert np.isfinite(actual).all() and np.isfinite(expected).all()
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)
