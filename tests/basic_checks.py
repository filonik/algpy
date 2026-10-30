"""Algebra identities shared by the explicit per-algebra test files."""
import numpy as np
import pytest


def _dim(algebra):
    return algebra.N

def _comp(algebra, x):
    return np.asarray(algebra.to_float(x), dtype=float).reshape(-1)

def _sample(algebra, seed):
    rng = np.random.default_rng(seed)
    return algebra(*rng.uniform(0.5, 2.0, size=_dim(algebra)))

def close(algebra, a, b, tol=1e-8):
    return np.allclose(_comp(algebra, a), _comp(algebra, b), rtol=tol, atol=tol)

def check_one(alg):
    x = _sample(alg, 0)
    one = alg.one()
    assert _comp(alg, one)[0] == 1 and not _comp(alg, one)[1:].any()
    assert close(alg, one * x, x) and close(alg, x * one, x)

def check_zero(alg):
    x = _sample(alg, 1)
    zero = alg.zero()
    assert not _comp(alg, zero).any()
    assert close(alg, x + zero, x) and close(alg, x * zero, zero) and close(alg, zero * x, zero)

def check_add(alg):
    x, y = _sample(alg, 2), _sample(alg, 3)
    assert close(alg, x + y, y + x)
    assert np.allclose(_comp(alg, x + y), _comp(alg, x) + _comp(alg, y))

def check_sub(alg):
    x, y = _sample(alg, 4), _sample(alg, 5)
    assert np.allclose(_comp(alg, x - y), _comp(alg, x) - _comp(alg, y))
    assert close(alg, x - x, alg.zero())

def check_mul(alg):
    x, y, z = _sample(alg, 6), _sample(alg, 7), _sample(alg, 8)
    assert close(alg, x * (y + z), x * y + x * z)
    assert close(alg, (x + y) * z, x * z + y * z)
    assert close(alg, x * 2, x + x)
    a = np.array([x, y], dtype=alg)
    assert close(alg, (a * a)[0], x * x) and close(alg, (a * a)[1], y * y)

def check_inv(alg):
    x = _sample(alg, 9)
    one = alg.one()
    inv = one / x
    assert close(alg, x * inv, one) and close(alg, inv * x, one)
    if hasattr(x, "inv"):
        assert close(alg, x.inv(), inv)

def check_div(alg):
    x, y = _sample(alg, 10), _sample(alg, 11)
    assert close(alg, x / x, alg.one())
    assert close(alg, x / y, x * (alg.one() / y))
    assert close(alg, x / 2, x * 0.5)

def check_norm(alg):
    x = _sample(alg, 12)
    n = x.norm()
    norm_type = complex if hasattr(alg, "K") or alg.__name__.startswith("splitcomplex") else float
    assert type(n) is norm_type
    assert type(alg.one().norm()) is norm_type
    assert alg.one().norm() == pytest.approx(1)
    assert alg.zero().norm() == 0
    assert (alg.one() * 2).norm() == pytest.approx(2)
    if isinstance(n, float):
        assert n >= 0

def check_minpoly_eivals(alg):
    z = _sample(alg, 13)
    assert np.allclose(
        np.polynomial.Polynomial(z.minpoly()).coef,
        np.polynomial.Polynomial.fromroots(z.eivals()).coef,
    )

def check_basis(alg):
    n = _dim(alg)
    basis = alg.basis()
    assert basis.shape == (n,)
    assert np.array_equal(alg.to_float(basis), np.eye(n))
    assert close(alg, basis[0], alg.one())

def check_coeffs(alg):
    n = _dim(alg)
    basis, coeffs = alg.basis(), alg.coeffs()
    assert coeffs.shape == (n, n)
    for i in range(n):
        for j in range(n):
            assert close(alg, coeffs[i, j], basis[i] * basis[j])

def check_float_conversion(alg):
    n = _dim(alg)
    x = _sample(alg, 14)
    assert alg.to_float(x).tolist() == [x[i] for i in range(n)]
    assert alg.from_float(alg.to_float(x)) == x
    components = np.arange(2 * n, dtype=float).reshape(2, n)
    values = alg.from_float(components)
    assert values.dtype == alg.dtype
    assert np.array_equal(alg.to_float(values), components)
    coeffs = alg.coeffs()
    assert (alg.from_float(alg.to_float(coeffs)) == coeffs).all()

def check_k_mutable(alg):
    if not hasattr(alg, "K"):
        pytest.skip("algebra has no parameters")
    original = list(alg.K)
    try:
        e1 = alg.basis()[1]
        alg.K = [k + 1.0 for k in original]
        assert alg.K == [k + 1.0 for k in original]
        assert close(alg, e1 * e1, alg.one() * alg.K[0])
    finally:
        alg.K = original
    assert alg.K == original
