"""Algebra identities shared by the explicit per-algebra test files."""
import inspect
import numpy as np
import pytest
from algebra_reference import FIXED_METRICS, expected_norm_type


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
    norm_type = expected_norm_type(alg)
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

def check_gens(alg):
    if alg.__name__.startswith("complex_jet"):
        _check_complex_jet_gens(alg)
        return
    if alg.__name__.startswith("jet"):
        _check_jet_gens(alg)
        return
    count = alg.N.bit_length() - 1
    indices = 1 << np.arange(count)
    gens = alg.gens()
    expected = np.eye(alg.N)[indices]
    assert inspect.isbuiltin(alg.gens) and alg.gens.__self__ is alg
    assert isinstance(gens, np.ndarray)
    assert gens.shape == (count,) and gens.dtype == alg.dtype
    assert gens.flags.writeable
    np.testing.assert_array_equal(alg.to_float(gens), expected)
    np.testing.assert_array_equal(alg.to_float(alg.one().gens()), expected)

    if hasattr(alg, "K"):
        squares = alg.K
    elif alg.__name__ in FIXED_METRICS:
        squares = FIXED_METRICS[alg.__name__]
    elif alg.__name__.startswith("dual"):
        squares = [0]*count
    elif alg.__name__.startswith("splitcomplex"):
        squares = [1]*count
    else:
        squares = [-1]*count
    for generator, square in zip(gens, squares):
        np.testing.assert_array_equal(alg.to_float(generator*generator),
                                      alg.to_float(alg(square)))

    # Ascending generator products, with left association for Cayley-Dickson,
    # reproduce every basis element, even when generator squares are zero.
    for index, basis_element in enumerate(alg.basis()):
        value = alg.one()
        for bit, generator in enumerate(gens):
            if index & (1 << bit):
                value = value*generator
        np.testing.assert_array_equal(alg.to_float(value), alg.to_float(basis_element))

    fresh = alg.gens()
    assert not np.shares_memory(gens, fresh)
    if count:
        gens[0] = alg.zero()
        np.testing.assert_array_equal(alg.to_float(fresh), expected)
        np.testing.assert_array_equal(alg.to_float(alg.basis()), np.eye(alg.N))

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


def _check_jet_gens(alg):
    """Jets R[e]/(e^N) have the single generator e; basis element k is e^k."""
    gens = alg.gens()
    assert inspect.isbuiltin(alg.gens) and alg.gens.__self__ is alg
    assert gens.shape == (1,) and gens.dtype == alg.dtype and gens.flags.writeable
    np.testing.assert_array_equal(alg.to_float(gens), np.eye(alg.N)[[1]])
    e = gens[0]
    for k, basis_element in enumerate(alg.basis()):
        np.testing.assert_array_equal(alg.to_float(e.ipow(k)), alg.to_float(basis_element))
    assert e.ipow(alg.N) == alg.zero()


def _check_complex_jet_gens(alg):
    """C[e]/(e^n) over the reals: generators i and e at positions 1 and 2,
    basis element 2k is e^k and 2k + 1 is i*e^k."""
    gens = alg.gens()
    assert inspect.isbuiltin(alg.gens) and alg.gens.__self__ is alg
    assert gens.shape == (2,) and gens.dtype == alg.dtype and gens.flags.writeable
    np.testing.assert_array_equal(alg.to_float(gens), np.eye(alg.N)[[1, 2]])
    i, e = gens
    assert i*i == alg(-1) and i*e == e*i
    for index, basis_element in enumerate(alg.basis()):
        value = e.ipow(index // 2)*(i if index % 2 else alg.one())
        np.testing.assert_array_equal(alg.to_float(value), alg.to_float(basis_element))
    assert e.ipow(alg.N // 2) == alg.zero()
