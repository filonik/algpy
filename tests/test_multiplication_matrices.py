"""Left matrices act on coefficient columns, right matrices on rows."""
from functools import lru_cache
import itertools

import numpy as np
import pytest

from algpy import algebras
from algpy.algebras import general
from test_api_consistency import CONCRETE, FAMILIES


@lru_cache(None)
def parameterized_type(family, order):
    return getattr(general, f"{family}{order}")(*([-1]*order))


@pytest.fixture(params=[
    *((family, n, True) for family in FAMILIES for n in range(1, 5)),
    *((name, n, False) for name, n in CONCRETE),
], ids=lambda p: f"{p[0]}{p[1] if p[2] else ''}")
def algebra(request):
    name, size, parameterized = request.param
    if parameterized:
        cls = parameterized_type(name, size)
        cls.K = [-1]*size
        return cls
    return getattr(algebras, name)


def assert_action(cls, a, b):
    structure = cls.to_float(cls.coeffs())
    left, right = a.lmat(), b.rmat()
    np.testing.assert_allclose(
        left, np.einsum("i,ijk->kj", cls.to_float(a), structure),
        rtol=2e-12, atol=2e-12,
    )
    np.testing.assert_allclose(
        right, np.einsum("j,ijk->ik", cls.to_float(b), structure),
        rtol=2e-12, atol=2e-12,
    )
    np.testing.assert_allclose(left @ cls.to_float(b), cls.to_float(a*b),
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(cls.to_float(a) @ right, cls.to_float(a*b),
                               rtol=2e-12, atol=2e-12)


def test_matrix_shape_dtype_and_action(algebra):
    cls = algebra
    rng = np.random.default_rng(71)
    a = cls(*rng.uniform(-1, 1, cls.N))
    b = cls(*rng.uniform(-1, 1, cls.N))
    for value in [cls.zero(), cls.one(), a]:
        for matrix in [value.lmat(), value.rmat()]:
            assert isinstance(matrix, np.ndarray)
            assert matrix.shape == (cls.N, cls.N)
            assert matrix.dtype == np.dtype(np.float64)
            assert matrix.flags.c_contiguous and matrix.flags.owndata
            assert matrix.flags.writeable
        assert_action(cls, value, b)
    np.testing.assert_array_equal(cls.zero().lmat(), np.zeros((cls.N, cls.N)))
    np.testing.assert_array_equal(cls.zero().rmat(), np.zeros((cls.N, cls.N)))
    np.testing.assert_array_equal(cls.one().lmat(), np.eye(cls.N))
    np.testing.assert_array_equal(cls.one().rmat(), np.eye(cls.N))


def test_linearity_and_independent_results(algebra):
    cls = algebra
    rng = np.random.default_rng(72)
    a = cls(*rng.uniform(-1, 1, cls.N))
    b = cls(*rng.uniform(-1, 1, cls.N))
    for name in ["lmat", "rmat"]:
        matrix = getattr(a, name)()
        other = getattr(a, name)()
        assert not np.shares_memory(matrix, other)
        np.testing.assert_allclose(getattr(a+b, name)(), matrix+getattr(b, name)(),
                                   rtol=2e-13, atol=2e-13)
        matrix[:] = 123
        np.testing.assert_array_equal(getattr(a, name)(), other)
        with pytest.raises(TypeError):
            getattr(a, name)(b)


PARAMETERS = [k for n in range(1, 5) for k in itertools.product([-2, 0, 0.5], repeat=n)]


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("parameters", PARAMETERS, ids=str)
def test_general_parameters(family, parameters):
    cls = parameterized_type(family, len(parameters))
    cls.K = parameters
    rng = np.random.default_rng(73)
    a = cls(*rng.uniform(-1, 1, cls.N))
    b = cls(*rng.uniform(-1, 1, cls.N))
    assert_action(cls, a, b)


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_matrices_read_current_parameters(family, order):
    cls = parameterized_type(family, order)
    cls.K = [-1]*order
    a = cls.basis()[1]
    before = [a.lmat(), a.rmat()]
    cls.K = [0.5]+[0]*(order-1)
    for matrix, old, index in zip([a.lmat(), a.rmat()], before, [(0, 1), (1, 0)]):
        assert matrix[index] == 0.5
        assert old[index] == -1
    assert_action(cls, a, cls(*np.linspace(-0.5, 0.5, cls.N)))


def test_quaternion_orientation():
    q = algebras.quaternion(1, 2, 3, 4)
    np.testing.assert_array_equal(q.lmat(), [
        [1, -2, -3, -4], [2, 1, -4, 3], [3, 4, 1, -2], [4, -3, 2, 1],
    ])
    np.testing.assert_array_equal(q.rmat(), [
        [1, 2, 3, 4], [-2, 1, -4, 3], [-3, 4, 1, -2], [-4, -3, 2, 1],
    ])


@pytest.mark.parametrize("cls,expected", [
    (algebras.real, [[2]]),
    (algebras.complex1, [[2, -3], [3, 2]]),
    (algebras.splitcomplex1, [[2, 3], [3, 2]]),
    (algebras.dual1, [[2, 0], [3, 2]]),
])
def test_small_commutative_matrices(cls, expected):
    value = cls(2) if cls.N == 1 else cls(2, 3)
    np.testing.assert_array_equal(value.lmat(), expected)
    np.testing.assert_array_equal(value.rmat(), np.array(expected).T)


@pytest.mark.parametrize("cls", [algebras.octonion, algebras.sedenion])
def test_nonassociative_products(cls):
    basis = cls.basis()
    a, b, c = basis[1], basis[2], basis[4]
    direct = (a*b).lmat() @ cls.to_float(c)
    nested = a.lmat() @ b.lmat() @ cls.to_float(c)
    np.testing.assert_array_equal(direct, cls.to_float((a*b)*c))
    np.testing.assert_array_equal(nested, cls.to_float(a*(b*c)))
    assert not np.array_equal(direct, nested)
    assert_action(cls, a*b, c)
