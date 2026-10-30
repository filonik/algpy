"""The same C entry points and NumPy operations exist for every algebra."""




import numpy as np

import pytest

from algpy import algebras

from algpy.algebras import general
from algebra_reference import CONCRETE




def assert_components(cls, a, b):
    np.testing.assert_allclose(cls.to_float(a), cls.to_float(b), atol=1e-12, equal_nan=True)

def check_unary_scalar_array_and_method(algebra, name, method):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    # Include a strided input and an explicitly supplied output buffer.
    values = np.array([x, x, x], dtype=cls)[::2]
    out = np.empty_like(values)
    with np.errstate(all="ignore"):
        result = getattr(np, name)(values, out=out)
        scalar = getattr(np, name)(x)
        assert result is out and result.dtype == cls.dtype
        assert_components(cls, result[0], scalar)
        assert_components(cls, result[1], scalar)
        if method:
            assert_components(cls, getattr(x, method)(), scalar)

def check_binary_cast_broadcast_and_operator(algebra, name):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    values = np.array([[x], [x]], dtype=cls)
    with np.errstate(all="ignore"):
        result = getattr(np, name)(values, np.array([1.0, 2.0]))
        assert result.shape == (2, 2) and result.dtype == cls.dtype
        scalar = {"add": lambda: x+2, "subtract": lambda: x-2,
                  "multiply": lambda: x*2, "divide": lambda: x/2,
                  "power": lambda: x**2}[name]()
        assert_components(cls, result[1, 1], scalar)
        if name == "power":
            assert_components(cls, 2**x, np.power(cls(2), x))

def check_predicates_parts_and_inverse(algebra):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    values = np.array([x, cls.zero()], dtype=cls)
    assert np.equal(values, x).tolist() == [True, False]
    assert np.not_equal(values, x).tolist() == [False, True]
    assert_components(cls, np.real(x) + np.imag(x), x)
    assert_components(cls, np.real(x), x.re())
    assert_components(cls, np.imag(x), x.im())
    assert x.to_re() == pytest.approx(0.1)
    assert x.is_real() == (cls.N == 1)
    assert cls.from_re(2).is_real()
    assert_components(cls, cls.from_re(2), cls(2))
    # Shift away from zero divisors in the split-complex specializations.
    invertible = x + 4
    assert_components(cls, invertible.inv()*invertible, cls.one())
    with np.errstate(all="ignore"):
        assert_components(cls, invertible.mix(2, 0.5),
                          invertible * np.power(invertible.inv()*2, 0.5))

def check_norm_static_return_type(algebra):
    cls = algebra
    norm_type = complex if hasattr(cls, "K") or cls.__name__.startswith("splitcomplex") else float
    for value, expected in ((cls.zero(), 0), (cls.one(), 1), (cls(-2), 2)):
        norm = value.norm()
        assert type(norm) is norm_type
        assert norm == pytest.approx(expected)
    value = cls(2, *(np.arange(1, cls.N)/100))
    norm = value.norm()
    assert type(norm) is norm_type
    if cls.N == 1:
        assert norm == 2
    else:
        polynomial = value.minpoly()
        expected = complex(polynomial[0])**(1/(len(polynomial)-1))
        assert norm == pytest.approx(expected)

def check_general_norm_type_does_not_follow_parameters(family, n):
    cls = getattr(general, f"{family}{n}")(*([-1]*n))
    generator = cls.basis()[1]
    assert type(generator.norm()) is complex
    assert generator.norm() == pytest.approx(1)
    cls.K = [1]*n
    norm = generator.norm()
    assert type(norm) is complex
    # Cayley-Dickson has a quadratic norm; the other families' higher
    # even-degree norms can remain real for this particular generator.
    if family == "cayleydickson" or n == 1:
        assert norm == pytest.approx(1j)
    assert type(cls.one().norm()) is complex

def check_splitcomplex_negative_norm(n):
    cls = getattr(algebras, f"splitcomplex{n}")
    # One character/eigenvalue is -1, the others +1: determinant -1.
    projector = cls.one()
    for generator in range(n):
        projector = projector*(cls.one()+cls.basis()[1 << generator])/2
    value = cls.one() - 2*projector
    norm = value.norm()
    assert type(norm) is complex
    assert norm == pytest.approx(np.exp(1j*np.pi/cls.N))

def check_clifford_negative_norm(n, parameters, central_index):
    cls = getattr(general, f"clifford{n}")(*parameters)
    basis = cls.basis()
    projector = (cls.one()+basis[central_index])/2
    if n == 3:
        value = basis[1]*projector + (cls.one()-projector)
    else:
        projector = projector*(cls.one()+basis[1])/2
        value = cls.one()-2*projector
    assert value.norm() == pytest.approx(np.exp(1j*np.pi/4))
    assert type(value.norm()) is complex

def check_specialization_structure(name, n):
    cls = getattr(algebras, name)
    basis, coeffs = cls.basis(), cls.coeffs()
    assert basis.shape == (n,) and coeffs.shape == (n, n)
    np.testing.assert_array_equal(cls.to_float(basis), np.eye(n))
    for i in range(n):
        for j in range(n):
            assert_components(cls, coeffs[i, j], basis[i]*basis[j])
    x = cls(*np.linspace(0.1, 0.3, n))
    np.testing.assert_allclose(
        np.polynomial.Polynomial.fromroots(x.eivals()).coef,
        x.minpoly(), atol=1e-10,
    )

def check_parameterized_ufunc_reads_current_parameters():
    cls = general.cayleydickson2(-1, -1)
    x = cls(1, 0.2, 0.1, 0.05)
    values = np.array([x], dtype=cls)
    before = np.exp(values)
    cls.K = [1, 0]
    assert_components(cls, np.exp(values)[0], x.exp())
    assert not np.allclose(cls.to_float(before), cls.to_float(np.exp(values)))
    assert_components(cls, np.reciprocal(values)[0], x.inv())
