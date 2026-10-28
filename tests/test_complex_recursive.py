"""Check the recursive specializations against their complex characters."""
import itertools
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy.algebras import complex3, complex4


FUNCTIONS = [
    ("exp", "exp"), ("expm1", "expm1"), ("log", "log"), ("log1p", "log1p"),
    ("sqrt", "sqrt"), ("cos", "cos"), ("sin", "sin"), ("tan", "tan"),
    ("arccos", "acos"), ("arcsin", "asin"), ("arctan", "atan"),
    ("cosh", "cosh"), ("sinh", "sinh"), ("tanh", "tanh"),
    ("arccosh", "acosh"), ("arcsinh", "asinh"), ("arctanh", "atanh"),
]


@pytest.fixture(scope="module", params=[(complex3, 3), (complex4, 4)],
                ids=["complex3", "complex4"])
def algebra(request):
    cls, order = request.param
    # Independent, nonrecursive representation: send i0 to i, and each other
    # generator to +/-i. Basis products are indexed by their generator bits.
    matrix = []
    for signs in itertools.product([-1, 1], repeat=order-1):
        row = []
        for basis in range(cls.N):
            value = 1j**basis.bit_count()
            for generator, sign in enumerate(signs, 1):
                if basis & (1 << generator):
                    value *= sign
            row.append(value)
        matrix.append(row)
    return cls, np.array(matrix)


def reconstruct(matrix, components):
    return (components @ matrix.conj()).real/len(matrix)


@pytest.mark.parametrize("operation,method", FUNCTIONS)
def test_elementary_functions(algebra, operation, method):
    cls, matrix = algebra
    coefficients = list(np.random.default_rng(34).uniform(-0.2, 0.2, (12, cls.N)))
    for real, first_imag, outer_imag in [(2, 0, 0), (-2, 0, 0), (0, 2, 0),
                                         (0, 0, 2), (0, 0, 0)]:
        row = np.zeros(cls.N)
        row[0], row[1], row[cls.N//2] = real, first_imag, outer_imag
        coefficients.append(row)
    # A nonzero zero divisor, whose characters are either zero or two.
    row = np.zeros(cls.N)
    row[0] = row[cls.N//2+1] = 1
    coefficients.append(row)
    coefficients = np.array(coefficients)
    components = coefficients @ matrix.T
    if operation == "log":
        # log is undefined on a zero complex component.
        valid = np.all(components != 0, axis=-1)
        coefficients, components = coefficients[valid], components[valid]
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(components))
    values = cls.from_float(coefficients)
    storage = np.empty(2*len(values), dtype=cls)
    out = storage[::2]
    assert function(values[::-1], out=out) is out
    np.testing.assert_allclose(cls.to_float(out), expected[::-1], rtol=2e-12, atol=2e-14)
    for value, reference in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(getattr(value, method)()), reference,
                                   rtol=2e-12, atol=2e-14)
        np.testing.assert_allclose(cls.to_float(function(value)), reference,
                                   rtol=2e-12, atol=2e-14)


@pytest.mark.parametrize("forward,inverse", [
    ("sin", "arcsin"), ("cos", "arccos"), ("tan", "arctan"),
    ("sinh", "arcsinh"), ("cosh", "arccosh"), ("tanh", "arctanh"),
    ("exp", "log"), ("square", "sqrt"),
])
def test_inverse_functions(algebra, forward, inverse):
    cls, _ = algebra
    value = cls(0.3, *np.linspace(-0.02, 0.02, cls.N-1))
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value),
                               rtol=2e-12, atol=2e-14)


def test_trigonometric_identities(algebra):
    cls, _ = algebra
    value = cls(0.3, *np.linspace(-0.02, 0.02, cls.N-1))
    circular = np.square(np.sin(value)) + np.square(np.cos(value))
    hyperbolic = np.square(np.cosh(value)) - np.square(np.sinh(value))
    np.testing.assert_allclose(cls.to_float(circular), cls.to_float(cls.one()), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(hyperbolic), cls.to_float(cls.one()), atol=2e-14)


def test_power_and_mix(algebra):
    cls, matrix = algebra
    x = cls(1.5, *np.linspace(-0.02, 0.02, cls.N-1))
    y = cls(1.3, *np.linspace(0.01, -0.01, cls.N-1))
    xp, yp = cls.to_float(x) @ matrix.T, cls.to_float(y) @ matrix.T
    expected = reconstruct(matrix, np.power(xp, yp))
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=2e-12, atol=2e-14)
    values = np.array([x, y], dtype=cls)
    np.testing.assert_allclose(cls.to_float(np.power(values, y))[0], expected,
                               rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(x**2), cls.to_float(np.square(x)), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), reconstruct(matrix, np.power(2, xp)),
                               rtol=2e-12, atol=2e-14)
    for t in [0, 0.4, 1]:
        expected = reconstruct(matrix, xp*np.power(yp/xp, t))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected,
                                   rtol=2e-12, atol=2e-14)


def test_arithmetic_and_metadata(algebra):
    cls, matrix = algebra
    rng = np.random.default_rng(35)
    x = cls(1.5, *rng.uniform(-0.02, 0.02, cls.N-1))
    y = cls(1.2, *rng.uniform(-0.01, 0.01, cls.N-1))
    xp, yp = cls.to_float(x) @ matrix.T, cls.to_float(y) @ matrix.T
    for actual, expected in ((x*y, xp*yp), (x/y, xp/yp), (x.inv(), 1/xp),
                             (np.square(x), xp*xp)):
        np.testing.assert_allclose(cls.to_float(actual), reconstruct(matrix, expected),
                                   rtol=2e-12, atol=2e-14)
    assert x.norm() == pytest.approx(np.prod(np.abs(xp))**(1/len(xp)))
    # Preserve the original sign-bit ordering of all generator characters.
    eigenvalues = []
    order = cls.N.bit_length()-1
    for index in range(cls.N):
        positive = index & 1
        row = 0
        for generator in range(1, order):
            bit = (index >> generator) & 1
            row = 2*row + (bit if positive else 1-bit)
        eigenvalues.append(xp[row] if positive else xp[row].conjugate())
    np.testing.assert_allclose(x.eivals(), eigenvalues, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(x.minpoly(), np.polynomial.polynomial.polyfromroots(eigenvalues),
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))


@pytest.mark.parametrize("order", [3, 4])
def test_standalone_header(order, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix, lower = f"complex{order}", f"complex{order-1}"
    dim, half = 2**order, 2**(order-1)
    source = tmp_path/"standalone.c"
    source.write_text(f'''
#include "src/algpy/algebras/multicomplex/{prefix}.h"
#include <assert.h>
#if defined(ALGPY_MULTICOMPLEX1_H) || defined(ALGPY_MULTICOMPLEX2_H) || \\
    defined(ALGPY_MULTICOMPLEX3_H) || defined(ALGPY_MULTICOMPLEX4_H)
#error Concrete complex specializations must not depend on general multicomplex headers
#endif
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[1] = 0.25;
    a.data[{half}] = 0.125;
    {lower} n2 = {prefix}_norm2(a);
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    for (int i = 0; i < {half}; i++) {{
        assert(fabs(product.data[i] - n2.data[i]) < 1e-12);
        assert(fabs(product.data[i+{half}]) < 1e-12);
    }}
    assert(fabs({prefix}_norm{dim}(a) - pow({prefix}_norm(a), {dim})) < 1e-9);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)
