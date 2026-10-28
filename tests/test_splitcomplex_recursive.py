"""Check split-complex specializations against independent real characters."""
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy.algebras import splitcomplex1, splitcomplex2, splitcomplex3, splitcomplex4


FUNCTIONS = [
    ("exp", "exp", -1, 1), ("expm1", "expm1", -1, 1),
    ("log", "log", 0.2, 2), ("log1p", "log1p", -0.8, 1),
    ("sqrt", "sqrt", 0.2, 2),
    ("cos", "cos", -1, 1), ("sin", "sin", -1, 1), ("tan", "tan", -1, 1),
    ("arccos", "acos", -0.9, 0.9), ("arcsin", "asin", -0.9, 0.9),
    ("arctan", "atan", -2, 2),
    ("cosh", "cosh", -1, 1), ("sinh", "sinh", -1, 1), ("tanh", "tanh", -1, 1),
    ("arccosh", "acosh", 1.1, 3), ("arcsinh", "asinh", -2, 2),
    ("arctanh", "atanh", -0.9, 0.9),
]


@pytest.fixture(scope="module", params=[splitcomplex1, splitcomplex2, splitcomplex3,
                                        splitcomplex4], ids=lambda cls: cls.__name__)
def algebra(request):
    cls = request.param
    # Each generator independently maps to -1 or +1. Evaluate each basis
    # product directly, without the recursive splitting used by the C code.
    order = cls.N.bit_length()-1
    matrix = []
    for index in range(cls.N):
        row = []
        for basis in range(cls.N):
            sign = 1
            for generator in range(order):
                if basis & (1 << generator):
                    sign *= 1 if index & (1 << generator) else -1
            row.append(sign)
        matrix.append(row)
    return cls, np.array(matrix)


def reconstruct(matrix, components):
    return components @ matrix/len(matrix)


@pytest.mark.parametrize("operation,method,low,high", FUNCTIONS)
def test_elementary_functions(algebra, operation, method, low, high):
    cls, matrix = algebra
    components = np.random.default_rng(41).uniform(low, high, (12, cls.N))
    coefficients = reconstruct(matrix, components)
    # Compute the reference from the actual input coefficients after rounding.
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(coefficients @ matrix.T))
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


@pytest.mark.parametrize("operation", ["exp", "expm1", "sqrt", "sin", "cos", "tanh"])
def test_zero_divisors(algebra, operation):
    cls, matrix = algebra
    # A nonzero idempotent, with real characters zero and one.
    x = cls(0.5, 0.5)
    components = cls.to_float(x) @ matrix.T
    function = getattr(np, operation)
    expected = reconstruct(matrix, function(components))
    np.testing.assert_allclose(cls.to_float(function(x)), expected, atol=2e-14)


@pytest.mark.parametrize("operation,value", [
    ("sqrt", -1), ("log", -1), ("log1p", -2),
    ("arcsin", 2), ("arccos", 2), ("arccosh", 0), ("arctanh", 2),
])
def test_real_function_domains(algebra, operation, value):
    cls, _ = algebra
    with np.errstate(all="ignore"):
        result = getattr(np, operation)(cls(value))
    assert np.isnan(cls.to_float(result)).all()


@pytest.mark.parametrize("forward,inverse,components", [
    ("sin", "arcsin", [-1, -0.5, 0.5, 1]),
    ("cos", "arccos", [-1, -0.5, 0.5, 1]),
    ("tan", "arctan", [-2, -0.5, 0.5, 2]),
    ("sinh", "arcsinh", [-2, -0.5, 0.5, 2]),
    ("cosh", "arccosh", [1, 1.5, 2, 3]),
    ("tanh", "arctanh", [-0.75, -0.5, 0.5, 0.75]),
    ("exp", "log", [0.25, 0.5, 1, 2]),
    ("square", "sqrt", [0, 0.5, 1, 2]),
])
def test_inverse_functions(algebra, forward, inverse, components):
    cls, matrix = algebra
    value = cls.from_float(reconstruct(matrix, np.resize(components, cls.N)))
    result = getattr(np, forward)(getattr(np, inverse)(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value), atol=2e-14)


def test_power_and_mix(algebra):
    cls, matrix = algebra
    rng = np.random.default_rng(42)
    xp = rng.uniform(0.5, 2, cls.N)
    yp = rng.uniform(0.5, 2, cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    expected = reconstruct(matrix, xp**yp)
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], expected,
                               rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), reconstruct(matrix, 2**xp), atol=2e-14)
    for t in [0, 0.4, 1]:
        expected = reconstruct(matrix, xp*(yp/xp)**t)
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected,
                                   rtol=2e-12, atol=2e-14)
    # Real pow handles negative components with integer exponents and zero**0.
    xp = np.resize([-2, -1, 0, 2], cls.N)
    yp = np.resize([2, 3, 0, 4], cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    np.testing.assert_allclose(cls.to_float(x**y), reconstruct(matrix, xp**yp), atol=2e-14)
    np.testing.assert_allclose(cls.to_float(x**2), cls.to_float(np.square(x)), atol=2e-14)
    np.testing.assert_allclose(cls.to_float((-2)**y), reconstruct(matrix, (-2)**yp),
                               atol=2e-14)


def test_arithmetic_and_metadata(algebra):
    cls, matrix = algebra
    xp = np.random.default_rng(43).uniform(0.5, 2, cls.N)
    xp[0] = -xp[0]
    yp = np.random.default_rng(44).uniform(0.5, 2, cls.N)
    x, y = cls.from_float(reconstruct(matrix, np.array([xp, yp])))
    for actual, expected in ((x*y, xp*yp), (x/y, xp/yp), (x.inv(), 1/xp),
                             (np.square(x), xp*xp)):
        np.testing.assert_allclose(cls.to_float(actual), reconstruct(matrix, expected),
                                   rtol=2e-12, atol=2e-14)
    assert type(x.norm()) is complex
    assert x.norm() == pytest.approx(complex(np.prod(xp))**(1/cls.N))
    # The original four-component API uses its own character ordering.
    eigenvalues = xp[[3, 1, 2, 0]] if cls.N == 4 else xp
    np.testing.assert_allclose(x.eivals(), eigenvalues, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(x.minpoly(), np.polynomial.polynomial.polyfromroots(xp),
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))
    # Products of basis elements are indexed by the XOR of their bit masks.
    expected = np.eye(cls.N)[np.arange(cls.N)[:, None] ^ np.arange(cls.N)]
    np.testing.assert_array_equal(cls.to_float(cls.coeffs()), expected)


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_standalone_header_and_apply(order, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix, lower = f"splitcomplex{order}", f"splitcomplex{order-1}"
    half = 2**(order-1)
    callback = "square_real" if order == 1 else f"{lower}_square"
    norm_check = (f"assert(fabs(product.data[0] - {prefix}_norm2(a)) < 1e-12);"
                  if order == 1 else f'''
    {lower} n2 = {prefix}_norm2(a);
    for (int i = 0; i < {half}; i++) {{
        assert(fabs(product.data[i] - n2.data[i]) < 1e-12);
        assert(fabs(product.data[i+{half}]) < 1e-12);
    }}''')
    source = tmp_path/"standalone.c"
    source.write_text(f'''
#include "src/algpy/algebras/multicomplex/{prefix}.h"
#include <assert.h>
#if defined(ALGPY_MULTICOMPLEX1_H) || defined(ALGPY_MULTICOMPLEX2_H) || \\
    defined(ALGPY_MULTICOMPLEX3_H) || defined(ALGPY_MULTICOMPLEX4_H)
#error Concrete split-complex specializations must not depend on general multicomplex headers
#endif
static inline double square_real(double a) {{ return a*a; }}
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[{half}] = 0.25;
    {prefix} applied = {prefix}_apply(a, {callback});
    {prefix} squared = {prefix}_square(a);
    for (int i = 0; i < {2**order}; i++) {{
        assert(fabs(applied.data[i] - squared.data[i]) < 1e-12);
    }}
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    {norm_check}
    {prefix} two = {prefix}_from_re(2);
    {prefix} three = {prefix}_from_re(3);
    assert({prefix}_rpow(2, three).data[0] == 8);
    assert({prefix}_lpow(two, 3).data[0] == 8);
    assert(fabs({prefix}_mix(two, {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)
