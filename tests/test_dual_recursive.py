"""Check recursive dual evaluation against finite Taylor series."""
import ctypes
from functools import lru_cache
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy.algebras import dual1, dual2, dual3, dual4
from test_field_derivatives import FUNCTIONS, reference


@pytest.fixture(scope="module", params=[dual1, dual2, dual3, dual4],
                ids=lambda cls: cls.__name__)
def algebra(request):
    return request.param


def product(left, right):
    # Direct subset convolution, independent of the recursive C arithmetic.
    # Overlapping generator bits multiply to zero; disjoint bits concatenate.
    result = np.zeros(len(left))
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if not i & j:
                result[i | j] += x*y
    return result


@lru_cache(None)
def derivative(function, x, order):
    return function(x) if order == 0 else reference(function, complex(x), order).real


def taylor(function, coefficients, offset=0):
    order = len(coefficients).bit_length()-1
    nilpotent = coefficients.copy()
    nilpotent[0] = 0
    power = np.zeros(len(coefficients))
    power[0] = 1
    result = np.zeros(len(coefficients))
    for degree in range(order+1):
        result += derivative(function, coefficients[0], degree+offset)*power/math.factorial(degree)
        power = product(power, nilpotent)
    return result


def real_points(name):
    if name in ["log", "sqrt", "inv"]:
        return [0.5, 1.25, 2]
    if name == "acosh":
        return [1.5, 2, 3]
    return [-0.4, 0, 0.4]


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_elementary_functions(algebra, name, function):
    cls = algebra
    coefficients = np.random.default_rng(51).uniform(-0.3, 0.3, (3, cls.N))
    coefficients[:, 0] = real_points(name)
    expected = np.array([taylor(function, row) for row in coefficients])
    values = cls.from_float(coefficients)
    storage = np.empty(2*len(values), dtype=cls)
    out = storage[::2]
    assert function(values[::-1], out=out) is out
    np.testing.assert_allclose(cls.to_float(out), expected[::-1], rtol=3e-10, atol=3e-10)
    for value, result in zip(values, expected):
        np.testing.assert_allclose(cls.to_float(getattr(value, name)()), result,
                                   rtol=3e-10, atol=3e-10)
        np.testing.assert_allclose(cls.to_float(function(value)), result, rtol=3e-10, atol=3e-10)


@pytest.fixture(scope="module")
def c_derivatives(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    directory = tmp_path_factory.mktemp("dual_derivatives")
    root = Path(__file__).resolve().parents[1]
    tables, cases = [], []
    for order in range(1, 4):
        prefix = f"dual{order}"
        rows = ["    {"+", ".join(f"{prefix}_{name}_deriv{d}" for d in range(1, 5-order))+"}"
                for name, _ in FUNCTIONS]
        tables.append(f"static {prefix} (*{prefix}_derivatives[][{4-order}])({prefix}) = {{\n"
                      +",\n".join(rows)+"\n};")
        cases.append(f'''
    case {order}: {{
        {prefix} a, result;
        for (int i = 0; i < {2**order}; i++) a.data[i] = input[i];
        result = {prefix}_derivatives[function][derivative-1](a);
        for (int i = 0; i < {2**order}; i++) output[i] = result.data[i];
        break;
    }}''')
    source = directory/"derivatives.c"
    source.write_text('''
#include "src/algpy/algebras/multicomplex/dual3.h"
'''+"\n".join(tables)+'''
void evaluate(int order, int function, int derivative, const double *input, double *output)
{
    switch (order) {
'''+"\n".join(cases)+'''
    }
}
''')
    path = directory/"derivatives.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-shared",
                    "-fPIC", "-I", str(root), str(source), "-lm", "-o", str(path)],
                   check=True, capture_output=True, text=True)
    library = ctypes.CDLL(str(path))
    pointer = ctypes.POINTER(ctypes.c_double)
    library.evaluate.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, pointer, pointer]
    library.evaluate.restype = None

    def evaluate(order, index, derivative, coefficients):
        storage = ctypes.c_double*(2**order)
        inputs, outputs = storage(*coefficients), storage()
        library.evaluate(order, index, derivative, inputs, outputs)
        return np.array(outputs)

    return evaluate


@pytest.mark.parametrize("order,offset", [(n, d) for n in range(1, 4) for d in range(1, 5-n)])
@pytest.mark.parametrize("index,name,function", [(i, name, function)
                         for i, (name, function) in enumerate(FUNCTIONS)],
                         ids=[name for name, _ in FUNCTIONS])
def test_derivative_helpers(c_derivatives, order, offset, index, name, function):
    coefficients = np.random.default_rng(52).uniform(-0.3, 0.3, 2**order)
    coefficients[0] = real_points(name)[0]
    expected = taylor(function, coefficients, offset)
    actual = c_derivatives(order, index, offset, coefficients)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-10)


def test_pure_nilpotents(algebra):
    cls = algebra
    coefficients = np.zeros(cls.N)
    for generator in range(cls.N.bit_length()-1):
        coefficients[1 << generator] = 1
    value = cls.from_float(coefficients)
    # Every subset occurs degree! times in the corresponding power.
    np.testing.assert_array_equal(cls.to_float(np.exp(value)), np.ones(cls.N))
    expected = [0]+[(-1)**(i.bit_count()-1)*math.factorial(i.bit_count()-1)
                    for i in range(1, cls.N)]
    np.testing.assert_array_equal(cls.to_float(np.log1p(value)), expected)


@pytest.mark.parametrize("function", [np.expm1, np.log1p])
def test_small_inputs(algebra, function):
    cls = algebra
    coefficients = np.linspace(-1, 1, cls.N)*1e-20
    coefficients[0] = 1e-20
    result = function(cls.from_float(coefficients))
    np.testing.assert_allclose(cls.to_float(result), coefficients, rtol=2e-15, atol=1e-35)


@pytest.mark.parametrize("forward,inverse", [
    (np.sin, np.arcsin), (np.cos, np.arccos), (np.tan, np.arctan),
    (np.sinh, np.arcsinh), (np.cosh, np.arccosh), (np.tanh, np.arctanh),
    (np.exp, np.log), (np.square, np.sqrt),
])
def test_inverse_functions(algebra, forward, inverse):
    cls = algebra
    value = cls(1.5 if inverse is np.arccosh else 0.4,
                *np.linspace(-0.1, 0.1, cls.N-1))
    result = forward(inverse(value))
    np.testing.assert_allclose(cls.to_float(result), cls.to_float(value), rtol=2e-12, atol=2e-14)


def test_power_and_mix(algebra):
    cls = algebra
    xp = np.r_[1.5, np.linspace(-0.1, 0.1, cls.N-1)]
    yp = np.r_[1.3, np.linspace(0.1, -0.1, cls.N-1)]
    x, y = cls.from_float(np.array([xp, yp]))
    expected = taylor(np.exp, product(taylor(np.log, xp), yp))
    np.testing.assert_allclose(cls.to_float(x**y), expected, rtol=3e-10, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], expected,
                               rtol=3e-10, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(x**2), product(xp, xp), rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cls.to_float(2**x), taylor(np.exp, math.log(2)*xp),
                               rtol=3e-10, atol=3e-10)
    for t in [0, 0.4, 1]:
        expected = product(xp, taylor(np.exp, t*(taylor(np.log, yp)-taylor(np.log, xp))))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), expected, rtol=3e-10, atol=3e-10)


def test_arithmetic_and_metadata(algebra):
    cls = algebra
    xp = np.r_[1.5, np.random.default_rng(53).uniform(-0.2, 0.2, cls.N-1)]
    yp = np.r_[1.25, np.random.default_rng(54).uniform(-0.2, 0.2, cls.N-1)]
    x, y = cls.from_float(np.array([xp, yp]))
    for actual, expected in [(x*y, product(xp, yp)),
                             (x/y, product(xp, taylor(np.reciprocal, yp))),
                             (x.inv(), taylor(np.reciprocal, xp))]:
        np.testing.assert_allclose(cls.to_float(actual), expected, rtol=3e-10, atol=3e-10)
    for real in [0, -0.75, 1.5]:
        value = cls(real, *xp[1:])
        assert type(value.norm()) is float
        assert value.norm() == abs(real)
        np.testing.assert_array_equal(value.eivals(), np.full(cls.N, real, dtype=complex))
        expected = np.polynomial.polynomial.polyfromroots([real]*cls.N)
        np.testing.assert_allclose(value.minpoly(), expected, rtol=2e-13, atol=1e-12)
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))
    expected = np.zeros((cls.N, cls.N, cls.N))
    for i in range(cls.N):
        for j in range(cls.N):
            if not i & j:
                expected[i, j, i | j] = 1
    np.testing.assert_array_equal(cls.to_float(cls.coeffs()), expected)


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_standalone_header_and_apply(order, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix = f"dual{order}"
    lower = f"dual{order-1}" if order > 1 else "double"
    half = 2**(order-1)
    cube = 'return a*a*a;' if order == 1 else f'return {lower}_mul({lower}_square(a), a);'
    cube_deriv = 'return 3*a*a;' if order == 1 else f'return {lower}_lmul({lower}_square(a), 3);'
    norm_check = (f'assert(fabs(product.data[0] - {prefix}_norm2(a)) < 1e-12);'
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
#error Concrete dual specializations must not depend on general multicomplex headers
#endif
static inline {lower} cube({lower} a) {{ {cube} }}
static inline {lower} cube_deriv1({lower} a) {{ {cube_deriv} }}
int main(void) {{
    {prefix} a = {prefix}_from_re(1.5);
    a.data[1] = 0.25;
    a.data[{half}] += 0.125;
    {prefix} applied = {prefix}_apply(a, cube, cube_deriv1);
    {prefix} expected = {prefix}_mul({prefix}_square(a), a);
    for (int i = 0; i < {2**order}; i++) {{
        assert(fabs(applied.data[i] - expected.data[i]) < 1e-12);
    }}
    {prefix} product = {prefix}_mul(a, {prefix}_conjugate(a));
    {norm_check}
    assert(fabs({prefix}_norm{2**order}(a) - pow({prefix}_norm(a), {2**order})) < 1e-9);
    {prefix} two = {prefix}_from_re(2), three = {prefix}_from_re(3);
    assert(fabs({prefix}_rpow(2, three).data[0] - 8) < 1e-12);
    assert(fabs({prefix}_lpow(two, 3).data[0] - 8) < 1e-12);
    assert(fabs({prefix}_mix(two, {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12);
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)
