"""Independent character decomposition with nilpotent Taylor coefficients."""
import ctypes
from functools import lru_cache
import itertools
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy import algebras
from algpy.algebras import general
from test_field_derivatives import FUNCTIONS, reference


SIGNATURES = [signs for n in range(1, 5) for signs in itertools.product([-1, 0, 1], repeat=n)]


@lru_cache(None)
def scalar_type(order):
    return getattr(general, f"multicomplex{order}")(*([-1]*order))


def scaled(signs):
    return tuple(sign*magnitude for sign, magnitude in zip(signs, [0.25, 4, 2.25, 9]))


@lru_cache(None)
def representation(parameters):
    n = len(parameters)
    nil = [i for i, k in enumerate(parameters) if k == 0]
    nonzero = [i for i, k in enumerate(parameters) if k != 0]
    negative = [i for i, k in enumerate(parameters) if k < 0]
    primary = negative[0] if negative else None
    choices = [i for i in nonzero if i != primary]
    rows = []
    for signs in itertools.product([-1, 1], repeat=len(choices)):
        eigenvalues = {i: np.sqrt(complex(parameters[i])) for i in nonzero}
        for i, sign in zip(choices, signs):
            eigenvalues[i] *= sign
        matrix = np.zeros((2**len(nil), 2**n), dtype=complex)
        for basis in range(2**n):
            nil_index = sum(1 << j for j, i in enumerate(nil) if basis & (1 << i))
            value = np.prod([eigenvalues[i] for i in nonzero if basis & (1 << i)])
            matrix[nil_index, basis] = value
        rows.append(matrix)
    transform = np.array(rows)
    flattened = transform.reshape(-1, 2**n)
    inverse = np.concatenate([flattened.real, flattened.imag]) if negative else flattened.real
    return transform, inverse, bool(negative)


def nil_product(a, b):
    result = np.zeros(len(a), dtype=complex)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if not i & j:
                result[i | j] += x*y
    return result


def expected(parameters, function, coefficients, offset=0):
    transform, inverse, complex_field = representation(tuple(parameters))
    values = transform @ coefficients
    result = np.zeros_like(values)
    for row, value in enumerate(values):
        nil = value.copy()
        nil[0] = 0
        power = np.zeros(len(nil), dtype=complex)
        power[0] = 1
        for degree in range(len(nil).bit_length()):
            order = degree+offset
            factor = reference(function, value[0], order) if order else function(value[0])
            result[row] += factor*power/math.factorial(degree)
            power = nil_product(power, nil)
    flattened = result.reshape(-1)
    rhs = np.r_[flattened.real, flattened.imag] if complex_field else flattened.real
    return np.linalg.solve(inverse, rhs)


def coefficients(parameters, name, seed=61):
    size = 2**len(parameters)
    result = np.random.default_rng(seed).uniform(-0.02, 0.02, size)
    # Bound each character's contribution despite unequal generator scales.
    for basis in range(size):
        weight = np.prod([math.sqrt(abs(k)) for i, k in enumerate(parameters)
                          if k and basis & (1 << i)])
        result[basis] /= weight
    result[0] = 2 if name == "acosh" else 1.5 if name in ["log", "sqrt", "inv"] else 0.3
    return result


@pytest.fixture(params=SIGNATURES, ids=lambda signs: str(signs))
def algebra(request):
    signs = request.param
    cls = scalar_type(len(signs))
    return cls, signs


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_functions_for_every_signature(algebra, name, function):
    cls, signs = algebra
    for parameters in dict.fromkeys([signs, scaled(signs)]):
        cls.K = parameters
        rows = np.array([coefficients(parameters, name, seed) for seed in [61, 62]])
        wanted = np.array([expected(parameters, function, row) for row in rows])
        values = cls.from_float(rows)
        output = np.empty(2*len(values), dtype=cls)[::2]
        assert function(values[::-1], out=output) is output
        np.testing.assert_allclose(cls.to_float(output), wanted[::-1], rtol=3e-9, atol=3e-10)
        for value, result in zip(values, wanted):
            np.testing.assert_allclose(cls.to_float(getattr(value, name)()), result,
                                       rtol=3e-9, atol=3e-10)
            np.testing.assert_allclose(cls.to_float(function(value)), result, rtol=3e-9, atol=3e-10)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("family,parameter", [("complex", -1), ("splitcomplex", 1), ("dual", 0)])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_specialization_agreement(name, function, family, parameter, order):
    cls = scalar_type(order)
    cls.K = [parameter]*order
    concrete = getattr(algebras, f"{family}{order}")
    row = coefficients(tuple(cls.K), name)
    actual = cls.to_float(function(cls.from_float(row)))
    wanted = concrete.to_float(function(concrete.from_float(row)))
    np.testing.assert_allclose(actual, wanted, rtol=2e-13, atol=2e-14)


def multiply(parameters, a, b):
    result = np.zeros(len(a))
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            weight = np.prod([k for bit, k in enumerate(parameters) if i & j & (1 << bit)])
            result[i ^ j] += weight*x*y
    return result


def test_power_and_mix(algebra):
    cls, signs = algebra
    parameters = scaled(signs)
    cls.K = parameters
    xp = coefficients(parameters, "log")
    yp = coefficients(parameters, "log", 62)
    yp[0] = 1.25
    x, y = cls.from_float(np.array([xp, yp]))
    logarithm = expected(parameters, np.log, xp)
    wanted = expected(parameters, np.exp, multiply(parameters, logarithm, yp))
    np.testing.assert_allclose(cls.to_float(x**y), wanted, rtol=3e-9, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(np.power(np.array([x, y]), y))[0], wanted,
                               rtol=3e-9, atol=3e-10)
    np.testing.assert_allclose(cls.to_float(x**2), multiply(parameters, xp, xp), atol=3e-12)
    wanted = expected(parameters, np.exp, math.log(2)*xp)
    np.testing.assert_allclose(cls.to_float(2**x), wanted, rtol=3e-9, atol=3e-10)
    for t in [0, 0.4, 1]:
        delta = t*(expected(parameters, np.log, yp)-logarithm)
        wanted = multiply(parameters, xp, expected(parameters, np.exp, delta))
        np.testing.assert_allclose(cls.to_float(x.mix(y, t)), wanted, rtol=3e-9, atol=3e-10)


def test_arithmetic_and_metadata(algebra):
    cls, signs = algebra
    parameters = scaled(signs)
    cls.K = parameters
    xp = coefficients(parameters, "log")
    yp = coefficients(parameters, "log", 62)
    x, y = cls.from_float(np.array([xp, yp]))
    for actual, wanted in [(x*y, multiply(parameters, xp, yp)),
                           (x/y, multiply(parameters, xp, expected(parameters, np.reciprocal, yp))),
                           (x.inv(), expected(parameters, np.reciprocal, xp))]:
        np.testing.assert_allclose(cls.to_float(actual), wanted, rtol=3e-9, atol=3e-10)
    eigenvalues = []
    for signs in itertools.product([-1, 1], repeat=len(parameters)):
        factors = [sign*np.sqrt(complex(k)) for sign, k in zip(signs, parameters)]
        eigenvalues.append(sum(xp[i]*np.prod([v for bit, v in enumerate(factors) if i & (1 << bit)])
                               for i in range(cls.N)))
    np.testing.assert_allclose(np.sort_complex(x.eivals()), np.sort_complex(eigenvalues),
                               rtol=2e-13, atol=2e-14)
    polynomial = np.polynomial.polynomial.polyfromroots(eigenvalues)
    np.testing.assert_allclose(x.minpoly(), polynomial, rtol=3e-10, atol=3e-10)
    assert type(x.norm()) is complex
    assert x.norm() == pytest.approx(complex(polynomial[0].real)**(1/cls.N))


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_real_component_powers(order):
    cls = scalar_type(order)
    cls.K = [4, 0.25, 9, 2][:order]
    row = coefficients(tuple(cls.K), "log")
    row[0] = -2
    x = cls.from_float(row)
    np.testing.assert_allclose(cls.to_float(x**2), multiply(cls.K, row, row), atol=3e-13)
    divisor = cls(1, 0.5)
    np.testing.assert_array_equal(cls.to_float(divisor**0), cls.to_float(cls.one()))


@pytest.mark.parametrize("parameters", [(-4,), (0, -4), (4, 0, -4), (0, 4, -4, 0), (0, -4, -2, 0)])
@pytest.mark.parametrize("name,function", [("log", np.log), ("sqrt", np.sqrt), ("asin", np.arcsin),
                                           ("acos", np.arccos), ("acosh", np.arccosh)])
def test_complex_continuation(parameters, name, function):
    cls = scalar_type(len(parameters))
    cls.K = parameters
    row = np.zeros(cls.N)
    row[0] = -2
    row[1 << next(i for i, k in enumerate(parameters) if k < 0)] = 0.4
    for bit, k in enumerate(parameters):
        if k == 0:
            row[1 << bit] = 0.1
    wanted = expected(parameters, function, row)
    np.testing.assert_allclose(cls.to_float(function(cls.from_float(row))), wanted,
                               rtol=3e-9, atol=3e-10)


@pytest.fixture(scope="module")
def c_derivatives(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    directory = tmp_path_factory.mktemp("multicomplex_derivatives")
    root = Path(__file__).resolve().parents[1]
    tables, cases = [], []
    for order in range(1, 4):
        prefix = f"multicomplex{order}"
        rows = ["    {"+", ".join(f"{prefix}_{name}_deriv{d}" for d in range(1, 5-order))+"}"
                for name, _ in FUNCTIONS]
        tables.append(f"static {prefix} (*{prefix}_derivatives[][{4-order}])(vec{order}, {prefix}) = {{\n"
                      +",\n".join(rows)+"\n};")
        cases.append(f'''
    case {order}: {{
        vec{order} k;
        {prefix} a, result;
        for (int i = 0; i < {order}; i++) k.data[i] = parameters[i];
        for (int i = 0; i < {2**order}; i++) a.data[i] = input[i];
        result = {prefix}_derivatives[function][derivative-1](k, a);
        for (int i = 0; i < {2**order}; i++) output[i] = result.data[i];
        break;
    }}''')
    source = directory/"derivatives.c"
    source.write_text('#include "src/algpy/algebras/general/multicomplex/multicomplex3.h"\n'
                      +"\n".join(tables)+'''
void evaluate(int order, int function, int derivative, const double *parameters,
              const double *input, double *output)
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
    library.evaluate.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, pointer, pointer, pointer]
    library.evaluate.restype = None

    def evaluate(parameters, index, derivative, row):
        order = len(parameters)
        storage = ctypes.c_double*(2**order)
        inputs, output = storage(*row), storage()
        k = (ctypes.c_double*order)(*parameters)
        library.evaluate(order, index, derivative, k, inputs, output)
        return np.array(output)

    return evaluate


@pytest.mark.parametrize("order,offset", [(n, d) for n in range(1, 4) for d in range(1, 5-n)])
@pytest.mark.parametrize("index,name,function", [(i, name, f) for i, (name, f) in enumerate(FUNCTIONS)],
                         ids=[name for name, _ in FUNCTIONS])
def test_derivative_helpers(c_derivatives, order, offset, index, name, function):
    for signs in itertools.product([-1, 0, 1], repeat=order):
        parameters = scaled(signs)
        row = coefficients(parameters, name)
        wanted = expected(parameters, function, row, offset)
        actual = c_derivatives(parameters, index, offset, row)
        np.testing.assert_allclose(actual, wanted, rtol=3e-9, atol=3e-10)


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_standalone_header_and_apply(order, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    prefix, lower = f"multicomplex{order}", f"multicomplex{order-1}"
    if order == 1:
        callbacks = '''
static double cube(double x) { return x*x*x; }
static double cube_deriv1(double x) { return 3*x*x; }
static complex cube_complex(complex x) { return complex_mul(complex_square(x), x); }
'''
        arguments = 'cube, cube_deriv1, cube_complex'
    else:
        callbacks = f'''
static {lower} cube(vec{order-1} k, {lower} x) {{ return {lower}_mul(k, {lower}_square(k, x), x); }}
static {lower} cube_deriv1(vec{order-1} k, {lower} x) {{ return vec{2**(order-1)}_scale({lower}_square(k, x), 3); }}
'''
        arguments = 'cube, cube_deriv1'
    source = tmp_path/"standalone.c"
    source.write_text(f'''
#include "src/algpy/algebras/general/multicomplex/{prefix}.h"
#include <assert.h>
{callbacks}
int main(void) {{
    for (int seed = 0; seed < {3**order}; seed++) {{
        int bits = seed;
        vec{order} k;
        for (int i = 0; i < {order}; i++) {{ k.data[i] = (bits%3-1)*4; bits /= 3; }}
        {prefix} a = {prefix}_from_re(1.5);
        for (int i = 1; i < {2**order}; i++) a.data[i] = 0.01*i;
        {prefix} actual = {prefix}_apply(k, a, {arguments});
        {prefix} wanted = {prefix}_mul(k, {prefix}_square(k, a), a);
        for (int i = 0; i < {2**order}; i++) assert(fabs(actual.data[i]-wanted.data[i]) < 1e-11);
        {prefix} two = {prefix}_from_re(2), three = {prefix}_from_re(3);
        assert(fabs({prefix}_rpow(k, 2, three).data[0]-8) < 1e-11);
        assert(fabs({prefix}_lpow(k, two, 3).data[0]-8) < 1e-11);
        assert(fabs({prefix}_mix(k, two, {prefix}_from_re(8), 0.5).data[0]-4) < 1e-11);
    }}
    return 0;
}}
''')
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)
