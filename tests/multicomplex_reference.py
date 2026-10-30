"""Independent multicomplex character and nilpotent reference calculations."""
import ctypes
from functools import lru_cache
import itertools
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy.algebras import general
from field_reference import FUNCTIONS, reference

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


def multiply(parameters, a, b):
    result = np.zeros(len(a))
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            weight = np.prod([k for bit, k in enumerate(parameters) if i & j & (1 << bit)])
            result[i ^ j] += weight*x*y
    return result


@pytest.fixture(scope="session")
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


def check_standalone_header_and_apply(order, tmp_path):
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
        assert(fabs({prefix}_rpows(k, 2, three).data[0]-8) < 1e-11);
        assert(fabs({prefix}_lpows(k, two, 3).data[0]-8) < 1e-11);
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
