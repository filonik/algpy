"""Independent Cayley-Dickson products and quadratic reference inputs."""
from functools import lru_cache
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy.algebras import general

@lru_cache(None)
def scalar_type(order):
    return getattr(general, f"cayleydickson{order}")(*([-1]*order))


@lru_cache(None)
def basis_coefficient(parameters, i, j):
    """Coefficient of e_(i xor j), from the defining doubling rule."""
    if not parameters:
        return 1
    half = 2**(len(parameters)-1)
    lower = parameters[:-1]
    if i < half:
        return (basis_coefficient(lower, i, j) if j < half else
                basis_coefficient(lower, j-half, i))
    if j < half:
        return (1 if j == 0 else -1)*basis_coefficient(lower, i-half, j)
    return parameters[-1]*(1 if j == half else -1)*basis_coefficient(lower, j-half, i-half)


def multiply(parameters, a, b):
    result = np.zeros(len(a))
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i ^ j] += basis_coefficient(parameters, i, j)*x*y
    return result


def sample(parameters, name, seed=71):
    row = np.random.default_rng(seed).uniform(-0.04, 0.04, 2**len(parameters))
    for i in range(1, len(row)):
        square = basis_coefficient(tuple(parameters), i, i)
        row[i] /= max(1, np.sqrt(abs(square)))
    row[0] = 2 if name == "acosh" else 1.5 if name in {"log", "sqrt", "inv"} else 0.3
    return row


@pytest.fixture(scope="session")
def matrix_functions():
    pytest.importorskip("scipy.linalg")
    from scipy_reference import FUNCTIONS as matrix_references
    return {name: function for name, function in matrix_references}


def check_general_apply_custom_callbacks(order, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    name = f"cayleydickson{order}"
    size = 2**order
    source = tmp_path/'apply.c'
    source.write_text(f'''
#include "src/algpy/algebras/general/cayleydickson/{name}.h"
#include <assert.h>
static double polynomial(double a) {{ return a*a*a + 2*a + 4; }}
static double derivative(double a) {{ return 3*a*a + 2; }}
static complex complex_polynomial(complex a) {{
    return complex_ladd(complex_add(complex_mul(complex_square(a), a),
                                   complex_lmul(a, 2)), 4);
}}
int main(void) {{
    for (int parameter = -1; parameter <= 1; parameter++) {{
        vec{order} k;
        {name} a = {name}_from_re(0.3);
        for (int i = 0; i < {order}; i++) k.data[i] = parameter;
        a.data[1] = 0.25;
        {name} actual = {name}_apply(k, a, polynomial, derivative, complex_polynomial);
        {name} want = {name}_add({name}_mul(k, {name}_square(k, a), a),
                                {name}_ladd({name}_lmul(k, a, 2), 4));
        for (int i = 0; i < {size}; i++) assert(fabs(actual.data[i] - want.data[i]) < 1e-14);
    }}
    return 0;
}}
''')
    binary = tmp_path/'apply'
    subprocess.run([compiler, '-std=c99', '-Wall', '-Wextra', '-Werror', '-I', str(root),
                    str(source), '-lm', '-o', str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)


def check_standalone_specialization_header(name, size, tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    source = tmp_path/'standalone.c'
    source.write_text(f'''
#include "src/algpy/algebras/cayleydickson/{name}.h"
#include <assert.h>
#if defined(ALGPY_CAYLEYDICKSON1_H) || defined(ALGPY_CAYLEYDICKSON2_H) || \\
    defined(ALGPY_CAYLEYDICKSON3_H) || defined(ALGPY_CAYLEYDICKSON4_H)
#error Specializations must not include the general Cayley-Dickson headers
#endif
static complex polynomial(complex a) {{
    return complex_ladd(complex_add(complex_mul(complex_square(a), a),
                                   complex_lmul(a, 2)), 4);
}}
int main(void) {{
    assert(ALGPY_{name.upper()}_N == {size});
    {name} a = {name}_from_re(0.3);
    a.data[{size-1}] = 0.25;
    {name} actual = {name}_apply(a, polynomial);
    {name} want = {name}_add({name}_mul({name}_square(a), a),
                            {name}_ladd({name}_lmul(a, 2), 4));
    for (int i = 0; i < {size}; i++) assert(fabs(actual.data[i] - want.data[i]) < 1e-14);
    return 0;
}}
''')
    binary = tmp_path/'standalone'
    subprocess.run([compiler, '-std=c99', '-Wall', '-Wextra', '-Werror', '-I', str(root),
                    str(source), '-lm', '-o', str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)
