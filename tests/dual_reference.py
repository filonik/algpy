"""Independent Taylor arithmetic and a compiled C derivative fixture."""
import ctypes
from functools import lru_cache
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from field_reference import FUNCTIONS, reference


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


@pytest.fixture(scope="session")
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
