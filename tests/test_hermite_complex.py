"""Complex quartics and Hermite calculus, independent of any algebra."""
import ctypes
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from field_reference import FUNCTIONS, center, operation


@pytest.fixture(scope="module")
def complex_engine(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    directory = tmp_path_factory.mktemp("hermite_complex")
    cases = []
    for index, (name, _) in enumerate(FUNCTIONS):
        callbacks = (f"double_{name}, complex_{name}, complex_{name}_deriv1, "
                     f"complex_{name}_deriv2, complex_{name}_deriv3, result")
        cases.append(f"case {index}: ok = solve ? poly4_complex_apply(p, w, {callbacks}) : "
                     f"hermite_complex_interpolate(n, w, nodes, {callbacks}); break;")
    source = directory/"engine.c"
    source.write_text('''
#include "src/algpy/polynomials/hermite_complex.h"
int evaluate(int function, int solve, int n, const double *offset,
             const double *input, const double *coefficients, double *output) {
    complex w = {{offset[0], offset[1]}}, nodes[4], p[5], result[4];
    int ok;
    for (int i = 0; i < n; i++) nodes[i] = complex_make(input + 2*i);
    for (int i = 0; i < 5; i++) p[i] = complex_make(coefficients + 2*i);
    switch (function) {
'''+"\n".join(cases)+'''
    default: return 0;
    }
    for (int i = 0; i < 4; i++) {
        output[2*i] = result[i].data[0]; output[2*i+1] = result[i].data[1];
    }
    return ok;
}
void roots(int n, const double *input, double *output) {
    complex p[5] = {0}, result[4];
    if (n < 2 || n > 4) return;
    for (int i = 0; i <= n; i++) p[i] = complex_make(input + 2*i);
    if (n == 2) poly2_complex_coeffs_to_roots(p, result);
    else if (n == 3) poly3_complex_coeffs_to_roots(p, result);
    else poly4_complex_coeffs_to_roots(p, result);
    for (int i = 0; i < n; i++) {
        output[2*i] = result[i].data[0]; output[2*i+1] = result[i].data[1];
    }
}
''')
    library_path = directory/"engine.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-shared", "-fPIC", "-I", str(root), str(source), "-lm", "-o", str(library_path)],
                   check=True, capture_output=True, text=True)
    library = ctypes.CDLL(str(library_path))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    library.evaluate.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, pointer, pointer, pointer, pointer]
    library.evaluate.restype = ctypes.c_int
    library.roots.argtypes = [ctypes.c_int, pointer, pointer]
    library.roots.restype = None
    return library


def interpolate(engine, name, nodes, solve):
    nodes = np.ascontiguousarray(nodes, dtype=np.complex128)
    p = np.zeros(5, dtype=np.complex128)
    p[:len(nodes)+1] = np.polynomial.polynomial.polyfromroots(nodes)
    output = np.empty(4, dtype=np.complex128)
    w = np.array([center(name) + .25j], dtype=np.complex128)
    success = engine.evaluate([name for name, _ in FUNCTIONS].index(name), solve, len(nodes),
                             w.view(float), nodes.view(float), p.view(float), output.view(float))
    assert success, (name, nodes, solve)
    return w[0], output


NODES = [
    [-.2+.1j, -.1-.05j, .1+.03j, .2+.07j],
    [0]*4, [.1+.2j]*4, [-.1+.2j]*2 + [.2-.1j]*2,
    [.1+.2j]*3 + [-.3-.1j], [0, 0, 0], [.1j, .1j], [.2j],
]


def check_complex_matrix(complex_engine, nodes, name):
    from scipy_reference import FUNCTIONS as references
    p = np.polynomial.polynomial.polyfromroots(nodes)
    matrix = np.diag(np.ones(len(nodes)-1), -1).astype(complex)
    matrix[:, -1] = -p[:-1]
    identity = np.eye(len(nodes))
    for solve in ([False, True] if len(nodes) == 4 else [False]):
        w, coefficients = interpolate(complex_engine, name, nodes, solve)
        actual = np.zeros_like(matrix)
        power = identity
        for coefficient in coefficients:
            actual += coefficient*power
            power = power@matrix
        expected = dict(references)[operation(name)](w*identity + matrix)
        np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)


@pytest.mark.parametrize("nodes", NODES, ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_complex_hermite(complex_engine, nodes, name, function):
    check_complex_matrix(complex_engine, nodes, name)


@pytest.mark.parametrize("delta", [1e-2, 1e-3, 1e-6, 1e-12, 1e-40, 0])
@pytest.mark.parametrize("name", [name for name, _ in FUNCTIONS])
def test_complex_nearly_confluent(complex_engine, name, delta):
    nodes = delta*np.array([-.2+.1j, -.1-.05j, .1+.03j, .2+.07j])
    check_complex_matrix(complex_engine, nodes, name)


@pytest.mark.parametrize("nodes", NODES[:5] + [[.1+.2j, -.3+.1j], [.1+.2j]*3], ids=str)
def test_complex_polynomial_roots(complex_engine, nodes):
    from itertools import permutations
    nodes = np.asarray(nodes, dtype=np.complex128)
    p = np.ascontiguousarray(np.polynomial.polynomial.polyfromroots(nodes))
    output = np.empty(len(nodes), dtype=np.complex128)
    complex_engine.roots(len(nodes), p.view(float), output.view(float))
    error = min(np.max(np.abs(output - candidate)) for candidate in permutations(nodes))
    # Coefficient rounding splits an m-fold root by O(eps**(1/m)). Check
    # that conditioning bound as well as the much stricter backward error.
    multiplicity = max(np.count_nonzero(nodes == node) for node in nodes)
    assert error <= 8*np.finfo(float).eps**(1/multiplicity)*max(1, max(abs(nodes)))
    np.testing.assert_allclose(np.polynomial.polynomial.polyfromroots(output), p, rtol=2e-12, atol=2e-14)


@pytest.mark.parametrize("scale", [1e-60, 1, 1e60])
def test_complex_quartic_scaling(complex_engine, scale):
    from itertools import permutations
    nodes = scale*np.array([-.2+.1j, -.1-.05j, .1+.03j, .2+.07j])
    p = np.ascontiguousarray(np.polynomial.polynomial.polyfromroots(nodes))
    output = np.empty(4, dtype=np.complex128)
    complex_engine.roots(4, p.view(float), output.view(float))
    error = min(np.max(np.abs(output/scale - candidate)) for candidate in permutations(nodes/scale))
    assert error < 2e-12
