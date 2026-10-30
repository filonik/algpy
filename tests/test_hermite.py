"""The polynomial engine has no algebra dependency, including at Jordan blocks."""
import ctypes
import itertools
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from field_reference import FUNCTIONS
from field_reference import center, operation


@pytest.fixture(scope="module")
def interpolate(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    directory = tmp_path_factory.mktemp("hermite")
    cases = []
    for index, (name, _) in enumerate(FUNCTIONS):
        callbacks = (f"double_{name}, double_{name}_deriv1, double_{name}_deriv2, "
                     f"double_{name}_deriv3, complex_{name}, complex_{name}_deriv1, result")
        cases.append(f"case {index}: return solve ? poly4_apply(p, w, {callbacks}) : "
                     f"hermite_interpolate(n, w, nodes, {callbacks});")
    source = directory/"hermite.c"
    source.write_text('''
#include "src/algpy/polynomials/hermite.h"
int evaluate(int function, int solve, int n, double w, const double *input,
             const double *p, double *result) {
    complex nodes[4];
    for (int i = 0; i < n; i++) nodes[i] = complex_make(input[2*i], input[2*i+1]);
    switch (function) {
'''+"\n".join(cases)+'''
    default: return 0;
    }
}
''')
    library_path = directory/"hermite.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-shared", "-fPIC", "-I", str(root), str(source), "-lm",
                    "-o", str(library_path)], check=True, capture_output=True, text=True)
    library = ctypes.CDLL(str(library_path))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    library.evaluate.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double,
                                 pointer, pointer, pointer]
    library.evaluate.restype = ctypes.c_int

    def evaluate(name, nodes, solve=False):
        nodes = np.ascontiguousarray(nodes, dtype=np.complex128)
        p = np.ascontiguousarray(np.polynomial.polynomial.polyfromroots(nodes).real)
        output = np.empty(4)
        success = library.evaluate([name for name, _ in FUNCTIONS].index(name), solve,
            len(nodes), center(name), nodes.view(float), p, output)
        assert success, (name, nodes, solve)
        return output

    return evaluate


def companion(nodes):
    p = np.polynomial.polynomial.polyfromroots(nodes).real
    result = np.diag(np.ones(len(nodes)-1), -1)
    result[:, -1] = -p[:-1]
    return result


def check_matrix(interpolate, name, nodes, solve):
    pytest.importorskip("scipy.linalg")
    from scipy_reference import FUNCTIONS as references
    c = interpolate(name, nodes, solve)
    matrix = companion(nodes)
    identity = np.eye(len(nodes))
    actual = c[0]*identity + c[1]*matrix + c[2]*(matrix@matrix) + c[3]*(matrix@matrix@matrix)
    expected = dict(references)[operation(name)](center(name)*identity + matrix)
    # Companion matrices retain a substantial nilpotent part even when the
    # nodes coalesce. Errors in the cubic coefficient cannot hide in tiny b^3.
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=3e-12)


NODES = [
    [-.2, -.1, .1, .2],
    [-.2, .2, .1+.2j, .1-.2j],
    [-.1+.2j, -.1-.2j, .1+.3j, .1-.3j],
    [0]*4, [.1]*4, [-.1, -.1, .2, .2], [.1, .1, .1, -.3],
    [.2j, .2j, -.2j, -.2j],
    [-.2, -.2, .2], [0, 0, 0], [.2j, -.2j], [0],
]


@pytest.mark.parametrize("nodes", NODES, ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_hermite_and_quartic(interpolate, nodes, name, function):
    check_matrix(interpolate, name, nodes, False)
    if len(nodes) == 4:
        check_matrix(interpolate, name, nodes, True)


@pytest.mark.parametrize("delta", [1e-2, 1e-3, 1e-6, 1e-12, 1e-40, 0])
@pytest.mark.parametrize("nodes", [[-.3, -.1, .2, .4], [.5, -.3, -.1+.2j, -.1-.2j],
                                  [1j, 1j, -1j, -1j]], ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_nearly_confluent(interpolate, delta, nodes, name, function):
    nodes = np.array(nodes)*delta
    check_matrix(interpolate, name, nodes, False)
    check_matrix(interpolate, name, nodes, True)


@pytest.mark.parametrize("nodes", [[1e-10+.02j, 1e-10-.02j, -1e-10+.02j, -1e-10-.02j],
                                  [-.1+1e-10j, -.1-1e-10j, .2, .3]], ids=str)
@pytest.mark.parametrize("name", ["exp", "log", "sqrt"])
def test_nearby_node_ordering(interpolate, nodes, name):
    for permutation in itertools.permutations(nodes):
        check_matrix(interpolate, name, permutation, False)
