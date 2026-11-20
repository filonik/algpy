"""Real degree-eight Hermite interpolation for degenerate Clifford5 centres."""
import ctypes
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest
from field_reference import FUNCTIONS, center, operation


@pytest.fixture(scope="module")
def interpolate8(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    directory = tmp_path_factory.mktemp("hermite8")
    cases = []
    for index, (name, _) in enumerate(FUNCTIONS):
        real = ', '.join(f'double_{name}_deriv{i}' for i in range(1, 8))
        analytic = ', '.join(f'complex_{name}_deriv{i}' for i in range(1, 4))
        cases.append(f'''case {index}: {{
        double (*derivatives[7])(double) = {{{real}}};
        complex (*cderivatives[3])(complex) = {{{analytic}}};
        return hermite8_interpolate(n, w, nodes, double_{name}, derivatives,
                                    complex_{name}, cderivatives, output);
    }}''')
    source = directory/"interpolate.c"
    source.write_text('''
#include "src/algpy/polynomials/hermite8.h"
int evaluate(int function, int n, double w, const double *input, double *output) {
    complex nodes[8];
    for (int i = 0; i < n; i++) nodes[i] = complex_make(input + 2*i);
    switch (function) {
'''+"\n".join(cases)+'''
    default: return 0;
    }
}
''')
    library_path = directory/"interpolate.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-shared", "-fPIC", "-I", str(root), str(source), "-lm", "-o", str(library_path)],
                   check=True, capture_output=True, text=True)
    library = ctypes.CDLL(str(library_path))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    library.evaluate.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double, pointer, pointer]
    library.evaluate.restype = ctypes.c_int

    def evaluate(name, nodes, value=None):
        nodes = np.ascontiguousarray(nodes, dtype=np.complex128)
        output = np.empty(8)
        assert library.evaluate([name for name, _ in FUNCTIONS].index(name), len(nodes),
                                 center(name) if value is None else value,
                                 nodes.view(float), output), (name, nodes)
        return output

    return evaluate


NODES = [[0]*8, [.1]*8, [-.1]*4 + [.2]*4, [.2j]*4 + [-.2j]*4,
         np.repeat([-.2, -.1, .1, .2], 2),
         np.repeat([-.1+.2j, -.1-.2j, .2+.1j, .2-.1j], 2)]


def check_matrix(interpolate8, name, nodes):
    from scipy_reference import FUNCTIONS as references
    p = np.polynomial.polynomial.polyfromroots(nodes).real
    matrix = np.diag(np.ones(len(nodes)-1), -1)
    matrix[:, -1] = -p[:-1]
    coefficients = interpolate8(name, nodes)
    identity = np.eye(len(nodes))
    power, actual = identity, np.zeros_like(matrix)
    for coefficient in coefficients:
        actual += coefficient*power
        power = power@matrix
    expected = dict(references)[operation(name)](center(name)*identity + matrix)
    # The seventh Newton difference of separated nodes is more sensitive
    # to value roundoff than the cubic interpolants used by smaller algebras.
    np.testing.assert_allclose(actual, expected, rtol=3e-9, atol=5e-10)


@pytest.mark.parametrize("nodes", NODES, ids=str)
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_real_and_complex_multiplicities(interpolate8, nodes, name, function):
    check_matrix(interpolate8, name, nodes)


@pytest.mark.parametrize("scale", [.2, .02, .002, 1e-6, 1e-12, 0])
@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
def test_eight_nearly_coincident_nodes(interpolate8, scale, name, function):
    nodes = np.repeat(np.array([-.7, -.2, .3, .8])*scale, 2)
    check_matrix(interpolate8, name, nodes)


def test_saturated_value_retains_small_derivatives(interpolate8):
    expm = pytest.importorskip("scipy.linalg").expm
    nodes = np.repeat([-.001, -.0003, .0002, .0008], 2)
    p = np.polynomial.polynomial.polyfromroots(nodes)
    matrix = np.diag(np.ones(7), -1)
    matrix[:, -1] = -p[:-1]
    coefficients = interpolate8("tanh", nodes, 100)
    power, actual = np.eye(8), np.zeros((8, 8))
    for coefficient in coefficients:
        actual += coefficient*power
        power = power@matrix
    # tanh(100+x) = 1 - 2*exp(-200-2*x) to relative O(exp(-200)).
    # Off-diagonal entries retain the tiny change after the value rounds to 1.
    expected = -2*np.exp(-200)*expm(-2*matrix)
    off_diagonal = ~np.eye(8, dtype=bool)
    np.testing.assert_allclose(actual[off_diagonal], expected[off_diagonal], rtol=3e-10, atol=0)


@pytest.mark.parametrize("name,function", [("cos", np.cos), ("cosh", np.cosh)])
def test_even_functions_at_zero_with_conjugate_nodes(interpolate8, name, function):
    # The odd coefficients nearly vanish. Their imaginary roundoff must be
    # measured before cancellation in the complete homogeneous moments.
    upper = np.array([-.01+.02j, .03+.04j, .02+.05j, -.04+.06j])
    nodes = np.r_[upper, upper.conj()]
    coefficients = interpolate8(name, nodes, 0)
    polynomial = np.polynomial.polynomial.polyfromroots(nodes).real
    matrix = np.diag(np.ones(7), -1)
    matrix[:, -1] = -polynomial[:-1]
    power, actual = np.eye(8), np.zeros((8, 8))
    for coefficient in coefficients:
        actual += coefficient*power
        power = power@matrix
    from scipy_reference import FUNCTIONS as references
    expected = dict(references)[name](matrix)
    np.testing.assert_allclose(actual, expected, rtol=3e-10, atol=2e-12)
