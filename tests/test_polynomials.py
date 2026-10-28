"""Exercise the public C polynomial helpers without adding Python bindings."""
import ctypes
import itertools
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest


@pytest.fixture(scope="module")
def polynomials(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    directory = tmp_path_factory.mktemp("polynomials")
    source = directory/"polynomials.c"
    lines = []
    for n in range(1, 5):
        lines.append(f'#include "src/algpy/polynomials/poly{n}.h"')
        lines.append(f'''
void roots{n}(const double *p, double *out) {{
    complex z[{n}];
    poly{n}_coeffs_to_roots(p, z);
    for (int i = 0; i < {n}; i++) {{
        out[2*i] = z[i].data[0]; out[2*i+1] = z[i].data[1];
    }}
}}
void coeffs{n}(const double *input, double *p) {{
    complex z[{n}];
    for (int i = 0; i < {n}; i++) {{
        z[i] = complex_make(input[2*i], input[2*i+1]);
    }}
    poly{n}_roots_to_coeffs(z, p);
}}
''')
    source.write_text("\n".join(lines))
    library = directory/"polynomials.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Werror",
                    "-shared", "-fPIC", "-I", str(root), str(source), "-lm",
                    "-o", str(library)], check=True, capture_output=True, text=True)
    lib = ctypes.CDLL(str(library))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    for n in range(1, 5):
        for name in ("roots", "coeffs"):
            fn = getattr(lib, f"{name}{n}")
            fn.argtypes = [pointer, pointer]
            fn.restype = None

    class Polynomials:
        def roots(self, p):
            p = np.ascontiguousarray(p, dtype=float)
            z = np.empty((len(p)-1, 2))
            getattr(lib, f"roots{len(p)-1}")(p, z)
            return z[:, 0] + 1j*z[:, 1]

        def coeffs(self, z):
            z = np.ascontiguousarray(z, dtype=np.complex128)
            p = np.empty(len(z)+1)
            getattr(lib, f"coeffs{len(z)}")(z.view(float), p)
            return p

    return Polynomials()


CASES = [
    [-3.0], [0.0], [2.5],
    [1.0, 2.0], [-2.0, 3.0], [0.0, 0.0], [2.0, 2.0], [1+2j, 1-2j],
    [-3.0, 1.0, 2.0], [-1.0, 0.0, 1.0], [2.0, 2.0, 2.0],
    [-2.0, 1.0, 1.0], [3.0, 1j, -1j], [0.0, 2+3j, 2-3j],
    [-3.0, -1.0, 2.0, 4.0], [0.0, 0.0, 0.0, 0.0], [2.0]*4,
    [-3.0, -3.0, 2.0, 2.0], [1.0, 1.0, 1.0, -2.0],
    [-3.0, 2.0, 1+2j, 1-2j], [1j, -1j, 2j, -2j],
    [1+2j, 1-2j, -3+4j, -3-4j], [1+2j, 1-2j]*2,
    # Shifted biquadratic with four non-real roots.
    [2+1j, 2-1j, 2+3j, 2-3j],
]


def assert_roots(actual, expected, rtol=2e-12, atol=2e-12):
    assert np.isfinite(actual).all(), actual
    # Order is deliberately unspecified; preserve multiplicities when matching.
    assert any(np.all(np.abs(actual-np.array(order)) <= atol+rtol*np.abs(order))
               for order in itertools.permutations(expected)), (actual, expected)


@pytest.mark.parametrize("roots", CASES)
def test_vieta_and_known_roots(polynomials, roots):
    coefficients = np.polynomial.polynomial.polyfromroots(roots).real
    for order in itertools.permutations(roots):
        np.testing.assert_allclose(polynomials.coeffs(order), coefficients,
                                   rtol=2e-15, atol=2e-15)
    # Repeated roots are intrinsically ill-conditioned when recovering from coefficients.
    repeated = len(set(roots)) != len(roots)
    assert_roots(polynomials.roots(coefficients), roots,
                 rtol=2e-5 if repeated else 2e-12, atol=2e-5 if repeated else 2e-12)


@pytest.mark.parametrize("n", range(1, 5))
@pytest.mark.parametrize("scale", [1.0, -7.0])
def test_nonmonic_and_coefficient_scaling(polynomials, n, scale):
    expected = np.arange(1, n+1, dtype=float)
    p = np.polynomial.polynomial.polyfromroots(expected)*scale
    assert_roots(polynomials.roots(p), expected)


def test_quadratic_small_root(polynomials):
    assert_roots(polynomials.roots([1.0, -1e6, 1.0]),
                 [1e-6, 1e6-1e-6], rtol=2e-12, atol=0)


@pytest.mark.parametrize("n", range(1, 5))
def test_random_coefficients(polynomials, n):
    rng = np.random.default_rng(400+n)
    for _ in range(100):
        p = np.r_[rng.uniform(-3, 3, n), rng.uniform(0.5, 3)]
        roots = polynomials.roots(p)
        assert_roots(roots, np.polynomial.polynomial.polyroots(p), rtol=2e-10, atol=2e-10)
        np.testing.assert_allclose(polynomials.coeffs(roots), p/p[-1], rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("parameters", [
    (-1, -2, -3, -4), (1, 1, 1, 1), (-1, 0, 1, 0), (0, 0, 0, 0),
])
def test_clifford4_eigenvalues(parameters):
    from algpy.algebras.general import clifford4

    cls = clifford4(*parameters)
    # Scalar and generator inputs give repeated eigenvalues.
    assert_roots(cls.one().eivals(), [1.0]*4)
    generator_root = np.sqrt(complex(parameters[0]))
    assert_roots(cls.basis()[1].eivals(), [generator_root, -generator_root]*2)
    rng = np.random.default_rng(44)
    for _ in range(25):
        value = cls(*rng.normal(size=16))
        np.testing.assert_allclose(
            np.polynomial.polynomial.polyfromroots(value.eivals()),
            value.minpoly(), rtol=2e-10, atol=2e-10,
        )
