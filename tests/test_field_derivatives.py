"""Compare the explicit C derivatives with Cauchy's integral formula."""
import ctypes
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest


FUNCTIONS = [
    ("exp", np.exp), ("expm1", np.expm1), ("log", np.log), ("log1p", np.log1p),
    ("inv", np.reciprocal), ("square", np.square), ("sqrt", np.sqrt),
    ("sin", np.sin), ("cos", np.cos), ("tan", np.tan),
    ("asin", np.arcsin), ("acos", np.arccos), ("atan", np.arctan),
    ("sinh", np.sinh), ("cosh", np.cosh), ("tanh", np.tanh),
    ("asinh", np.arcsinh), ("acosh", np.arccosh), ("atanh", np.arctanh),
]
INDICES = {name: i for i, (name, _) in enumerate(FUNCTIONS)}


@pytest.fixture(scope="module")
def derivatives(tmp_path_factory):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    directory = tmp_path_factory.mktemp("field_derivatives")
    root = Path(__file__).resolve().parents[1]
    tables = []
    for field in ["double", "complex"]:
        rows = ["    {"+", ".join(f"{field}_{name}_deriv{n}" for n in range(1, 5))+"}"
                for name, _ in FUNCTIONS]
        tables.append(f"static {field} (*{field}_derivatives[][4])({field}) = {{\n"
                      +",\n".join(rows)+"\n};")
    source = directory/"derivatives.c"
    source.write_text('''
#include "src/algpy/fields/double.h"
#include "src/algpy/fields/complex.h"
'''+"\n".join(tables)+'''
double evaluate_double(int function, int order, double x)
{
    return double_derivatives[function][order-1](x);
}
void evaluate_complex(int function, int order, double x, double y, double out[2])
{
    complex result = complex_derivatives[function][order-1](complex_make(x, y));
    out[0] = result.data[0];
    out[1] = result.data[1];
}
''')
    library_path = directory/"derivatives.so"
    subprocess.run([compiler, "-std=c99", "-O2", "-Wall", "-Wextra", "-Werror",
                    "-shared", "-fPIC", "-I", str(root), str(source), "-lm",
                    "-o", str(library_path)], check=True, capture_output=True, text=True)
    library = ctypes.CDLL(str(library_path))
    library.evaluate_double.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double]
    library.evaluate_double.restype = ctypes.c_double
    library.evaluate_complex.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_double,
                                        ctypes.c_double, ctypes.POINTER(ctypes.c_double)]
    library.evaluate_complex.restype = None

    def evaluate(field, name, order, value):
        if field == "double":
            return library.evaluate_double(INDICES[name], order, value)
        result = (ctypes.c_double*2)()
        library.evaluate_complex(INDICES[name], order, value.real, value.imag, result)
        return complex(result[0], result[1])

    return evaluate


def reference(function, value, order):
    # Extract a Taylor coefficient with a discrete contour integral. This
    # checks the derivatives using only the function, without their formulas
    # or the cancellation of a fourth-order finite difference.
    radius = 0.125
    circle = np.exp(2j*np.pi*np.arange(256)/256)
    return (math.factorial(order)/radius**order
            * np.mean(function(value+radius*circle)*circle**(-order)))


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_real_derivatives(derivatives, name, function, order):
    points = [-2, -0.25, 0, 0.4, 2]
    if name in ["log", "sqrt"]:
        points = [0.25, 0.75, 1.5, 3]
    elif name == "log1p":
        points = [-0.75, -0.25, 0, 0.4, 2]
    elif name == "inv":
        points = [-2, -0.5, 0.5, 2]
    elif name in ["asin", "acos", "atanh"]:
        points = [-0.75, -0.25, 0, 0.4, 0.75]
    elif name == "tan":
        points = [-1, -0.25, 0, 0.4, 1]
    elif name == "acosh":
        points = [1.25, 1.5, 2, 3]
    for value in points:
        actual = derivatives("double", name, order, value)
        expected = reference(function, complex(value), order).real
        assert actual == pytest.approx(expected, rel=3e-10, abs=3e-10)


@pytest.mark.parametrize("name,function", FUNCTIONS, ids=[name for name, _ in FUNCTIONS])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_complex_derivatives(derivatives, name, function, order):
    points = [0.3+0.4j, -0.5+0.6j, 1.6+1.2j, -2+0.7j, 0.5j]
    points += [z.conjugate() for z in points]
    for value in points:
        actual = derivatives("complex", name, order, value)
        expected = reference(function, value, order)
        assert actual == pytest.approx(expected, rel=3e-10, abs=3e-10)


@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_complex_branch_cuts(derivatives, order):
    # Analytic one-sided limits at points on the cuts, including signed zero
    # on either axis. These detect square-root branch errors independently
    # of the contour tests, whose circles must avoid the cuts.
    sqrt = [-1j/(2*np.sqrt(2)), -1j/(8*np.sqrt(2)),
            -3j/(32*np.sqrt(2)), -15j/(128*np.sqrt(2))][order-1]
    asin = [1j/np.sqrt(3), -2j/(3*np.sqrt(3)),
            1j/np.sqrt(3), -22j/(9*np.sqrt(3))][order-1]
    asinh = [-1j/np.sqrt(3), 2/(3*np.sqrt(3)),
             1j/np.sqrt(3), -22/(9*np.sqrt(3))][order-1]
    acosh = [-1/np.sqrt(3), -2/(3*np.sqrt(3)),
             -1/np.sqrt(3), -22/(9*np.sqrt(3))][order-1]
    for zero in [0.0, -0.0]:
        upper = not np.signbit(zero)
        for name, value, expected in [
            ("sqrt", complex(-2, zero), sqrt if upper else sqrt.conjugate()),
            ("asin", complex(2, zero), asin if upper else asin.conjugate()),
            ("acos", complex(2, zero), -asin if upper else -asin.conjugate()),
            ("asinh", complex(zero, 2), asinh if upper else -asinh),
            ("asinh", complex(zero, -2), np.conj(asinh if upper else -asinh)),
            ("acosh", complex(-2, zero), acosh),
        ]:
            actual = derivatives("complex", name, order, value)
            assert actual == pytest.approx(expected, rel=2e-14, abs=2e-14)


@pytest.mark.parametrize("name", ["log", "inv", "sqrt", "atan", "asinh", "acosh"])
def test_large_real_inputs(derivatives, name):
    x = 1e60
    # Leading terms suffice at this magnitude, with relative corrections at
    # most 1e-120. Derivatives up to order four are still representable here.
    expected = {
        "log": [1/x, -1/x**2, 2/x**3, -6/x**4],
        "inv": [-1/x**2, 2/x**3, -6/x**4, 24/x**5],
        "sqrt": [c*x**(0.5-n) for n, c in enumerate([0.5, -0.25, 0.375, -0.9375], 1)],
        "atan": [1/x**2, -2/x**3, 6/x**4, -24/x**5],
        "asinh": [1/x, -1/x**2, 2/x**3, -6/x**4],
        "acosh": [1/x, -1/x**2, 2/x**3, -6/x**4],
    }[name]
    for order, value in enumerate(expected, 1):
        assert derivatives("double", name, order, x) == pytest.approx(value, rel=5e-15, abs=0)


def test_tanh_derivatives_after_saturation(derivatives):
    # tanh rounds to one, but its derivatives are nonzero at this input.
    for order, factor in enumerate([4, -8, 16, -32], 1):
        assert derivatives("double", "tanh", order, 100) == pytest.approx(
            factor*np.exp(-200), rel=5e-15, abs=0)


@pytest.mark.parametrize("name,function", FUNCTIONS[:4]+FUNCTIONS[6:],
                         ids=[name for name, _ in FUNCTIONS[:4]+FUNCTIONS[6:]])
def test_dual_callers(name, function):
    from algpy.algebras import dual1, dual2

    x = 1.5 if name in ["log", "sqrt", "acosh"] else 0.4
    first = reference(function, complex(x), 1).real
    second = reference(function, complex(x), 2).real
    for cls, coefficients, expected in [
        (dual1, [x, 0.3], [function(x), 0.3*first]),
        (dual2, [x, 0.3, -0.2, 0.1],
         [function(x), 0.3*first, -0.2*first, 0.1*first-0.06*second]),
    ]:
        actual = getattr(cls(*coefficients), name)()
        np.testing.assert_allclose(cls.to_float(actual), expected, rtol=3e-10, atol=3e-10)
