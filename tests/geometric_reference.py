"""Independent Clifford products and checks for fixed geometric algebras."""
import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algebra_reference import FIXED_METRICS
from field_reference import reference


def pga_matrix_basis(dimension):
    """Faithful complex blocks, independent of the Clifford product formulas."""
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    generators = [np.kron(x, np.eye(2)), np.kron(z, x),
                  np.kron(z, y), np.kron(z, z)]
    parity = generators[0]@generators[1]@generators[2]@generators[3]
    if dimension == 5:
        generators.append(parity)
        generators = [np.block([[g, np.zeros_like(g)], [np.zeros_like(g), -g]])
                      for g in generators]
        parity = np.block([[np.zeros((4, 4)), np.eye(4)],
                           [np.eye(4), np.zeros((4, 4))]])
    size = len(parity)
    zero = np.zeros_like(parity)
    null = np.block([[zero, parity], [zero, zero]])
    generators = [null, *[np.block([[g, zero], [zero, g]]) for g in generators]]
    basis = []
    for mask in range(2**len(generators)):
        value = np.eye(2*size, dtype=complex)
        for bit, generator in enumerate(generators):
            if mask & (1 << bit):
                value = value@generator
        basis.append(value)
    return np.array(basis)


def multiply(parameters, a, b):
    result = np.zeros(len(a))
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            sign = (-1)**sum(bool(i & (1 << p)) and bool(j & (1 << q))
                              for p in range(len(parameters)) for q in range(p))
            for p, square in enumerate(parameters):
                if i & j & (1 << p):
                    sign *= square
            result[i ^ j] += sign*x*y
    return result


def check_arithmetic(algebra):
    parameters = FIXED_METRICS[algebra.__name__]
    rng = np.random.default_rng(142)
    for _ in range(12):
        x, y, z = rng.uniform(-.2, .2, (3, algebra.N))
        a, b, c = algebra(*x), algebra(*y), algebra(*z)
        np.testing.assert_allclose(algebra.to_float(a*b), multiply(parameters, x, y), atol=2e-15)
        np.testing.assert_allclose(algebra.to_float(np.square(a)), multiply(parameters, x, x), atol=2e-15)
        np.testing.assert_allclose(algebra.to_float((a*b)*c), algebra.to_float(a*(b*c)), atol=2e-15)
        inv = (a + 2).inv()
        np.testing.assert_allclose(algebra.to_float((a + 2)*inv), algebra.to_float(algebra.one()), atol=2e-14)


def check_nilpotent(algebra, name, function):
    parameters = FIXED_METRICS[algebra.__name__]
    nilpotent = np.zeros(algebra.N)
    if algebra.__name__ == "pga1":
        nilpotent[1] = .25
        order = 2
    elif algebra.__name__.startswith("pga"):
        nilpotent[[1, 2, 6]] = .25
        order = 4
    elif algebra.__name__ == "cga1":
        nilpotent[[1, 2]] = .25
        order = 2
    else:
        # e0+e3 is null. With e1,e2 it spans a PGA2 subalgebra.
        nilpotent[[1, 8, 2, 6]] = .25
        order = 4
    w = 2 if name == "acosh" else 1.5 if name in {"inv", "log", "sqrt"} else .3
    expected = np.zeros(algebra.N)
    expected[0] = function(w)
    power = nilpotent.copy()
    factorial = 1
    for derivative in range(1, order):
        assert power.any()
        factorial *= derivative
        expected += reference(function, complex(w), derivative).real/factorial*power
        power = multiply(parameters, power, nilpotent)
    np.testing.assert_array_equal(power, 0)
    value = nilpotent.copy()
    value[0] = w
    actual = function(algebra(*value))
    np.testing.assert_allclose(algebra.to_float(actual), expected, rtol=3e-12, atol=3e-13)


def check_standalone_header(algebra, tmp_path, degree=3):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    name = algebra.__name__
    root = Path(__file__).resolve().parents[1]
    code = f'''
#include "src/algpy/algebras/clifford/{name}.h"
#include <assert.h>
#if defined(ALGPY_CLIFFORD1_H) || defined(ALGPY_CLIFFORD2_H) || \\
    defined(ALGPY_CLIFFORD3_H) || defined(ALGPY_CLIFFORD4_H) || \\
    defined(ALGPY_CLIFFORD5_H) || \\
    defined(ALGPY_CAYLEYDICKSON1_H) || defined(ALGPY_CAYLEYDICKSON2_H) || \\
    defined(ALGPY_MULTICOMPLEX1_H)
#error Fixed specializations must not include general algebra headers
#endif
'''
    real_order = 7 if name == "pga4" else 3
    complex_order = 3 if name in ("cga3", "pga4") else 1
    for order in range(real_order+1):
        expression = (f"{math.factorial(degree)//math.factorial(degree-order)}.*double_ipow(x, {degree-order})"
                      if order <= degree else "0")
        if order == 0:
            expression += " + 2*x + 4"
        elif order == 1:
            expression += " + 2"
        code += f"static double d{order}(double x) {{ (void)x; return {expression}; }}\n"
    for order in range(complex_order+1):
        expression = (f"complex_lmul(complex_ipow(x, {degree-order}), {math.factorial(degree)//math.factorial(degree-order)}.)"
                      if order <= degree else "complex_zero()")
        if order == 0:
            expression = f"complex_ladd(complex_add({expression}, complex_lmul(x, 2)), 4)"
        elif order == 1:
            expression = f"complex_ladd({expression}, 2)"
        code += f"static complex c{order}(complex x) {{ (void)x; return {expression}; }}\n"
    callbacks = ", ".join([*(f"d{i}" for i in range(real_order+1)),
                           *(f"c{i}" for i in range(complex_order+1))])
    code += f'''
static void check_polynomial({name} a) {{
    {name} actual = {name}_apply(a, {callbacks});
    {name} expected = {name}_add({name}_ipow(a, {degree}), {name}_ladd({name}_lmul(a, 2), 4));
    for (int i = 0; i < {algebra.N}; i++) assert(fabs(actual.data[i] - expected.data[i]) < 2e-12);
}}
int main(void) {{
    {name} a = {name}_from_re(.3);
    for (int i = 1; i < {algebra.N}; i++) a.data[i] = (i%3 - 1)*.04;
    check_polynomial(a);
    {name} near = {name}_from_re(.3);
    near.data[1] = .25;
    near.data[2] = 1e-10;
    check_polynomial(near);
'''
    if algebra.N >= 8:
        code += f"    near.data[6] = .5e-10;\n    check_polynomial(near);\n"
    code += f"    {name} sum = {name}_zero();\n"
    for grade in range(algebra.N.bit_length()):
        code += f"    sum = {name}_add(sum, {name}_grade{grade}(a));\n"
    code += f"    assert({name}_equal(sum, a));\n    return 0;\n}}\n"
    source = tmp_path/"standalone.c"
    source.write_text(code)
    executable = tmp_path/"standalone"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Wextra", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(executable)], check=True, capture_output=True, text=True)
    subprocess.run([str(executable)], check=True)


def check_projector_sqrt(algebra):
    generator = algebra.gens()[1]
    projector = (algebra.one() + generator)/2
    np.testing.assert_allclose(algebra.to_float(np.sqrt(projector)), algebra.to_float(projector),
                               rtol=3e-13, atol=3e-14)


def check_pga_radical(algebra):
    row = np.random.default_rng(143).uniform(-.04, .04, algebra.N)
    row[0] = 1.5
    row[1::2] = 0
    base = algebra(*row)
    row[1::2] = np.linspace(-.2, .3, algebra.N//2)
    value = algebra(*row)
    np.testing.assert_allclose(value.minpoly(), base.minpoly(), rtol=2e-14, atol=2e-14)
    assert value.norm() == pytest.approx(base.norm())
    for function in [np.exp, np.log, np.sin]:
        np.testing.assert_allclose(algebra.to_float(function(value))[::2],
                                   algebra.to_float(function(base))[::2], rtol=3e-12, atol=3e-13)
    translation = .25*algebra.gens()[0]*algebra.gens()[1]
    np.testing.assert_allclose(algebra.to_float(np.exp(translation)),
                               algebra.to_float(algebra.one() + translation), atol=2e-14)


def check_pga_projected_nilpotent_sqrt(algebra):
    generators = algebra.gens()
    projector = (algebra.one() + generators[1])/2
    center = generators[0]*generators[1]*generators[2]
    nilpotent = .25*center*projector
    with np.errstate(all="raise"):
        actual = np.sqrt(projector + nilpotent)
    expected = projector + nilpotent/2
    np.testing.assert_allclose(algebra.to_float(actual), algebra.to_float(expected), atol=2e-14)


def check_cga_points(algebra):
    generators = algebra.gens()
    origin = (generators[0] - generators[1])/2
    infinity = generators[0] + generators[1]
    points = []
    coordinates = [[.25]*len(generators[2:]), [-.5]*len(generators[2:])]
    for row in coordinates:
        value = origin + .5*np.dot(row, row)*infinity
        for coefficient, generator in zip(row, generators[2:]):
            value = value + coefficient*generator
        np.testing.assert_array_equal(algebra.to_float(value*value), np.zeros(algebra.N))
        points.append(value)
    inner = (points[0]*points[1] + points[1]*points[0])/2
    expected = -.5*np.sum((np.array(coordinates[0]) - coordinates[1])**2)
    np.testing.assert_allclose(algebra.to_float(inner), algebra.to_float(algebra(expected)), atol=2e-14)
