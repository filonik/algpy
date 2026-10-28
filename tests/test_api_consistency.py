"""The same C entry points and NumPy operations exist for every algebra."""
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from algpy import algebras
from algpy.algebras import general


FAMILIES = ("cayleydickson", "clifford", "multicomplex")
CONCRETE = (
    ("real", 1), ("complex1", 2), ("quaternion", 4),
    ("octonion", 8), ("sedenion", 16),
    *((f"{family}{n}", 2**n)
      for family in ("complex", "dual", "splitcomplex") for n in range(1, 5)
      if (family, n) != ("complex", 1)),
)
UNARY = (
    ("positive", None), ("negative", None), ("conjugate", "conjugate"),
    ("reciprocal", "inv"), ("square", "square"),
    ("exp", "exp"), ("expm1", "expm1"), ("log", "log"),
    ("log1p", "log1p"), ("sqrt", "sqrt"),
    ("sin", "sin"), ("cos", "cos"), ("tan", "tan"),
    ("arcsin", "asin"), ("arccos", "acos"), ("arctan", "atan"),
    ("sinh", "sinh"), ("cosh", "cosh"), ("tanh", "tanh"),
    ("arcsinh", "asinh"), ("arccosh", "acosh"), ("arctanh", "atanh"),
)


@pytest.fixture(scope="module", params=[
    *((family, n, True) for family in FAMILIES for n in range(1, 5)),
    *((name, n, False) for name, n in CONCRETE),
], ids=lambda p: f"{p[0]}{p[1] if p[2] else ''}")
def algebra(request):
    name, n, parameterized = request.param
    if parameterized:
        return getattr(general, f"{name}{n}")(*([-1.0]*n))
    return getattr(algebras, name)


def assert_components(cls, a, b):
    np.testing.assert_allclose(cls.to_float(a), cls.to_float(b), atol=1e-12, equal_nan=True)


@pytest.mark.parametrize("name,method", UNARY)
def test_unary_scalar_array_and_method(algebra, name, method):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    # Include a strided input and an explicitly supplied output buffer.
    values = np.array([x, x, x], dtype=cls)[::2]
    out = np.empty_like(values)
    with np.errstate(all="ignore"):
        result = getattr(np, name)(values, out=out)
        scalar = getattr(np, name)(x)
        assert result is out and result.dtype == cls.dtype
        assert_components(cls, result[0], scalar)
        assert_components(cls, result[1], scalar)
        if method:
            assert_components(cls, getattr(x, method)(), scalar)


@pytest.mark.parametrize("name", ("add", "subtract", "multiply", "divide", "power"))
def test_binary_cast_broadcast_and_operator(algebra, name):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    values = np.array([[x], [x]], dtype=cls)
    with np.errstate(all="ignore"):
        result = getattr(np, name)(values, np.array([1.0, 2.0]))
        assert result.shape == (2, 2) and result.dtype == cls.dtype
        scalar = {"add": lambda: x+2, "subtract": lambda: x-2,
                  "multiply": lambda: x*2, "divide": lambda: x/2,
                  "power": lambda: x**2}[name]()
        assert_components(cls, result[1, 1], scalar)
        if name == "power":
            assert_components(cls, 2**x, np.power(cls(2), x))


def test_predicates_parts_and_inverse(algebra):
    cls = algebra
    x = cls(*np.linspace(0.1, 0.3, cls.N))
    values = np.array([x, cls.zero()], dtype=cls)
    assert np.equal(values, x).tolist() == [True, False]
    assert np.not_equal(values, x).tolist() == [False, True]
    assert_components(cls, np.real(x) + np.imag(x), x)
    assert_components(cls, np.real(x), x.re())
    assert_components(cls, np.imag(x), x.im())
    assert x.to_re() == pytest.approx(0.1)
    assert x.is_real() == (cls.N == 1)
    assert cls.from_re(2).is_real()
    assert_components(cls, cls.from_re(2), cls(2))
    # Shift away from zero divisors in the split-complex specializations.
    invertible = x + 4
    assert_components(cls, invertible.inv()*invertible, cls.one())
    with np.errstate(all="ignore"):
        assert_components(cls, invertible.mix(2, 0.5),
                          invertible * np.power(invertible.inv()*2, 0.5))


def test_norm_static_return_type(algebra):
    cls = algebra
    norm_type = complex if hasattr(cls, "K") or cls.__name__.startswith("splitcomplex") else float
    for value, expected in ((cls.zero(), 0), (cls.one(), 1), (cls(-2), 2)):
        norm = value.norm()
        assert type(norm) is norm_type
        assert norm == pytest.approx(expected)
    value = cls(2, *(np.arange(1, cls.N)/100))
    norm = value.norm()
    assert type(norm) is norm_type
    if cls.N == 1:
        assert norm == 2
    else:
        polynomial = value.minpoly()
        expected = complex(polynomial[0])**(1/(len(polynomial)-1))
        assert norm == pytest.approx(expected)


@pytest.mark.parametrize("family", FAMILIES)
@pytest.mark.parametrize("n", range(1, 5))
def test_general_norm_type_does_not_follow_parameters(family, n):
    cls = getattr(general, f"{family}{n}")(*([-1]*n))
    generator = cls.basis()[1]
    assert type(generator.norm()) is complex
    assert generator.norm() == pytest.approx(1)
    cls.K = [1]*n
    norm = generator.norm()
    assert type(norm) is complex
    # Cayley-Dickson has a quadratic norm; the other families' higher
    # even-degree norms can remain real for this particular generator.
    if family == "cayleydickson" or n == 1:
        assert norm == pytest.approx(1j)
    assert type(cls.one().norm()) is complex


@pytest.mark.parametrize("n", range(1, 5))
def test_splitcomplex_negative_norm(n):
    cls = getattr(algebras, f"splitcomplex{n}")
    # One character/eigenvalue is -1, the others +1: determinant -1.
    projector = cls.one()
    for generator in range(n):
        projector = projector*(cls.one()+cls.basis()[1 << generator])/2
    value = cls.one() - 2*projector
    norm = value.norm()
    assert type(norm) is complex
    assert norm == pytest.approx(np.exp(1j*np.pi/cls.N))


@pytest.mark.parametrize("n,parameters,central_index", [
    (3, (1, 1, -1), 7), (4, (1, 1, 1, -1), 11),
])
def test_clifford_negative_norm(n, parameters, central_index):
    cls = getattr(general, f"clifford{n}")(*parameters)
    basis = cls.basis()
    projector = (cls.one()+basis[central_index])/2
    if n == 3:
        value = basis[1]*projector + (cls.one()-projector)
    else:
        projector = projector*(cls.one()+basis[1])/2
        value = cls.one()-2*projector
    assert value.norm() == pytest.approx(np.exp(1j*np.pi/4))
    assert type(value.norm()) is complex


@pytest.mark.parametrize("name,n", CONCRETE)
def test_specialization_structure(name, n):
    cls = getattr(algebras, name)
    basis, coeffs = cls.basis(), cls.coeffs()
    assert basis.shape == (n,) and coeffs.shape == (n, n)
    np.testing.assert_array_equal(cls.to_float(basis), np.eye(n))
    for i in range(n):
        for j in range(n):
            assert_components(cls, coeffs[i, j], basis[i]*basis[j])
    x = cls(*np.linspace(0.1, 0.3, n))
    np.testing.assert_allclose(
        np.polynomial.Polynomial.fromroots(x.eivals()).coef,
        x.minpoly(), atol=1e-10,
    )


def test_parameterized_ufunc_reads_current_parameters():
    cls = general.cayleydickson2(-1, -1)
    x = cls(1, 0.2, 0.1, 0.05)
    values = np.array([x], dtype=cls)
    before = np.exp(values)
    cls.K = [1, 0]
    assert_components(cls, np.exp(values)[0], x.exp())
    assert not np.allclose(cls.to_float(before), cls.to_float(np.exp(values)))
    assert_components(cls, np.reciprocal(values)[0], x.inv())


def test_todo_stubs_have_defined_nan_results():
    # Unfinished operations in the remaining general families are safe to call.
    for cls, operation in ((general.cayleydickson3(-1, -1, -1), "exp"),
                           (general.clifford4(-1, -1, -1, -1), "exp")):
        values = np.array([cls.one()], dtype=cls)
        assert np.isnan(cls.to_float(getattr(np, operation)(values))).all()


def test_common_c_api_compiles_and_runs(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    headers = sorted((root/"src/algpy/algebras").rglob("*.h"))
    lines = ['#include "src/algpy/fields/complex.h"']
    lines += [f'#include "{p.relative_to(root)}"' for p in headers]
    lines.append("int main(void) {")
    for p in headers:
        prefix = p.stem
        parameterized = "general" in p.parts
        dim = 2**int(prefix[-1]) if parameterized else dict(CONCRETE).get(
            prefix, 2 if prefix == "complex" else 0)
        assert dim
        degree = (2 if prefix.startswith("cayleydickson") or prefix in
                  ("complex", "complex1", "dual1", "splitcomplex1", "quaternion", "octonion", "sedenion", "clifford1", "clifford2")
                  else 4 if prefix.startswith("clifford") else dim)
        k = "k, " if parameterized else ""
        lines += ["{", f"{prefix} a = {prefix}_one(), b = {prefix}_zero();",
                  f"{prefix} basis[{dim}], coeffs[{dim*dim}];",
                  f"double poly[{degree+1}]; complex ev[{degree}];",
                  f"double left_matrix[{dim*dim}], right_matrix[{dim*dim}];"]
        if parameterized:
            order = int(prefix[-1])
            lines.append(f"vec{order} k = {{{{" + ", ".join(["-1"]*order) + "}};")
        lines += [f"(void){prefix}_make("+", ".join(["0"]*dim)+");",
                  f"(void){prefix}_from_re(1);", f"(void){prefix}_to_re(a);",
                  f"{prefix}_basis(basis);", f"{prefix}_coeffs({k}coeffs);",
                  f"{prefix}_minpoly({k}a, poly);", f"{prefix}_eivals({k}a, ev);",
                  f"{prefix}_lmat({k}a, left_matrix);", f"{prefix}_rmat({k}a, right_matrix);",
                  f"for (int i = 0; i < {dim}; i++) {{",
                  f"for (int j = 0; j < {dim}; j++) {{",
                  f"if (left_matrix[i*{dim}+j] != (i == j) || right_matrix[i*{dim}+j] != (i == j)) return 1;",
                  "}", "}"]
        for op in ("re", "im", "is_real", "pos", "neg", "conjugate"):
            lines.append(f"(void){prefix}_{op}(a);")
        for op in ("add", "sub", "equal", "not_equal"):
            lines.append(f"(void){prefix}_{op}(a, b);")
        for op in ("mul", "div", "pow"):
            lines.append(f"(void){prefix}_{op}({k}a, a);")
        norm_type = "complex" if parameterized or prefix.startswith("splitcomplex") else "double"
        lines.append(f"{norm_type} norm = {prefix}_norm({k}a); (void)norm;")
        for op in ("inv", "square", "exp", "expm1", "log", "log1p", "sqrt",
                   "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
                   "asinh", "acosh", "atanh"):
            lines.append(f"(void){prefix}_{op}({k}a);")
        for op in ("add", "sub", "mul", "div", "pow"):
            params = "" if op in ("add", "sub") else k
            lines += [f"(void){prefix}_r{op}({params}1, a);",
                      f"(void){prefix}_l{op}({params}a, 1);"]
        lines += [f"(void){prefix}_mix({k}a, a, 0.5);", "}"]
        if prefix in ("complex3", "complex4"):
            lines += [
                f"if (!(fabs({prefix}_rpow(2, {prefix}_from_re(3)).data[0] - 8) < 1e-12)) return 1;",
                f"if (!(fabs({prefix}_lpow({prefix}_from_re(2), 3).data[0] - 8) < 1e-12)) return 1;",
                f"if (!(fabs({prefix}_mix({prefix}_from_re(2), {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12)) return 1;",
            ]
    lines.append("return 0; }")
    source = tmp_path/"api.c"
    source.write_text("\n".join(lines))
    binary = tmp_path/"api"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)
