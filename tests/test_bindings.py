"""Check the complete C API and extension loading across the library."""
from pathlib import Path
import importlib
import json
import shutil
import subprocess
import sys
import numpy as np
import pytest
from algebra_reference import CONCRETE, FAMILIES, FIXED_METRICS, COMPLEX_NORMS
from algpy import algebras
from algpy.algebras import general


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
        jet = prefix.startswith("jet")
        complex_jet = prefix.startswith("complex_jet")
        generators = 1 if jet else 2 if complex_jet else dim.bit_length()-1
        degree = (8 if prefix in ("cga3", "clifford5", "pga4") else 2 if prefix.startswith("cayleydickson") or prefix in
                  ("complex", "complex1", "dual1", "splitcomplex1", "quaternion", "octonion", "sedenion", "trigintaduonion", "clifford1", "clifford2", "pga1")
                  else 4 if prefix.startswith("clifford") or prefix in FIXED_METRICS else dim)
        k = "k, " if parameterized else ""
        lines += ["{", f"{prefix} a = {prefix}_one(), b = {prefix}_zero();",
                  f"{prefix} basis[{dim}], coeffs[{dim*dim}];",
                  f"{prefix} gens[{generators+1}];",
                  f"double poly[{degree+1}]; complex ev[{degree}];",
                  f"double left_matrix[{dim*dim}], right_matrix[{dim*dim}];"]
        if parameterized:
            order = int(prefix[-1])
            lines.append(f"vec{order} k = {{{{" + ", ".join(["-1"]*order) + "}};")
        lines += [f"double input[{dim+2}];",
                  f"for (int i = 0; i < {dim+2}; i++) input[i] = i + .25;",
                  "input[1] = -0.0;",
                  f"{prefix} (*make)(const double *) = {prefix}_make;",
                  f"{prefix} made = make(input + 1);",
                  "if (memcmp(made.data, input + 1, sizeof(made.data)) != 0) return 1;",
                  "input[1] = 7;",
                  "if (!signbit(made.data[0])) return 1;",
                  "made = make(made.data);",
                  "if (!signbit(made.data[0])) return 1;",
                  f"(void){prefix}_from_re(1);", f"(void){prefix}_to_re(a);",
                  f"{prefix}_basis(basis);", f"{prefix}_coeffs({k}coeffs);",
                  f"{prefix}_minpoly({k}a, poly);", f"{prefix}_eivals({k}a, ev);",
                  f"{prefix}_lmat({k}a, left_matrix);", f"{prefix}_rmat({k}a, right_matrix);",
                  f"for (int i = 0; i < {dim}; i++) {{",
                  f"for (int j = 0; j < {dim}; j++) {{",
                  f"if (left_matrix[i*{dim}+j] != (i == j) || right_matrix[i*{dim}+j] != (i == j)) return 1;",
                  "}", "}"]
        lines += [f"gens[{generators}] = {prefix}_from_re(7);",
                  f"{prefix}_gens(gens);",
                  f"if (!{prefix}_equal(gens[{generators}], {prefix}_from_re(7))) return 1;",
                  f"for (int i = 0; i < {generators}; i++) {{",
                  f"for (int j = 0; j < {dim}; j++) {{",
                  "if (gens[i].data[j] != (j == (1 << i))) return 1;",
                  "}", "}"]
        if not generators:
            lines.append(f"{prefix}_gens(NULL);")
        for op in ("re", "im", "is_real", "pos", "neg", "conjugate"):
            lines.append(f"(void){prefix}_{op}(a);")
        lines += [f"{prefix} graded = {prefix}_zero();",
                  f"for (int i = 0; i < {dim}; i++) input[i] = i + .5;",
                  f"made = {prefix}_make(input);"]
        # Grade k keeps the basis indices with k generator bits, or for jets
        # the coefficient of e^k.
        grade_of = "i" if jet else "i/2" if complex_jet else "__builtin_popcount(i)"
        for order in range(dim if jet else dim//2 if complex_jet else generators + 1):
            lines += [f"graded = {prefix}_add(graded, {prefix}_grade{order}(made));",
                      f"for (int i = 0; i < {dim}; i++) {{",
                      f"if ({prefix}_grade{order}(made).data[i] != ({grade_of} == {order} ? made.data[i] : 0)) return 1;",
                      "}"]
        lines.append(f"if (!{prefix}_equal(graded, made)) return 1;")
        centers = {"clifford1": "clifford1", "clifford3": "clifford1", "clifford5": "clifford1",
                   "pga2": "dual1", "pga4": "dual1", "cga1": "splitcomplex1", "cga3": "complex1"}
        if prefix in centers:
            count = {2: 2, 8: 3, 32: 5}[dim]
            lines += [f"{centers[prefix]} central[{count}];",
                      f"{prefix}_center_minpoly({k}made, central);",
                      f"if (central[{count - 1}].data[0] != 1 || central[{count - 1}].data[1] != 0) return 1;"]
        for op in ("add", "sub", "equal", "not_equal"):
            lines.append(f"(void){prefix}_{op}(a, b);")
        for op in ("mul", "div", "pow", "lpow", "rpow"):
            lines.append(f"(void){prefix}_{op}({k}a, a);")
        lines += [f"if ({prefix}_ipow({k}{prefix}_from_re(-2), 3).data[0] != -8) return 1;",
                  f"if ({prefix}_ipow({k}a, INT64_MIN).data[0] != 1) return 1;",
                  f"if ({prefix}_ipow({k}{prefix}_zero(), 0).data[0] != 1) return 1;"]
        norm_type = "complex" if parameterized or prefix in COMPLEX_NORMS else "double"
        lines.append(f"{norm_type} norm = {prefix}_norm({k}a); (void)norm;")
        for op in ("inv", "square", "exp", "expm1", "log", "log1p", "sqrt",
                   "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
                   "asinh", "acosh", "atanh"):
            lines.append(f"(void){prefix}_{op}({k}a);")
        for op in ("add", "sub", "mul", "div", "pow"):
            params = "" if op in ("add", "sub") else k
            suffix = "pows" if op == "pow" else op
            lines += [f"(void){prefix}_r{suffix}({params}1, a);",
                      f"(void){prefix}_l{suffix}({params}a, 1);"]
        lines += [f"(void){prefix}_mix({k}a, a, 0.5);", "}"]
        if prefix in ("complex3", "complex4"):
            lines += [
                f"if (!(fabs({prefix}_rpows(2, {prefix}_from_re(3)).data[0] - 8) < 1e-12)) return 1;",
                f"if (!(fabs({prefix}_lpows({prefix}_from_re(2), 3).data[0] - 8) < 1e-12)) return 1;",
                f"if (!(fabs({prefix}_mix({prefix}_from_re(2), {prefix}_from_re(8), 0.5).data[0] - 4) < 1e-12)) return 1;",
            ]
    lines.append("return 0; }")
    source = tmp_path/"api.c"
    source.write_text("\n".join(lines))
    binary = tmp_path/"api"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)

def test_extensions_are_complete_without_python_package_initialization():
    modules = []
    for name, _ in CONCRETE:
        family = ("clifford" if name in FIXED_METRICS else
                  "jet" if name.startswith(("jet", "complex_jet")) else
                  "cayleydickson" if name in ["real", "quaternion", "octonion", "sedenion", "trigintaduonion"] else "multicomplex")
        fullname = f"algpy.algebras.{family}._{name}"
        path = importlib.import_module(fullname).__file__
        modules.append((fullname, path, name, 0))
    for family in FAMILIES:
        for order in range(1, 6):
            name = f"{family}{order}"
            fullname = f"algpy.algebras.general.{family}._{name}"
            path = importlib.import_module(fullname).__file__
            modules.append((fullname, path, name, order))
    script = '''
import importlib.util
import inspect
import json
import sys
import numpy as np
for fullname, path, name, order in json.loads(sys.argv[1]):
    spec = importlib.util.spec_from_file_location(fullname, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[fullname] = module
    spec.loader.exec_module(module)
    factory = getattr(module, name)
    cls = factory(*([-1]*order)) if order else factory
    jet = name.startswith("jet")
    complex_jet = name.startswith("complex_jet")
    assert cls.N in ([2, 3, 4, 5, 6] if jet else [4, 6, 8, 10, 12] if complex_jet else [1, 2, 4, 8, 16, 32])
    assert inspect.isbuiltin(cls.to_float) and inspect.isbuiltin(cls.from_float)
    assert inspect.isbuiltin(cls.one().lpow) and inspect.isbuiltin(cls.one().rpow)
    assert inspect.isbuiltin(cls.gens)
    assert cls.gens().shape == ((1,) if jet else (2,) if complex_jet else (cls.N.bit_length()-1,))
    assert cls.gens().dtype == cls.dtype
    data = np.arange(cls.N, dtype=float)
    np.testing.assert_array_equal(cls.to_float(cls.from_float(data)), data)
assert "algpy.algebras" not in sys.modules
assert "algpy.algebras._conversions" not in sys.modules
'''
    subprocess.run([sys.executable, "-c", script, json.dumps(modules)],
                   check=True, capture_output=True, text=True)


def test_lattice_c_api_compiles_and_runs(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    headers = sorted((root/"src/algpy/lattices").rglob("*.h"))
    assert [p.stem for p in headers] == [f"exterior{n}" for n in range(1, 6)]
    lines = [f'#include "{p.relative_to(root)}"' for p in headers]
    lines.append("int main(void) {")
    for p in headers:
        prefix = p.stem
        n = int(prefix[-1])
        dim = 2**n
        lines += ["{", f"{prefix} a = {prefix}_bottom(), b = {prefix}_top(), z = {prefix}_zero();",
                  f"{prefix} basis[{dim}], gens[{n+1}];",
                  f"double input[{dim}];",
                  f"for (int i = 0; i < {dim}; i++) input[i] = i + .25;",
                  f"{prefix} made = {prefix}_make(input);",
                  "if (memcmp(made.data, input, sizeof(made.data)) != 0) return 1;",
                  f"if (!{prefix}_equal({prefix}_join(a, made), made)) return 1;",
                  f"if (!{prefix}_equal({prefix}_meet(b, made), made)) return 1;",
                  f"if (!{prefix}_equal({prefix}_rdual(a), b)) return 1;",
                  f"if (!{prefix}_equal({prefix}_ldual(b), a)) return 1;",
                  f"if (!{prefix}_equal({prefix}_ldual({prefix}_rdual(made)), made)) return 1;",
                  f"if ({prefix}_not_equal({prefix}_join(z, made), z)) return 1;",
                  f"{prefix}_basis(basis);",
                  f"gens[{n}] = {prefix}_from_re(7);",
                  f"{prefix}_gens(gens);",
                  f"if (!{prefix}_equal(gens[{n}], {prefix}_from_re(7))) return 1;",
                  f"for (int i = 0; i < {n}; i++) {{",
                  f"if (!{prefix}_equal(gens[i], basis[1 << i])) return 1;",
                  "a = " + f"{prefix}_join(a, gens[i]);",
                  "}",
                  f"if (!{prefix}_equal(a, b)) return 1;",
                  f"(void){prefix}_to_re(made); (void){prefix}_is_real(made); (void){prefix}_norm2(made);",
                  f"(void){prefix}_pos(made); (void){prefix}_neg(made);",
                  f"(void){prefix}_add(made, made); (void){prefix}_sub(made, made);",
                  f"(void){prefix}_rmul(2, made); (void){prefix}_lmul(made, 2); (void){prefix}_ldiv(made, 2);"]
        lines += [f"(void){prefix}_grade{k}(made);" for k in range(n+1)]
        if n % 2:
            lines.append(f"if (!{prefix}_equal({prefix}_dual(made), {prefix}_ldual(made))) return 1;")
        lines.append("}")
    lines.append("return 0; }")
    source = tmp_path/"lattice_api.c"
    source.write_text("\n".join(lines))
    binary = tmp_path/"lattice_api"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)


def test_fixed_metric_lattice_c_api_compiles_and_runs(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler unavailable")
    root = Path(__file__).resolve().parents[1]
    # Fixed geometric algebras and their generator counts.
    algebras = [(f"pga{k}", k + 1) for k in range(1, 5)] + [(f"cga{k}", k + 2) for k in range(1, 4)]
    lines = [f'#include "src/algpy/algebras/clifford/{prefix}.h"' for prefix, _ in algebras]
    lines.append("int main(void) {")
    for prefix, n in algebras:
        exterior = f"exterior{n}"
        dim = 2**n
        lines += ["{", f"double input[{dim}];",
                  f"for (int i = 0; i < {dim}; i++) input[i] = i + .25;",
                  f"{prefix} a = {prefix}_make(input), b = {prefix}_rdual(a);",
                  f"if (!{prefix}_equal({prefix}_bottom(), {prefix}_one())) return 1;",
                  f"if (!{prefix}_equal({prefix}_top(), {exterior}_top())) return 1;",
                  f"if (!{prefix}_equal({prefix}_join(a, b), {exterior}_join(a, b))) return 1;",
                  f"if (!{prefix}_equal({prefix}_meet(a, b), {exterior}_meet(a, b))) return 1;",
                  f"if (!{prefix}_equal({prefix}_join({prefix}_bottom(), a), a)) return 1;",
                  f"if (!{prefix}_equal({prefix}_meet({prefix}_top(), a), a)) return 1;",
                  f"if (!{prefix}_equal({prefix}_ldual(b), a)) return 1;"]
        if n % 2:
            lines.append(f"if (!{prefix}_equal({prefix}_dual(a), {prefix}_ldual(a))) return 1;")
        lines.append("}")
    lines.append("return 0; }")
    source = tmp_path/"fixed_metric_lattice_api.c"
    source.write_text("\n".join(lines))
    binary = tmp_path/"fixed_metric_lattice_api"
    subprocess.run([compiler, "-std=c99", "-Wall", "-Werror", "-I", str(root),
                    str(source), "-lm", "-o", str(binary)], check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], check=True)
