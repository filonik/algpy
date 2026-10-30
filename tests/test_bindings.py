"""Check the complete C API and extension loading across the library."""
from pathlib import Path
import importlib
import json
import shutil
import subprocess
import sys
import numpy as np
import pytest
from algebra_reference import CONCRETE, FAMILIES
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
        for op in ("mul", "div", "pow", "lpow", "rpow"):
            lines.append(f"(void){prefix}_{op}({k}a, a);")
        lines += [f"if ({prefix}_ipow({k}{prefix}_from_re(-2), 3).data[0] != -8) return 1;",
                  f"if ({prefix}_ipow({k}a, INT64_MIN).data[0] != 1) return 1;",
                  f"if ({prefix}_ipow({k}{prefix}_zero(), 0).data[0] != 1) return 1;"]
        norm_type = "complex" if parameterized or prefix.startswith("splitcomplex") else "double"
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
        family = "cayleydickson" if name in ["real", "quaternion", "octonion", "sedenion"] else "multicomplex"
        fullname = f"algpy.algebras.{family}._{name}"
        path = importlib.import_module(fullname).__file__
        modules.append((fullname, path, name, 0))
    for family in FAMILIES:
        for order in range(1, 5):
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
    assert cls.N in [1, 2, 4, 8, 16]
    assert inspect.isbuiltin(cls.to_float) and inspect.isbuiltin(cls.from_float)
    assert inspect.isbuiltin(cls.one().lpow) and inspect.isbuiltin(cls.one().rpow)
    data = np.arange(cls.N, dtype=float)
    np.testing.assert_array_equal(cls.to_float(cls.from_float(data)), data)
assert "algpy.algebras" not in sys.modules
assert "algpy.algebras._conversions" not in sys.modules
'''
    subprocess.run([sys.executable, "-c", script, json.dumps(modules)],
                   check=True, capture_output=True, text=True)
