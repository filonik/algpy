"""Run a benchmark entry point and check its output contract."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]

def check_benchmark_entry_point(tmp_path, algebra, reference):
    pytest.importorskip("scipy")
    pytest.importorskip("threadpoolctl")
    command = [sys.executable, str(ROOT/"scripts"/f"benchmark_{algebra}.py"),
               "--output", str(tmp_path), "--cases", "generic-negative",
               "--functions", "exp", "log", "--samples", "2",
               "--repeats", "1", "--target", ".001"]
    if algebra == "cayleydickson4":
        command += ["--reference", reference]
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
    results = json.loads((tmp_path/"results.json").read_text())
    metadata = results["metadata"]
    assert metadata["algebra"] == algebra
    assert metadata["matrix_shape"] == ([2, 2] if reference == "slice" else [16, 16])
    assert metadata["reference_kind"] == ("quadratic-slice" if reference == "slice" else
        "identity-action" if algebra == "cayleydickson4" else "full-matrix")
    assert metadata["output_kind"] == ("coefficients" if algebra == "cayleydickson4" else "matrix")
    assert metadata["samples"] == 2
    case, = results["cases"]
    assert case["parameters"] == [-1, -2, -3, -4]
    assert [f["name"] for f in case["functions"]] == ["exp", "log"]
    for function in case["functions"]:
        assert function["max_abs_error"] < 2e-10
        assert function["scipy_batch_mode"] in {"native", "loop"}
        for timing in function["timings"].values():
            assert np.isfinite(timing["median_us"]) and timing["median_us"] > 0
            assert len(timing["samples_us"]) == 1
            assert timing["loops_per_repeat"] >= 1
    with (tmp_path/"results.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))
    assert len(rows) == 2
    for row in rows:
        assert row["algebra"] == algebra
        assert float(row["scalar_speedup"]) > 0
        assert float(row["batch_speedup"]) > 0
        if algebra == "cayleydickson4":
            assert "matrix_output_speedup" not in row
            assert float(row["coefficients_output_speedup"]) > 0
        else:
            assert float(row["matrix_output_speedup"]) > 0
