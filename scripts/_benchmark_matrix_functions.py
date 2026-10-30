"""Shared runner for the 16-component algebra benchmark entry points.

Install dependencies with python -m pip install -e ".[benchmark]".
Run benchmark_multicomplex4.py, benchmark_clifford4.py or
benchmark_cayleydickson4.py. The process uses one BLAS thread.
Construction and correctness checks are untimed. Results go to benchmark-results/.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import timeit

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import scipy
from scipy import linalg
from threadpoolctl import threadpool_info, threadpool_limits

from algpy.algebras import general


FUNCTIONS = [
    ("exp", linalg.expm), ("log", linalg.logm), ("sqrt", linalg.sqrtm),
    ("sin", linalg.sinm), ("cos", linalg.cosm), ("tan", linalg.tanm),
    ("sinh", linalg.sinhm), ("cosh", linalg.coshm), ("tanh", linalg.tanhm),
]
MULTICOMPLEX_CASES = {
    "generic-negative": (-1, -2, -3, -4),
    "complex": (-1, -1, -1, -1),
    "mixed": (-2, 0, 0.5, -1.5),
    "split": (1, 1, 1, 1),
    "dual": (0, 0, 0, 0),
}
GENERAL_CASES = {
    "generic-negative": (-1, -2, -3, -4),
    "unit-negative": (-1, -1, -1, -1),
    "mixed": (-2, 0, 0.5, -1.5),
    "unit-positive": (1, 1, 1, 1),
    "nilpotent": (0, 0, 0, 0),
}
ROOT = Path(__file__).resolve().parents[1]


def measure(function, size, repeats, target):
    timer = timeit.Timer(function)
    function()
    elapsed = timer.timeit(number=1)
    number = max(1, int(target/max(elapsed, 1e-9)))
    samples = np.array(timer.repeat(repeat=repeats, number=number))/(number*size)*1e6
    return {
        "median_us": float(np.median(samples)),
        "q25_us": float(np.percentile(samples, 25)),
        "q75_us": float(np.percentile(samples, 75)),
        "samples_us": samples.tolist(),
        "loops_per_repeat": number,
    }


def input_coefficients(cls, count):
    rng = np.random.default_rng(106)
    coefficients = rng.uniform(-1, 1, (count, cls.N))
    coefficients[:, 0] = 0
    for row in coefficients:
        # Keep every eigenvalue in a disk of radius <= 0.2 about the center.
        matrix_norm = np.linalg.norm(cls(*row).lmat(), 1)
        row *= 0.2/max(1, matrix_norm)
    return coefficients


def scipy_batch_evaluator(function, matrices):
    try:
        reference = function(matrices)
    except ValueError:
        # Older SciPy releases require one matrix per call. Keep batch output
        # shapes the same, and record the fallback in the results.
        evaluate = lambda: np.stack([function(a) for a in matrices])
        return evaluate, evaluate(), "loop"
    return lambda: function(matrices), reference, "native"


def slice_coefficients(matrix, vector):
    """Map f([[w,q],[1,w]]) back from span{1,v} to the 16 coefficients."""
    result = vector*matrix[..., 1, 0, None]
    result[..., 0] = matrix[..., 0, 0]
    return result


def scipy_evaluators(function, matrices, identity_action=False, vectors=None):
    batch, reference, mode = scipy_batch_evaluator(function, matrices)
    if vectors is not None:
        scalar = lambda: [slice_coefficients(function(a), v) for a, v in zip(matrices, vectors)]
        return scalar, lambda: slice_coefficients(batch(), vectors), slice_coefficients(reference, vectors), mode
    if identity_action:
        scalar = lambda: [function(a)[:, 0] for a in matrices]
        return scalar, lambda: batch()[..., :, 0], reference[..., :, 0], mode
    return lambda: [function(a) for a in matrices], batch, reference, mode


def benchmark(args, algebra, cases):
    coefficients_output = algebra == "cayleydickson4"
    reference_kind = ("quadratic-slice" if args.reference == "slice" else
                      "identity-action" if coefficients_output else "full-matrix")
    matrix_dimension = 2 if args.reference == "slice" else 16
    output_key = ("algebra_scalar_coefficients_output" if coefficients_output else
                  "algebra_scalar_matrix_output")
    speedup_key = "coefficients_output_speedup" if coefficients_output else "matrix_output_speedup"
    cpuinfo = Path("/proc/cpuinfo")
    cpu = next((line.split(":", 1)[1].strip() for line in
                (cpuinfo.read_text().splitlines() if cpuinfo.exists() else [])
                if line.startswith("model name")), platform.processor())
    try:
        revision = subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    metadata = {
        "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "cpu": cpu,
        "git_commit": revision,
        "threadpools": threadpool_info(), "samples": args.samples,
        "repeats": args.repeats, "target_seconds": args.target,
        "seed": 106, "algebra": algebra, "reference_kind": reference_kind,
        "matrix_shape": [matrix_dimension, matrix_dimension],
        "output_kind": "coefficients" if coefficients_output else "matrix",
    }
    results = {"metadata": metadata, "cases": []}
    summary = []
    print(f"{algebra}: {matrix_dimension}x{matrix_dimension} SciPy {reference_kind} reference", flush=True)
    for label, parameters in cases.items():
        if label not in args.cases:
            continue
        cls = getattr(general, algebra)(*parameters)
        algebra_output = cls.to_float if coefficients_output else lambda value: value.lmat()
        coefficients = input_coefficients(cls, args.samples)
        conversion_values = [cls(*row) for row in coefficients]
        conversion = measure(lambda: [v.lmat() for v in conversion_values],
                             args.samples, args.repeats, args.target)
        case = {"name": label, "parameters": parameters,
                "lmat": conversion, "functions": []}
        for name, matrix_function in FUNCTIONS:
            if name not in args.functions:
                continue
            rows = coefficients.copy()
            rows[:, 0] = 1.25 if name in {"log", "sqrt"} else 0.35
            values = cls.from_float(rows)
            scalars = list(values)
            matrices = np.stack([v.lmat() for v in scalars])
            vectors = None
            if args.reference == "slice":
                vectors = rows.copy()
                vectors[:, 0] = 0
                # q=<v*v>_0 is a row contraction of the multiplication matrix.
                # Computing it from the centred element avoids subtracting w^2
                # from its norm. Input construction and reconstruction checks
                # are outside the timers.
                q = np.einsum("ni,ni->n", matrices[:, 0, 1:], vectors[:, 1:])
                matrices = np.zeros((args.samples, 2, 2))
                matrices[:, 0, 0] = matrices[:, 1, 1] = rows[:, 0]
                matrices[:, 0, 1], matrices[:, 1, 0] = q, 1
            function = getattr(np, name)
            scipy_scalar, scipy_batch, reference, batch_mode = scipy_evaluators(
                matrix_function, matrices, coefficients_output, vectors)
            actual = np.stack([algebra_output(v) for v in function(values)])
            np.testing.assert_allclose(actual, reference, rtol=2e-10, atol=2e-12,
                                       equal_nan=False)
            maximum_error = float(np.max(np.abs(actual-reference)))
            timings = {
                "algebra_scalar": measure(lambda: [function(v) for v in scalars],
                                          args.samples, args.repeats, args.target),
                "scipy_scalar": measure(scipy_scalar, args.samples, args.repeats, args.target),
                output_key: measure(
                    lambda: [algebra_output(function(v)) for v in scalars],
                    args.samples, args.repeats, args.target),
                "algebra_batch": measure(lambda: function(values),
                                         args.samples, args.repeats, args.target),
                "scipy_batch": measure(scipy_batch, args.samples, args.repeats, args.target),
            }
            medians = {key: value["median_us"] for key, value in timings.items()}
            record = {"name": name, "max_abs_error": maximum_error,
                      "scipy_batch_mode": batch_mode, "timings": timings}
            case["functions"].append(record)
            row = {"algebra": algebra, "reference_kind": reference_kind,
                   "matrix_dimension": matrix_dimension,
                   "case": label, "parameters": str(parameters), "function": name,
                   "scipy_batch_mode": batch_mode,
                   **{key+"_us": value for key, value in medians.items()},
                   "scalar_speedup": medians["scipy_scalar"]/medians["algebra_scalar"],
                   speedup_key: medians["scipy_scalar"]/medians[output_key],
                   "batch_speedup": medians["scipy_batch"]/medians["algebra_batch"],
                   "lmat_us": conversion["median_us"], "max_abs_error": maximum_error}
            summary.append(row)
            print(f"{label:17} {name:5} algebra {medians['algebra_scalar']:8.3f} us; "
                  f"SciPy {medians['scipy_scalar']:9.3f} us; "
                  f"speedup {row['scalar_speedup']:7.1f}x; "
                  f"{metadata['output_kind']}-output {medians[output_key]:8.3f} us; "
                  f"batch {row['batch_speedup']:7.1f}x", flush=True)
        results["cases"].append(case)
        (args.output/"results.json").write_text(json.dumps(results, indent=2)+"\n", encoding="utf-8")
    with (args.output/"results.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print("Maximum reference error:", max(row["max_abs_error"] for row in summary), flush=True)
    print("Results:", args.output.resolve(), flush=True)


def main(algebra):
    cases = MULTICOMPLEX_CASES if algebra == "multicomplex4" else GENERAL_CASES
    parser = argparse.ArgumentParser(description=f"Compare {algebra} with SciPy matrix functions.")
    parser.add_argument("--output", type=Path, help="Directory for results.json and results.csv")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--target", type=float, default=0.025,
                        help="Target seconds per timed repeat")
    parser.add_argument("--cases", nargs="+", choices=list(cases), default=list(cases))
    parser.add_argument("--functions", nargs="+", choices=[name for name, _ in FUNCTIONS],
                        default=[name for name, _ in FUNCTIONS])
    if algebra == "cayleydickson4":
        parser.add_argument("--reference", choices=["matrix", "slice"], default="matrix",
                            help="SciPy 16x16 matrix action or 2x2 quadratic slice")
    else:
        parser.set_defaults(reference="matrix")
    args = parser.parse_args()
    if args.samples < 1 or args.repeats < 1 or not np.isfinite(args.target) or args.target <= 0:
        parser.error("samples, repeats and target must be positive and finite")
    if args.output is None:
        args.output = ROOT/"benchmark-results"/algebra
        if args.reference == "slice":
            args.output /= "slice"
    args.output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1, user_api="blas"):
        benchmark(args, algebra, cases)
