"""Exercise benchmark CLIs, reference projections, and saved result formats."""
import csv
import json
from pathlib import Path
import runpy
import subprocess
import sys
import numpy as np
import pytest


pytest.importorskip("scipy")
pytest.importorskip("threadpoolctl")

ROOT = Path(__file__).resolve().parents[1]




@pytest.mark.parametrize("projection", ["matrix", "identity-action", "slice"])
def test_unbatched_scipy_fallback(projection):
    from scipy.linalg import expm
    runner = runpy.run_path(str(ROOT/"scripts"/"_benchmark_matrix_functions.py"))
    matrices = np.array([[[.4, -.04], [1, .4]], [[.8, .01], [1, .8]]])
    vectors = None
    expected = np.stack([expm(a) for a in matrices])
    if projection == "slice":
        vectors = np.zeros((2, 16))
        vectors[0, 1], vectors[1, 3] = .2, -.1
        coefficients = expected[:, 1, 0, None]*vectors
        coefficients[:, 0] = expected[:, 0, 0]
        expected = coefficients
    elif projection == "identity-action":
        expected = expected[:, :, 0]

    def unbatched(matrix):
        if matrix.ndim != 2:
            raise ValueError("This simulated older SciPy accepts only one matrix")
        return expm(matrix)

    scalar, batch, reference, mode = runner["scipy_evaluators"](
        unbatched, matrices, projection != "matrix", vectors)
    assert mode == "loop"
    for actual in [reference, scalar(), batch()]:
        np.testing.assert_allclose(actual, expected, rtol=2e-14, atol=2e-14)
