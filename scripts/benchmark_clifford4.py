"""Compare clifford4 elementary functions with SciPy matrix functions.

Run python scripts/benchmark_clifford4.py --help for options.
"""
from _benchmark_matrix_functions import main


if __name__ == "__main__":
    main("clifford4")
