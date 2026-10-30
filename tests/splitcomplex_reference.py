"""Independent real characters used by split-complex specialization tests."""
import numpy as np


FUNCTIONS = [
    ("exp", "exp", -1, 1), ("expm1", "expm1", -1, 1),
    ("log", "log", 0.2, 2), ("log1p", "log1p", -0.8, 1),
    ("sqrt", "sqrt", 0.2, 2),
    ("cos", "cos", -1, 1), ("sin", "sin", -1, 1), ("tan", "tan", -1, 1),
    ("arccos", "acos", -0.9, 0.9), ("arcsin", "asin", -0.9, 0.9),
    ("arctan", "atan", -2, 2),
    ("cosh", "cosh", -1, 1), ("sinh", "sinh", -1, 1), ("tanh", "tanh", -1, 1),
    ("arccosh", "acosh", 1.1, 3), ("arcsinh", "asinh", -2, 2),
    ("arctanh", "atanh", -0.9, 0.9),
]


def character_matrix(cls):
    # Each generator independently maps to -1 or +1. Evaluate each basis
    # product directly, without the recursive splitting used by the C code.
    order = cls.N.bit_length()-1
    matrix = []
    for index in range(cls.N):
        row = []
        for basis in range(cls.N):
            sign = 1
            for generator in range(order):
                if basis & (1 << generator):
                    sign *= 1 if index & (1 << generator) else -1
            row.append(sign)
        matrix.append(row)
    return np.array(matrix)


def reconstruct(matrix, components):
    return components @ matrix/len(matrix)
