"""Independent complex character coordinates used by specialization tests."""
import itertools

import numpy as np


FUNCTIONS = [
    ("exp", "exp"), ("expm1", "expm1"), ("log", "log"), ("log1p", "log1p"),
    ("sqrt", "sqrt"), ("cos", "cos"), ("sin", "sin"), ("tan", "tan"),
    ("arccos", "acos"), ("arcsin", "asin"), ("arctan", "atan"),
    ("cosh", "cosh"), ("sinh", "sinh"), ("tanh", "tanh"),
    ("arccosh", "acosh"), ("arcsinh", "asinh"), ("arctanh", "atanh"),
]


def character_matrix(cls):
    order = cls.N.bit_length()-1
    # Independent, nonrecursive representation: send i0 to i, and each other
    # generator to +/-i. Basis products are indexed by their generator bits.
    matrix = []
    for signs in itertools.product([-1, 1], repeat=order-1):
        row = []
        for basis in range(cls.N):
            value = 1j**basis.bit_count()
            for generator, sign in enumerate(signs, 1):
                if basis & (1 << generator):
                    value *= sign
            row.append(value)
        matrix.append(row)
    return np.array(matrix)


def reconstruct(matrix, components):
    return (components @ matrix.conj()).real/len(matrix)
