"""Independent bitmask reference for exterior lattices.

Basis blade ``a`` is the product of the generators whose bits are set in
``a``, in increasing order. This is the coefficient order of the C types.
"""
import numpy as np
from algpy import lattices


def lattice_type(name):
    return getattr(lattices, name)


def generators(cls):
    return cls.N.bit_length()-1


def reorder_sign(a, b):
    """Sign of e_a e_b reordered into canonical order."""
    a >>= 1
    swaps = 0
    while a:
        swaps += bin(a & b).count("1")
        a >>= 1
    return -1 if swaps & 1 else 1


def join(x, y):
    n = len(x)
    z = np.zeros(n)
    for a in range(n):
        for b in range(n):
            if not a & b:
                z[a | b] += reorder_sign(a, b)*x[a]*y[b]
    return z


def rdual(x):
    """rdual(e_a) join e_a = top."""
    n = len(x)
    full = n-1
    z = np.zeros(n)
    for a in range(n):
        z[full ^ a] = reorder_sign(full ^ a, a)*x[a]
    return z


def ldual(x):
    """e_a join ldual(e_a) = top."""
    n = len(x)
    full = n-1
    z = np.zeros(n)
    for a in range(n):
        z[full ^ a] = reorder_sign(a, full ^ a)*x[a]
    return z


def meet(x, y):
    return rdual(join(ldual(x), ldual(y)))


def grade(x, k):
    return np.array([c if bin(a).count("1") == k else 0 for a, c in enumerate(x)])


def samples(cls, count=3, seed=0):
    rng = np.random.default_rng(seed + cls.N)
    return rng.normal(size=(count, cls.N))
