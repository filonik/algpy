"""Custom number types for numpy, implemented in C.

The general families live in ``algpy.algebras.general``.
"""
from .cayleydickson import quaternion
from .multicomplex import complex1, complex2, dual1, dual2, splitcomplex1, splitcomplex2

__all__ = [
    "quaternion",
    "complex1",
    "complex2",
    "dual1",
    "dual2",
    "splitcomplex1",
    "splitcomplex2",
]
