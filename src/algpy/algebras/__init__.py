"""Custom number types for NumPy, implemented in C.

The general families live in ``algpy.algebras.general``.
"""
from .cayleydickson import real, quaternion, octonion, sedenion
from .multicomplex import (
    complex1,
    complex2,
    complex3,
    complex4,
    dual1,
    dual2,
    dual3,
    dual4,
    splitcomplex1,
    splitcomplex2,
    splitcomplex3,
    splitcomplex4,
)

# The complex C algebra and complex1 share one concrete Python type.
complex = complex1

__all__ = [
    "real",
    "complex",
    "quaternion",
    "octonion",
    "sedenion",
    "complex1",
    "complex2",
    "complex3",
    "complex4",
    "dual1",
    "dual2",
    "dual3",
    "dual4",
    "splitcomplex1",
    "splitcomplex2",
    "splitcomplex3",
    "splitcomplex4",
]
