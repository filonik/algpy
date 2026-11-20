"""Custom number types for NumPy, implemented in C.

The general families live in ``algpy.algebras.general``.
"""
from .clifford import pga1, pga2, pga3, pga4, cga1, cga2, cga3
from .cayleydickson import real, quaternion, octonion, sedenion, trigintaduonion
from .jet import (
    jet1, jet2, jet3, jet4, jet5,
    complex_jet1, complex_jet2, complex_jet3, complex_jet4, complex_jet5,
)
from .multicomplex import (
    complex1,
    complex2,
    complex3,
    complex4,
    complex5,
    dual1,
    dual2,
    dual3,
    dual4,
    dual5,
    splitcomplex1,
    splitcomplex2,
    splitcomplex3,
    splitcomplex4,
    splitcomplex5,
)

# The complex C algebra and complex1 share one concrete Python type.
complex = complex1

__all__ = [
    "pga1", "pga2", "pga3", "pga4", "cga1", "cga2", "cga3",
    "real",
    "complex",
    "quaternion",
    "octonion",
    "sedenion",
    "trigintaduonion",
    "complex1",
    "complex2",
    "complex3",
    "complex4",
    "complex5",
    "dual1",
    "dual2",
    "dual3",
    "dual4",
    "dual5",
    "splitcomplex1",
    "splitcomplex2",
    "splitcomplex3",
    "splitcomplex4",
    "splitcomplex5",
    "jet1",
    "jet2",
    "jet3",
    "jet4",
    "jet5",
    "complex_jet1",
    "complex_jet2",
    "complex_jet3",
    "complex_jet4",
    "complex_jet5",
]
