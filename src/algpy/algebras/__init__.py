"""Custom number types for numpy, implemented in Cython."""
from .general.multicomplex import multicomplex, multicomplex1, multicomplex2, multicomplex3, multicomplex4
from .general.clifford import clifford, clifford1, clifford2
from .general.cayleydickson import (
    cayleydickson,
    cayleydickson1,
    cayleydickson2,
    cayleydickson3,
    cayleydickson4,
)
from .quaternion import quaternion
from .complex import complex1
from .dual import dual1
from .splitcomplex import splitcomplex1

__all__ = [
    "multicomplex",
    "multicomplex1",
    "multicomplex2",
    "multicomplex3",
    "multicomplex4",
    "clifford",
    "clifford1",
    "clifford2",
    "cayleydickson",
    "cayleydickson1",
    "cayleydickson2",
    "cayleydickson3",
    "cayleydickson4",
    "quaternion",
    "complex1",
    "dual1",
    "splitcomplex1",
]
