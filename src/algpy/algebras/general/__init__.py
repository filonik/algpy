"""General algebra families: multicomplex, Cayley-Dickson and Clifford."""
from .multicomplex import multicomplex, multicomplex1, multicomplex2, multicomplex3, multicomplex4
from .clifford import clifford, clifford1, clifford2, clifford3, clifford4
from .cayleydickson import (
    cayleydickson,
    cayleydickson1,
    cayleydickson2,
    cayleydickson3,
    cayleydickson4,
)

__all__ = [
    "multicomplex",
    "multicomplex1",
    "multicomplex2",
    "multicomplex3",
    "multicomplex4",
    "clifford",
    "clifford1",
    "clifford2",
    "clifford3",
    "clifford4",
    "cayleydickson",
    "cayleydickson1",
    "cayleydickson2",
    "cayleydickson3",
    "cayleydickson4",
]
