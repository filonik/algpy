"""Jet algebras R[e]/(e^(N+1)) and C[e]/(e^(N+1)) as numpy user-defined dtypes.

jetN has N+1 coefficients, those of 1, e, ..., e^N. Applying a function to
x + e gives its normalized derivatives f^(k)(x)/k!, k = 0..N. complex_jetN
has N+1 complex coefficients, 2(N+1) real components, and gives the
normalized complex derivatives of holomorphic functions.
"""
from ._jet1 import jet1
from ._jet2 import jet2
from ._jet3 import jet3
from ._jet4 import jet4
from ._jet5 import jet5
from ._complex_jet1 import complex_jet1
from ._complex_jet2 import complex_jet2
from ._complex_jet3 import complex_jet3
from ._complex_jet4 import complex_jet4
from ._complex_jet5 import complex_jet5

__all__ = [
    "jet1", "jet2", "jet3", "jet4", "jet5",
    "complex_jet1", "complex_jet2", "complex_jet3", "complex_jet4", "complex_jet5",
]
