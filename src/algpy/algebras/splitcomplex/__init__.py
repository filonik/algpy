"""Split-complex numbers registered as a NumPy user-defined dtype."""
from ._splitcomplex1 import splitcomplex1
from .._conversions import attach as _attach_conversions

_attach_conversions(splitcomplex1, 2)

__all__ = ["splitcomplex1"]
