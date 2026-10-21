"""Complex1 number type registered as a NumPy user-defined dtype."""
from ._complex1 import complex1
from .._conversions import attach as _attach_conversions

_attach_conversions(complex1, 2)

__all__ = ["complex1"]
