"""Dual numbers registered as a NumPy user-defined dtype."""
from ._dual1 import dual1
from .._conversions import attach as _attach_conversions

_attach_conversions(dual1, 2)

__all__ = ["dual1"]
