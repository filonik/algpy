"""Quaternion number type registered as a numpy user-defined dtype.

Use ``from algpy.algebras import quaternion`` and pass ``quaternion`` as the
dtype, e.g. ``numpy.array([...], dtype=quaternion)``.
"""
from ._quaternion import quaternion
from .._conversions import attach as _attach_conversions

_attach_conversions(quaternion, 4)

__all__ = ["quaternion"]
