"""Concrete cayleydickson algebras registered as NumPy dtypes."""
from ._real import real
from ._quaternion import quaternion
from ._octonion import octonion
from ._sedenion import sedenion

__all__ = [
    "real",
    "quaternion",
    "octonion",
    "sedenion",
]
