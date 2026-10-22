"""Multicomplex specializations (complex, dual, split-complex) registered as NumPy user-defined dtypes."""
from ._complex1 import complex1
from ._complex2 import complex2
from ._dual1 import dual1
from ._dual2 import dual2
from ._splitcomplex1 import splitcomplex1
from ._splitcomplex2 import splitcomplex2
from .._conversions import attach as _attach_conversions

_attach_conversions(complex1, 2)
_attach_conversions(complex2, 4)
_attach_conversions(dual1, 2)
_attach_conversions(dual2, 4)
_attach_conversions(splitcomplex1, 2)
_attach_conversions(splitcomplex2, 4)

__all__ = ["complex1", "complex2", "dual1", "dual2", "splitcomplex1", "splitcomplex2"]
