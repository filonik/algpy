"""Multicomplex number types registered as numpy user-defined dtypes.

``multicomplex1(k)`` returns the scalar type (usable as a numpy dtype) of the
numbers ``x + y*e`` with ``e*e = k``; ``multicomplex2(k0, k1)`` that of the
numbers ``w + x*e1 + y*e2 + z*e1*e2`` with ``e1*e1 = k0`` and ``e2*e2 = k1``;
``multicomplex3(k0, k1, k2)`` the 8-component numbers over ``e0, e1, e2``; ``multicomplex4(k0, k1, k2, k3)`` the 16-component
numbers over ``e0, e1, e2, e3``.
"""
from ._multicomplex1 import multicomplex1
from ._multicomplex2 import multicomplex2
from ._multicomplex3 import multicomplex3
from ._multicomplex4 import multicomplex4

_CONSTRUCTORS = {1: multicomplex1, 2: multicomplex2, 3: multicomplex3, 4: multicomplex4}


def multicomplex(*args):
    """Dispatch to ``multicomplexN`` where N is the number of arguments."""
    constructor = _CONSTRUCTORS.get(len(args))
    if constructor is None:
        supported = ", ".join(str(n) for n in sorted(_CONSTRUCTORS))
        raise TypeError(
            f"multicomplex expects {supported} argument(s), got {len(args)}"
        )
    return constructor(*args)


__all__ = ["multicomplex", "multicomplex1", "multicomplex2", "multicomplex3", "multicomplex4"]
