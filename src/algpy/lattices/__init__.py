"""Lattice types for NumPy, implemented in C.

Exterior algebras are the Clifford algebras with every generator square zero.
They expose a lattice interface: join (``|``) with identity ``bottom()``,
meet (``&``) with identity ``top()``, and the Poincaré duals ``ldual`` and
``rdual``, which coincide in odd dimensions and are then also ``~``.
The ufuncs ``ldual`` and ``rdual`` apply the duals to arrays of any lattice
dtype.
"""
from .._ufuncs import ldual, rdual
from .clifford import exterior1, exterior2, exterior3, exterior4, exterior5

__all__ = ["ldual", "rdual", "exterior1", "exterior2", "exterior3", "exterior4", "exterior5"]
