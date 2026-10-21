import numpy as _np


def attach(scalar_type, n):
    def _to_float(a, dtype=float):
        a = _np.asarray(a, dtype=scalar_type)
        flat = _np.ascontiguousarray(a.reshape(-1))
        return flat.view(float).reshape(a.shape + (n,)).astype(dtype, copy=False)

    def _from_float(a):
        a = _np.asarray(a, dtype=float)
        if a.ndim == 0 or a.shape[-1] != n:
            raise ValueError(f"last axis must have length {n}")
        flat = _np.ascontiguousarray(a.reshape(-1))
        return flat.view(scalar_type).reshape(a.shape[:-1])

    scalar_type.N = n
    scalar_type.to_float = staticmethod(_to_float)
    scalar_type.from_float = staticmethod(_from_float)
