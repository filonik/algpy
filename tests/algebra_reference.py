"""Algebra catalogue and cached general types used by shared checks."""
from functools import lru_cache
from algpy import algebras
from algpy.algebras import general

FAMILIES = ("cayleydickson", "clifford", "multicomplex")

FIXED_METRICS = {
    "pga1": (0, 1),
    "pga2": (0, 1, 1),
    "pga3": (0, 1, 1, 1),
    "pga4": (0, 1, 1, 1, 1),
    "cga1": (-1, 1, 1),
    "cga2": (-1, 1, 1, 1),
    "cga3": (-1, 1, 1, 1, 1),
}

COMPLEX_NORMS = {"pga1", "cga1", "cga2", *(f"splitcomplex{n}" for n in range(1, 6))}

def expected_norm_type(algebra):
    return complex if hasattr(algebra, "K") or algebra.__name__ in COMPLEX_NORMS else float

CONCRETE = (
    *((name, 2**len(metric)) for name, metric in FIXED_METRICS.items()),
    ("real", 1), ("complex1", 2), ("quaternion", 4),
    ("octonion", 8), ("sedenion", 16), ("trigintaduonion", 32),
    *((f"{family}{n}", 2**n)
      for family in ("complex", "dual", "splitcomplex") for n in range(1, 6)
      if (family, n) != ("complex", 1)),
    # Jets R[e]/(e^(n+1)): one generator, n+1 coefficients.
    *((f"jet{n}", n + 1) for n in range(1, 6)),
    # Complex jets C[e]/(e^(n+1)): generators i and e, 2(n+1) real components.
    *((f"complex_jet{n}", 2*(n + 1)) for n in range(1, 6)),
)

UNARY = (
    ("positive", None), ("negative", None), ("conjugate", "conjugate"),
    ("reciprocal", "inv"), ("square", "square"),
    ("exp", "exp"), ("expm1", "expm1"), ("log", "log"),
    ("log1p", "log1p"), ("sqrt", "sqrt"),
    ("sin", "sin"), ("cos", "cos"), ("tan", "tan"),
    ("arcsin", "asin"), ("arccos", "acos"), ("arctan", "atan"),
    ("sinh", "sinh"), ("cosh", "cosh"), ("tanh", "tanh"),
    ("arcsinh", "asinh"), ("arccosh", "acosh"), ("arctanh", "atanh"),
)

@lru_cache(None)
def algebra_type(name):
    if any(name.startswith(family) for family in FAMILIES):
        return getattr(general, name)(*([-1]*int(name[-1])))
    return getattr(algebras, name)
