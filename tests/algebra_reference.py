"""Algebra catalogue and cached general types used by shared checks."""
from functools import lru_cache
from algpy import algebras
from algpy.algebras import general

FAMILIES = ("cayleydickson", "clifford", "multicomplex")

CONCRETE = (
    ("real", 1), ("complex1", 2), ("quaternion", 4),
    ("octonion", 8), ("sedenion", 16),
    *((f"{family}{n}", 2**n)
      for family in ("complex", "dual", "splitcomplex") for n in range(1, 5)
      if (family, n) != ("complex", 1)),
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
