"""Field function catalogue and independent Cauchy derivative references."""
import math
import numpy as np

FUNCTIONS = [
    ("exp", np.exp), ("expm1", np.expm1), ("log", np.log), ("log1p", np.log1p),
    ("inv", np.reciprocal), ("square", np.square), ("sqrt", np.sqrt),
    ("sin", np.sin), ("cos", np.cos), ("tan", np.tan),
    ("asin", np.arcsin), ("acos", np.arccos), ("atan", np.arctan),
    ("sinh", np.sinh), ("cosh", np.cosh), ("tanh", np.tanh),
    ("asinh", np.arcsinh), ("acosh", np.arccosh), ("atanh", np.arctanh),
]

INDICES = {name: i for i, (name, _) in enumerate(FUNCTIONS)}

def reference(function, value, order):
    # Extract a Taylor coefficient with a discrete contour integral. This
    # checks the derivatives using only the function, without their formulas
    # or the cancellation of a fourth-order finite difference.
    radius = 0.125
    circle = np.exp(2j*np.pi*np.arange(256)/256)
    return (math.factorial(order)/radius**order
            * np.mean(function(value+radius*circle)*circle**(-order)))

def operation(name):
    return {"inv": "reciprocal", "asin": "arcsin", "acos": "arccos",
            "atan": "arctan", "asinh": "arcsinh", "acosh": "arccosh",
            "atanh": "arctanh"}.get(name, name)

def center(name):
    return 2 if name == "acosh" else 1.5 if name in {"log", "sqrt", "inv"} else 0.3
