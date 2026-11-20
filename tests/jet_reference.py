"""Independent truncated-polynomial references for the jet algebras."""
import math

import numpy as np

from dual_reference import derivative


def product(left, right):
    """Truncated Cauchy product of coefficient vectors in R[e]/(e^n)."""
    result = np.zeros(len(left))
    for i, x in enumerate(left):
        for j, y in enumerate(right[:len(left) - i]):
            result[i + j] += x*y
    return result


def taylor(function, coefficients):
    """f(a) for a = x + n: the sum of f^(k)(x)/k! * n^k."""
    nilpotent = np.array(coefficients, float)
    nilpotent[0] = 0
    power = np.zeros(len(coefficients))
    power[0] = 1
    result = np.zeros(len(coefficients))
    for k in range(len(coefficients)):
        result += derivative(function, coefficients[0], k)/math.factorial(k)*power
        power = product(power, nilpotent)
    return result


def complex_derivative(function, z, order, radius=0.125):
    """k-th complex derivative at z by a discrete contour integral."""
    if order == 0:
        return complex(function(z))
    circle = np.exp(2j*np.pi*np.arange(1024)/1024)
    return math.factorial(order)/radius**order*np.mean(function(z + radius*circle)*circle**(-order))


def complex_product(left, right):
    result = np.zeros(len(left), complex)
    for i, x in enumerate(left):
        for j, y in enumerate(right[:len(left) - i]):
            result[i + j] += x*y
    return result


def complex_taylor(function, coefficients):
    """f(a) for complex coefficients a = z + n: the sum of f^(k)(z)/k! * n^k."""
    nilpotent = np.array(coefficients, complex)
    nilpotent[0] = 0
    power = np.zeros(len(coefficients), complex)
    power[0] = 1
    result = np.zeros(len(coefficients), complex)
    for k in range(len(coefficients)):
        result += complex_derivative(function, coefficients[0], k)/math.factorial(k)*power
        power = complex_product(power, nilpotent)
    return result


def to_complex(cls, values):
    """Complex coefficients of complex jets: real components 2k, 2k+1."""
    data = cls.to_float(values)
    return data[..., 0::2] + 1j*data[..., 1::2]


def from_complex(cls, coefficients):
    coefficients = np.asarray(coefficients, complex)
    data = np.empty(coefficients.shape[:-1] + (2*coefficients.shape[-1],))
    data[..., 0::2], data[..., 1::2] = coefficients.real, coefficients.imag
    return cls.from_float(data)
