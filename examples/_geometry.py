"""Constructors and decoders for geometric objects in the example notebooks.

Everything else comes from the algebra types: the geometric product `*`,
grade projections `x.grade[k]`, the exterior and regressive products `|` and
`&`, and the Poincaré duals `ldual()`/`rdual()`. For vectors, the inner
product is the scalar part of the geometric product, `(a*b)[0]`.
"""
from functools import lru_cache

import numpy as np


# Projective geometric algebra, dual model: n is the null generator.

def pga_point(algebra, coordinates):
    """Finite point: the left dual of the homogeneous vector n + x."""
    generators = algebra.gens()
    if len(coordinates) != len(generators) - 1:
        raise ValueError("Coordinate dimension does not match this PGA")
    homogeneous = generators[0]
    for coordinate, generator in zip(coordinates, generators[1:]):
        homogeneous = homogeneous + float(coordinate)*generator
    return homogeneous.ldual()


def pga_hyperplane(algebra, normal, offset):
    """Line (PGA2) or plane (PGA3) normal·x = offset: -offset*n + normal."""
    generators = algebra.gens()
    if len(normal) != len(generators) - 1:
        raise ValueError("Normal dimension does not match this PGA")
    result = -float(offset)*generators[0]
    for coordinate, generator in zip(normal, generators[1:]):
        result = result + float(coordinate)*generator
    return result


def pga_coordinates(point):
    """Euclidean coordinates of a finite point."""
    algebra = type(point)
    indices = 1 << np.arange(algebra.N.bit_length() - 1)
    values = algebra.to_float(point.rdual())[indices]
    if abs(values[0]) < 1e-12:
        raise ValueError("This example expects a finite point")
    return values[1:]/values[0]


# Conformal geometric algebra: generators e-, e+, then Euclidean axes.

@lru_cache(None)
def cga_frame(algebra):
    """The origin and infinity points, and the Euclidean axes."""
    negative, positive, *euclidean = algebra.gens()
    return (negative - positive)/2, negative + positive, tuple(euclidean)


def cga_point(algebra, coordinates):
    """Null vector origin + x + |x|^2*infinity/2."""
    origin, infinity, euclidean = cga_frame(algebra)
    if len(coordinates) != len(euclidean):
        raise ValueError("Coordinate dimension does not match this CGA")
    result = origin + float(np.dot(coordinates, coordinates))/2*infinity
    for coordinate, generator in zip(coordinates, euclidean):
        result = result + float(coordinate)*generator
    return result


def cga_plane(algebra, normal, offset):
    """Outer-product plane normal·x = offset, the undual of normal + offset*infinity."""
    _, infinity, euclidean = cga_frame(algebra)
    if len(normal) != len(euclidean):
        raise ValueError("Normal dimension does not match this CGA")
    vector = float(offset)*infinity
    for coordinate, axis in zip(normal, euclidean):
        vector = vector + float(coordinate)*axis
    return vector*algebra.top()


def cga_coordinates(point):
    """Euclidean coordinates of a finite conformal point."""
    _, infinity, euclidean = cga_frame(type(point))
    weight = -(point*infinity)[0]
    if abs(weight) < 1e-12:
        raise ValueError("This example expects a finite conformal point")
    return np.array([(point*axis)[0]/weight for axis in euclidean])


def _round_data(sphere):
    """Center and radius of an inner-product circle/sphere vector."""
    _, infinity, euclidean = cga_frame(type(sphere))
    normalized = sphere/(-(sphere*infinity)[0])
    radius2 = (normalized*normalized)[0]
    if radius2 < -1e-10:
        raise ValueError("This example expects a real circle or sphere")
    center = np.array([(normalized*axis)[0] for axis in euclidean])
    return center, np.sqrt(max(0, radius2))


def cga_sphere_data(sphere):
    """Center and radius of an outer-product circle (CGA2) or sphere (CGA3)."""
    return _round_data(sphere*type(sphere).top().inv())


def cga_point_pair(pair):
    """The two null points of an outer-product point pair.

    If B=P^Q, D=infinity contraction B, then (B +/- sqrt(B^2))*D^{-1}
    are proportional to the two null points. Tangencies need a separate limit.
    """
    _, infinity, _ = cga_frame(type(pair))
    discriminant = (pair*pair)[0]
    if discriminant <= 1e-12:
        raise ValueError("Choose two distinct real intersections for this demo")
    direction = (infinity*pair).grade[1]
    root = np.sqrt(discriminant)
    return tuple(((pair + sign*root)*direction.inv()).grade[1] for sign in [-1, 1])


def cga_circle_data(circle):
    """Center, radius, and unit plane normal of an outer-product CGA3 circle."""
    _, infinity, euclidean = cga_frame(type(circle))
    dual = circle*type(circle).top().inv()
    plane = -(infinity*dual).grade[1]
    centered_sphere = -(plane*dual).grade[1]/(plane*plane)[0]
    center, radius = _round_data(centered_sphere)
    normal = np.array([(plane*axis)[0] for axis in euclidean])
    return center, radius, normal/np.linalg.norm(normal)
