"""Shared assertions for exterior lattices: join, meet, duals and bindings."""
import pickle
import numpy as np
import pytest
from algpy.algebras.general import clifford
import lattice_reference as ref


def components(cls, value):
    return cls.to_float(value)


def assert_components(cls, value, expected):
    np.testing.assert_allclose(components(cls, value), expected, rtol=1e-12, atol=1e-12)


def check_constants(cls):
    n = ref.generators(cls)
    assert cls.N == 2**n
    np.testing.assert_array_equal(cls.to_float(cls.basis()), np.eye(cls.N))
    gens = cls.gens()
    assert gens.shape == (n,) and gens.dtype == cls.dtype
    np.testing.assert_array_equal(cls.to_float(gens), np.eye(cls.N)[[1 << i for i in range(n)]])
    assert_components(cls, cls.zero(), np.zeros(cls.N))
    assert_components(cls, cls.bottom(), np.eye(cls.N)[0])
    assert_components(cls, cls.top(), np.eye(cls.N)[-1])


def check_join_and_meet_reference(cls):
    for x, y in zip(ref.samples(cls), ref.samples(cls, seed=1)):
        a, b = cls(*x), cls(*y)
        assert_components(cls, a | b, ref.join(x, y))
        assert_components(cls, a & b, ref.meet(x, y))


def check_join_is_degenerate_clifford_product(cls):
    """The exterior algebra is the Clifford algebra with all squares zero."""
    algebra = clifford(*([0]*ref.generators(cls)))
    for x, y in zip(ref.samples(cls), ref.samples(cls, seed=1)):
        np.testing.assert_allclose(
            cls.to_float(cls(*x) | cls(*y)),
            algebra.to_float(algebra(*x)*algebra(*y)), rtol=1e-12, atol=1e-12)


def check_identities_and_absorbing_elements(cls):
    for x in ref.samples(cls):
        a = cls(*x)
        assert a | cls.bottom() == a and cls.bottom() | a == a
        assert a & cls.top() == a and cls.top() & a == a
        assert a | cls.zero() == cls.zero()
        assert a & cls.zero() == cls.zero()


def check_associativity(cls):
    x, y, z = (cls(*v) for v in ref.samples(cls))
    assert_components(cls, (x | y) | z, components(cls, x | (y | z)))
    assert_components(cls, (x & y) & z, components(cls, x & (y & z)))


def check_generators(cls):
    gens = list(cls.gens())
    top = cls.bottom()
    for g in gens:
        assert g | g == cls.zero()
        top = top | g
    assert top == cls.top()
    for i, g in enumerate(gens):
        for h in gens[i+1:]:
            assert g | h == -(h | g)


def check_duals_on_basis(cls):
    for blade in cls.basis():
        assert blade.rdual() | blade == cls.top()
        assert blade | blade.ldual() == cls.top()
        assert blade.ldual().rdual() == blade
        assert blade.rdual().ldual() == blade


def check_duals_reference_and_inverse(cls):
    for x in ref.samples(cls):
        a = cls(*x)
        assert_components(cls, a.rdual(), ref.rdual(x))
        assert_components(cls, a.ldual(), ref.ldual(x))
        assert a.ldual().rdual() == a and a.rdual().ldual() == a
        assert cls.bottom().rdual() == cls.top()
        assert cls.bottom().ldual() == cls.top()
        assert cls.top().rdual() == cls.bottom()


def check_poincare_duality(cls):
    """Meet is join conjugated by either dual."""
    for x, y in zip(ref.samples(cls), ref.samples(cls, seed=1)):
        a, b = cls(*x), cls(*y)
        assert_components(cls, a & b, components(cls, (a.ldual() | b.ldual()).rdual()))
        assert_components(cls, a & b, components(cls, (a.rdual() | b.rdual()).ldual()))
        assert_components(cls, (a | b).ldual(), components(cls, a.ldual() & b.ldual()))


def check_odd_dual(cls):
    for x in ref.samples(cls):
        a = cls(*x)
        assert a.ldual() == a.rdual() == a.dual() == ~a
        assert ~~a == a
    values = cls.from_float(ref.samples(cls))
    np.testing.assert_array_equal(cls.to_float(~values),
                                  cls.to_float(np.invert(values)))
    np.testing.assert_array_equal(cls.to_float(~values),
                                  [components(cls, v.dual()) for v in values])


def check_even_duals_differ(cls):
    assert not hasattr(cls, "dual")
    gens = cls.gens()
    assert gens[0].ldual() == -gens[0].rdual()
    with pytest.raises(TypeError, match="ldual"):
        ~cls.bottom()
    with pytest.raises(TypeError):
        np.invert(cls.from_float(ref.samples(cls)))


def check_vector_space_operations(cls):
    x, y = ref.samples(cls, 2)
    a, b = cls(*x), cls(*y)
    assert_components(cls, a + b, x + y)
    assert_components(cls, a - b, x - y)
    assert_components(cls, -a, -x)
    assert_components(cls, +a, x)
    assert_components(cls, 2*a, 2*x)
    assert_components(cls, a*np.float32(2), 2*x)
    assert_components(cls, a/4, x/4)
    assert a + 1 == a + cls.bottom()
    # Reals embed as multiples of bottom, so joining with one scales.
    assert a | 3 == 3*a and 3 | a == a*3


def check_no_lattice_products(cls):
    """Exterior lattices only scale; lattice values combine by join and meet."""
    a, b = (cls(*x) for x in ref.samples(cls, 2))
    with pytest.raises(TypeError):
        a*b
    with pytest.raises(TypeError):
        a/b
    with pytest.raises(TypeError):
        1/a
    values = cls.from_float(ref.samples(cls))
    with pytest.raises(TypeError):
        np.multiply(values, values)


def check_numpy_ufuncs(cls):
    xs, ys = ref.samples(cls, 4), ref.samples(cls, 4, seed=1)
    a, b = cls.from_float(xs), cls.from_float(ys)
    np.testing.assert_allclose(cls.to_float(a | b), [ref.join(x, y) for x, y in zip(xs, ys)], atol=1e-12)
    np.testing.assert_allclose(cls.to_float(np.bitwise_and(a, b)), [ref.meet(x, y) for x, y in zip(xs, ys)], atol=1e-12)
    np.testing.assert_allclose(cls.to_float(a | b[0]), [ref.join(x, ys[0]) for x in xs], atol=1e-12)
    np.testing.assert_allclose(cls.to_float(a[:, None] & b[None, :]),
                               [[ref.meet(x, y) for y in ys] for x in xs], atol=1e-12)
    np.testing.assert_allclose(cls.to_float(a + b), xs + ys)
    np.testing.assert_allclose(cls.to_float(-(a - b)), ys - xs)
    np.testing.assert_allclose(cls.to_float(a*np.arange(4)), xs*np.arange(4)[:, None])
    np.testing.assert_allclose(cls.to_float(np.arange(4)*a), xs*np.arange(4)[:, None])
    np.testing.assert_allclose(cls.to_float(a/2), xs/2)
    np.testing.assert_array_equal(a == a, [True]*4)
    np.testing.assert_array_equal(a != b, [True]*4)
    out = np.empty_like(a)
    assert np.bitwise_or(a, b, out=out) is out
    np.testing.assert_array_equal(cls.to_float(out), cls.to_float(a | b))
    strided = cls.from_float(np.repeat(xs, 2, axis=0))[::2]
    np.testing.assert_array_equal(cls.to_float(strided & b), cls.to_float(a & b))


def check_scalar_protocols(cls):
    x = ref.samples(cls)[0]
    a = cls(*x)
    assert len(a) == cls.N and a[0] == x[0] and a[-1] == x[-1]
    with pytest.raises(IndexError):
        a[cls.N]
    assert pickle.loads(pickle.dumps(a)) == a
    assert hash(cls(2)) == hash(2.0) and cls(2) == 2
    assert hash(a) == hash(cls(*x))
    assert repr(cls.top()).startswith(f"{cls.__name__}(0, ")
    assert bool(a) and not bool(cls.zero())
    assert a.norm2() == pytest.approx(float(x @ x))
    assert len(a.grade) == ref.generators(cls)+1
    total = cls.zero()
    for k in range(ref.generators(cls)+1):
        assert_components(cls, a.grade[k], ref.grade(x, k))
        total = total + a.grade[k]
    assert total == a
    with pytest.raises(IndexError):
        a.grade[ref.generators(cls)+1]
    with pytest.raises(TypeError):
        a.grade[0] = a
    keywords = cls(x=2, **{"x" + "".join(map(str, range(ref.generators(cls)))): 3})
    assert keywords == 2*cls.bottom() + 3*cls.top()


def check_dual_ufuncs(cls):
    from algpy.lattices import ldual, rdual
    assert isinstance(ldual, np.ufunc) and isinstance(rdual, np.ufunc)
    xs = ref.samples(cls, 4)
    values = cls.from_float(xs)
    np.testing.assert_array_equal(cls.to_float(ldual(values)), [ref.ldual(x) for x in xs])
    np.testing.assert_array_equal(cls.to_float(rdual(values)), [ref.rdual(x) for x in xs])
    np.testing.assert_array_equal(cls.to_float(ldual(values)),
                                  [components(cls, v.ldual()) for v in values])
    np.testing.assert_array_equal(cls.to_float(rdual(ldual(values))), xs)
    grid = cls.from_float(np.repeat(xs, 2, axis=0).reshape(2, 4, cls.N))[:, ::2]
    np.testing.assert_array_equal(cls.to_float(ldual(grid)),
                                  [[ref.ldual(x) for x in xs[:2]], [ref.ldual(x) for x in xs[2:]]])
    out = np.empty_like(values)
    assert rdual(values, out=out) is out
    np.testing.assert_array_equal(cls.to_float(out), cls.to_float(rdual(values)))
    scalar = ldual(values[0])
    assert isinstance(scalar, cls) and scalar == values[0].ldual()
    assert rdual(cls.bottom()) == cls.top()
    with pytest.raises(TypeError):
        ldual(np.arange(3.))


def check_matches_exterior(cls, exterior):
    """Join, meet and duals do not depend on the metric."""
    assert cls.N == exterior.N
    assert cls.bottom() == cls.one()
    np.testing.assert_array_equal(cls.to_float(cls.top()), exterior.to_float(exterior.top()))
    for x, y in zip(ref.samples(cls), ref.samples(cls, seed=1)):
        a, b = cls(*x), cls(*y)
        c, d = exterior(*x), exterior(*y)
        np.testing.assert_array_equal(cls.to_float(a | b), exterior.to_float(c | d))
        np.testing.assert_array_equal(cls.to_float(a & b), exterior.to_float(c & d))
        np.testing.assert_array_equal(cls.to_float(a.ldual()), exterior.to_float(c.ldual()))
        np.testing.assert_array_equal(cls.to_float(a.rdual()), exterior.to_float(c.rdual()))


def pga_point(cls, coordinates):
    """Finite point: the left dual of the homogeneous vector (1, x, ...)."""
    gens = cls.gens()
    homogeneous = gens[0]
    for coordinate, generator in zip(coordinates, gens[1:]):
        homogeneous = homogeneous + float(coordinate)*generator
    return homogeneous.ldual()


def pga_coordinates(point):
    cls = type(point)
    homogeneous = cls.to_float(point.rdual())[[1 << i for i in range(ref.generators(cls))]]
    return homogeneous[1:]/homogeneous[0]


def check_pga_incidence(cls):
    """Meets of points give hyperplanes; joins of hyperplanes give points."""
    dimension = ref.generators(cls) - 1
    rng = np.random.default_rng(cls.N)
    target = rng.normal(size=dimension)
    point = pga_point(cls, target)
    np.testing.assert_allclose(pga_coordinates(point), target, rtol=1e-12)
    hyperplanes = []
    for _ in range(dimension):
        others = [pga_point(cls, rng.normal(size=dimension)) for _ in range(dimension - 1)]
        hyperplane = point
        for other in others:
            hyperplane = hyperplane & other
        assert cls.to_float(hyperplane)[[1 << i for i in range(dimension + 1)]].any()
        for incident in [point, *others]:
            np.testing.assert_allclose(cls.to_float(hyperplane | incident), 0, atol=1e-12)
        hyperplanes.append(hyperplane)
    intersection = cls.bottom()
    for hyperplane in hyperplanes:
        intersection = intersection | hyperplane
    np.testing.assert_allclose(pga_coordinates(intersection), target, rtol=1e-9, atol=1e-9)
    elsewhere = pga_point(cls, target + 1)
    assert abs(cls.to_float(hyperplanes[0] | elsewhere)[-1]) > 1e-6


def cga_point(cls, coordinates):
    """Null vector origin + x + |x|^2 infinity/2, with the -1, +1 generators first."""
    negative, positive, *euclidean = cls.gens()
    origin, infinity = (negative - positive)/2, negative + positive
    result = origin + float(np.dot(coordinates, coordinates))/2*infinity
    for coordinate, generator in zip(coordinates, euclidean):
        result = result + float(coordinate)*generator
    return result


def check_cga_incidence_and_metric_meet(cls):
    """Joins of points are round objects through them; meets agree with the
    metric meet undual(dual(a) | dual(b)), with dual(a) = a*I^-1, up to sign."""
    dimension = ref.generators(cls) - 2
    rng = np.random.default_rng(cls.N)
    points = [cga_point(cls, rng.normal(size=dimension)) for _ in range(dimension + 1)]
    for point in points:
        assert abs(cls.to_float(point*point)).max() < 1e-12
    sphere = points[0]
    for point in points[1:]:
        sphere = sphere | point
    for point in points:
        np.testing.assert_allclose(cls.to_float(sphere | point), 0, atol=1e-12)
    elsewhere = cga_point(cls, rng.normal(size=dimension))
    assert abs(cls.to_float(sphere | elsewhere)).max() > 1e-6

    pseudoscalar = cls.top()
    inverse = pseudoscalar.inv()
    vectors = [cls(*np.where([bin(i).count("1") == 1 for i in range(cls.N)],
                             rng.normal(size=cls.N), 0)) for _ in range(2*ref.generators(cls))]
    for ga in range(1, ref.generators(cls)):
        for gb in range(ref.generators(cls) - ga + 1, ref.generators(cls)):
            a = b = cls.bottom()
            for v in vectors[:ga]:
                a = a | v
            for v in vectors[ga:ga + gb]:
                b = b | v
            metric = ((a*inverse) | (b*inverse))*pseudoscalar
            lattice = cls.to_float(a & b)
            assert abs(lattice).max() > 1e-9
            sign = np.sign(lattice @ cls.to_float(metric))
            np.testing.assert_allclose(cls.to_float(metric), sign*lattice, rtol=1e-9, atol=1e-9)
