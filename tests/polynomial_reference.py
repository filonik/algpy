"""Root matching independent of ordering, preserving multiplicity."""
import numpy as np
import itertools

def assert_roots(actual, expected, rtol=2e-12, atol=2e-12):
    assert np.isfinite(actual).all(), actual
    # Order is deliberately unspecified; preserve multiplicities when matching.
    assert any(np.all(np.abs(actual-np.array(order)) <= atol+rtol*np.abs(order))
               for order in itertools.permutations(expected)), (actual, expected)
