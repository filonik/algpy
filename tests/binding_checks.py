"""Shared regression checks for the scalar and dtype wrapper boundary."""
import operator

import numpy as np


def check_nan_hash_stable(algebra):
    for component in range(algebra.N):
        coefficients = [0.]*algebra.N
        coefficients[component] = float("nan")
        value = algebra(*coefficients)
        expected_hash = hash(value)
        marker = object()
        mapping = {value: marker}
        # Keep new float allocations alive so hashing cannot accidentally
        # reuse the same temporary PyFloat address for a NaN component.
        occupied = [i + .5 for i in range(512)]
        assert hash(value) == expected_hash
        assert mapping[value] is marker
        assert len({value, value}) == 1
        del occupied


def check_subclass_mixed_arithmetic(algebra):
    class Derived(algebra):
        def __new__(cls, *components):
            return super().__new__(cls, *components)

    components = [2.] + [.125]*(algebra.N - 1)
    value = algebra(*components)
    derived = Derived(*components)
    for function in [operator.add, operator.sub, operator.mul, operator.truediv, operator.pow]:
        np.testing.assert_allclose(algebra.to_float(function(derived, 1)),
                                   algebra.to_float(function(value, 1)), rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(algebra.to_float(function(1, derived)),
                                   algebra.to_float(function(1, value)), rtol=2e-14, atol=2e-14)
    assert (derived == 2) == (value == 2)
    assert (derived != 2) == (value != 2)
    assert (2 == derived) == (2 == value)


def check_dtype_byteorder(algebra):
    rows = np.arange(1., 2*algebra.N + 1).reshape(2, algebra.N)
    native = algebra.from_float(rows)
    swapped = native.byteswap().view(native.dtype.newbyteorder())
    np.testing.assert_array_equal(algebra.to_float(swapped), rows)
    for value, expected in zip(swapped, rows):
        np.testing.assert_array_equal(algebra.to_float(value), expected)

    replacement = -rows[0] - .25
    swapped[::-1][0] = algebra(*replacement)
    np.testing.assert_array_equal(algebra.to_float(swapped)[1], replacement)
    np.testing.assert_array_equal(algebra.to_float(swapped[1]), replacement)
    swapped[0] = tuple(replacement)
    np.testing.assert_array_equal(algebra.to_float(swapped)[0], replacement)

    negative_zeros = algebra.from_float(-np.zeros((2, algebra.N)))
    swapped_zeros = negative_zeros.byteswap().view(negative_zeros.dtype.newbyteorder())
    assert np.count_nonzero(swapped_zeros) == 0
    assert not bool(swapped_zeros[0])


def check_real_input_dtypes(algebra):
    for dtype in [np.float16, np.longdouble]:
        source = np.array([-1.25, 0, 2.5], dtype=dtype)[::-1]
        expected = np.zeros((len(source), algebra.N))
        expected[:, 0] = source.astype(float)
        np.testing.assert_array_equal(algebra.to_float(source.astype(algebra)), expected)
        np.testing.assert_array_equal(algebra.to_float(source), expected)
        swapped = source.byteswap().view(source.dtype.newbyteorder())
        np.testing.assert_array_equal(algebra.to_float(swapped), expected)

        rows = np.full((len(source), algebra.N), .125)
        rows[:, 0] = .5
        values = algebra.from_float(rows)
        added = rows.copy()
        added[:, 0] += source
        multiplied = rows*source.astype(float)[:, None]
        for left, right in [(values, source), (source, values)]:
            np.testing.assert_array_equal(algebra.to_float(np.add(left, right)), added)
            np.testing.assert_array_equal(algebra.to_float(np.multiply(left, right)), multiplied)

    # Exercise binary16 storage encodings directly, including normal and
    # subnormal boundaries, signed zero, infinities and signed NaNs.
    bits = np.array([0, 0x8000, 1, 0x3ff, 0x400, 0x3c00, 0x7bff,
                     0x7c00, 0xfc00, 0x7e00, 0xfe00], dtype=np.uint16)
    if algebra.N == 1:
        bits = np.arange(65536, dtype=np.uint16)
    source = bits.view(np.float16)
    with np.errstate(all="ignore"):
        expected = source.astype(float)
        actual = algebra.to_float(source.astype(algebra))
    np.testing.assert_array_equal(actual[:, 0], expected)
    np.testing.assert_array_equal(np.signbit(actual[:, 0]), np.signbit(expected))
    np.testing.assert_array_equal(actual[:, 1:], 0)
