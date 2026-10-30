"""Native coefficient conversions preserve array shapes and view semantics."""


import inspect




import numpy as np

import pytest

from algpy import algebras

from algpy.algebras import general


def check_native_methods_and_scalar_shape(algebra):
    cls = algebra
    for method in [cls.to_float, cls.from_float]:
        assert inspect.isbuiltin(method)
        assert method.__self__ is cls
    coefficients = np.arange(cls.N, dtype=float)+0.25
    scalar = cls(*coefficients)
    np.testing.assert_array_equal(cls.to_float(scalar), coefficients)
    values = cls.from_float(a=coefficients.tolist())
    assert isinstance(values, np.ndarray)
    assert values.ndim == 0 and values.dtype == cls.dtype
    np.testing.assert_array_equal(cls.to_float(a=values), coefficients)
    np.testing.assert_array_equal(cls.to_float(values.item()), coefficients)
    # Class methods retain the old ability to be called through an instance.
    np.testing.assert_array_equal(scalar.to_float(scalar), coefficients)
    numeric = cls.to_float(np.array([1, 2], dtype=np.int32))
    expected = np.zeros((2, cls.N))
    expected[:, 0] = [1, 2]
    np.testing.assert_array_equal(numeric, expected)

def check_shapes_and_contiguous_sharing(algebra):
    cls = algebra
    data = np.arange(6*cls.N, dtype=float).reshape(2, 3, cls.N)
    values = cls.from_float(data)
    assert values.shape == (2, 3) and values.dtype == cls.dtype
    coefficients = cls.to_float(values)
    assert np.shares_memory(values, data)
    assert np.shares_memory(coefficients, values)
    np.testing.assert_array_equal(coefficients, data)
    coefficients[1, 2, 0] = -5
    assert data[1, 2, 0] == values[1, 2][0] == -5
    values[0, 0] = cls(7)
    np.testing.assert_array_equal(data[0, 0], [7]+[0]*(cls.N-1))

def check_strided_fortran_and_reversed_inputs(algebra):
    cls = algebra
    storage = np.arange(12*cls.N, dtype=float).reshape(2, 3, 2*cls.N)
    for data in [storage[..., ::2], storage[..., :cls.N].transpose(1, 0, 2),
                 storage[..., :cls.N][::-1, :, ::-1],
                 np.asfortranarray(storage[..., :cls.N])]:
        values = cls.from_float(data)
        assert values.shape == data.shape[:-1]
        np.testing.assert_array_equal(cls.to_float(values), data)
        if not data.flags.c_contiguous:
            assert not np.shares_memory(values, data)
    data = np.arange(6*cls.N, dtype=float).reshape(2, 3, cls.N)
    values = cls.from_float(data)
    for strided, expected in [(values[::-1, ::2], data[::-1, ::2]),
                              (values.T, data.transpose(1, 0, 2))]:
        output = cls.to_float(strided)
        np.testing.assert_array_equal(output, expected)
        if not strided.flags.c_contiguous:
            assert not np.shares_memory(output, strided)

def check_optional_output_dtype(algebra, dtype):
    cls = algebra
    data = np.arange(2*cls.N, dtype=float).reshape(2, cls.N)+0.25
    values = cls.from_float(data)
    result = cls.to_float(values, dtype)
    expected = data.astype(dtype, copy=False)
    assert result.dtype == expected.dtype
    np.testing.assert_array_equal(result, expected)
    if expected.dtype == np.dtype(float):
        assert np.shares_memory(result, values)
    else:
        assert not np.shares_memory(result, values)
    np.testing.assert_array_equal(cls.to_float(a=values, dtype=dtype), expected)

def check_readonly_and_byte_swapped_inputs(algebra):
    cls = algebra
    data = np.arange(2*cls.N, dtype=float).reshape(2, cls.N)
    data.flags.writeable = False
    values = cls.from_float(data)
    output = cls.to_float(values)
    assert not values.flags.writeable and not output.flags.writeable
    assert np.shares_memory(values, data) and np.shares_memory(output, data)
    swapped = np.asarray(data, dtype=np.dtype(float).newbyteorder("S"))
    values = cls.from_float(swapped)
    np.testing.assert_array_equal(cls.to_float(values), data)
    assert not np.shares_memory(values, swapped)

def check_empty_inputs_and_invalid_shapes(algebra):
    cls = algebra
    for shape in [(0, cls.N), (2, 0, 3, cls.N)]:
        values = cls.from_float(np.empty(shape))
        assert values.shape == shape[:-1] and values.dtype == cls.dtype
        assert cls.to_float(values).shape == shape
    for bad in [1, np.array(1.0), np.empty((2, cls.N+1)), np.empty((2, 0))]:
        with pytest.raises(ValueError, match=f"last axis must have length {cls.N}"):
            cls.from_float(bad)
    with pytest.raises(TypeError):
        cls.to_float()
    with pytest.raises(TypeError):
        cls.from_float()
    with pytest.raises(TypeError):
        cls.to_float(cls.one(), unknown=True)
    with pytest.raises(TypeError):
        cls.from_float([0]*cls.N, dtype=float)

def check_numpy_dimension_limit():
    cls = algebras.complex1
    # Match NumPy's dimension limit without writing past the shape buffer.
    max_dimensions = 64 if np.lib.NumpyVersion(np.__version__) >= "2.0.0" else 32
    values = np.array([cls.one()], dtype=cls).reshape((1,)*max_dimensions)
    with pytest.raises(ValueError):
        cls.to_float(values)
