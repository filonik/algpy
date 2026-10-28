#ifndef ALGPY_CONVERSIONS_H
#define ALGPY_CONVERSIONS_H

#include <Python.h>
#include <numpy/arrayobject.h>

/* Wrapper helpers. Each extension imports its own NumPy C API before these
 * functions are called. Algebra storage is N consecutive native doubles. */
static inline int
algpy_set_dimension(PyObject *type, int n)
{
    PyObject *size = PyLong_FromLong(n);
    int result;
    if (size == NULL) {
        return -1;
    }
    result = PyObject_SetAttrString(type, "N", size);
    Py_DECREF(size);
    return result;
}

static inline PyObject *
algpy_to_float(PyObject *args, PyObject *kwargs, int typenum, int n)
{
    static char *keywords[] = {"a", "dtype", NULL};
    PyObject *input, *dtype = Py_None;
    PyArrayObject *array = NULL, *flat = NULL, *view = NULL, *reshaped = NULL;
    PyArray_Descr *descr;
    PyObject *result = NULL;
    PyObject *astype = NULL, *cast_args = NULL, *cast_kwargs = NULL;
    npy_intp dimensions[NPY_MAXDIMS+1];
    PyArray_Dims shape;
    int i;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|O:to_float", keywords, &input, &dtype)) {
        return NULL;
    }
    descr = PyArray_DescrFromType(typenum);
    if (descr == NULL) {
        return NULL;
    }
    /* FromAny and View steal their descriptor arguments. Ravel uses C order,
     * sharing contiguous data and copying strided data as needed. */
    array = (PyArrayObject *)PyArray_FromAny(input, descr, 0, 0,
                                           NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST, NULL);
    if (array == NULL) {
        goto done;
    }
    flat = (PyArrayObject *)PyArray_Ravel(array, NPY_CORDER);
    if (flat == NULL) {
        goto done;
    }
    descr = PyArray_DescrFromType(NPY_DOUBLE);
    if (descr == NULL) {
        goto done;
    }
    view = (PyArrayObject *)PyArray_View(flat, descr, NULL);
    if (view == NULL) {
        goto done;
    }
    for (i = 0; i < PyArray_NDIM(array); i++) {
        dimensions[i] = PyArray_DIM(array, i);
    }
    dimensions[i] = n;
    shape.ptr = dimensions;
    shape.len = PyArray_NDIM(array)+1;
    reshaped = (PyArrayObject *)PyArray_Newshape(view, &shape, NPY_CORDER);
    if (reshaped == NULL) {
        goto done;
    }
    if (dtype == Py_None) {
        result = (PyObject *)reshaped;
        reshaped = NULL;
        goto done;
    }
    /* Use NumPy's astype for its dtype inference (including unsized strings
     * and subarray dtypes), with the same copy=False behavior as the API. */
    astype = PyObject_GetAttrString((PyObject *)reshaped, "astype");
    if (astype == NULL) {
        goto done;
    }
    cast_args = PyTuple_Pack(1, dtype);
    if (cast_args == NULL) {
        goto done;
    }
    cast_kwargs = Py_BuildValue("{s:O}", "copy", Py_False);
    if (cast_kwargs == NULL) {
        goto done;
    }
    result = PyObject_Call(astype, cast_args, cast_kwargs);
done:
    Py_XDECREF(cast_kwargs);
    Py_XDECREF(cast_args);
    Py_XDECREF(astype);
    Py_XDECREF(reshaped);
    Py_XDECREF(view);
    Py_XDECREF(flat);
    Py_XDECREF(array);
    return result;
}

static inline PyObject *
algpy_from_float(PyObject *args, PyObject *kwargs, int typenum, int n)
{
    static char *keywords[] = {"a", NULL};
    PyObject *input, *result = NULL;
    PyArrayObject *array = NULL, *flat = NULL, *view = NULL;
    PyArray_Descr *descr;
    npy_intp dimensions[NPY_MAXDIMS];
    PyArray_Dims shape;
    int i;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O:from_float", keywords, &input)) {
        return NULL;
    }
    descr = PyArray_DescrFromType(NPY_DOUBLE);
    if (descr == NULL) {
        return NULL;
    }
    array = (PyArrayObject *)PyArray_FromAny(input, descr, 0, 0,
                                           NPY_ARRAY_ENSUREARRAY | NPY_ARRAY_FORCECAST, NULL);
    if (array == NULL) {
        goto done;
    }
    if (PyArray_NDIM(array) == 0 || PyArray_DIM(array, PyArray_NDIM(array)-1) != n) {
        PyErr_Format(PyExc_ValueError, "last axis must have length %d", n);
        goto done;
    }
    flat = (PyArrayObject *)PyArray_Ravel(array, NPY_CORDER);
    if (flat == NULL) {
        goto done;
    }
    descr = PyArray_DescrFromType(typenum);
    if (descr == NULL) {
        goto done;
    }
    view = (PyArrayObject *)PyArray_View(flat, descr, NULL);
    if (view == NULL) {
        goto done;
    }
    for (i = 0; i < PyArray_NDIM(array)-1; i++) {
        dimensions[i] = PyArray_DIM(array, i);
    }
    shape.ptr = dimensions;
    shape.len = PyArray_NDIM(array)-1;
    result = PyArray_Newshape(view, &shape, NPY_CORDER);
done:
    Py_XDECREF(view);
    Py_XDECREF(flat);
    Py_XDECREF(array);
    return result;
}

#endif /* ALGPY_CONVERSIONS_H */
