#ifndef ALGPY_NORM_H
#define ALGPY_NORM_H

/* Python conversion of a norm: a float if real, an algpy complex1 otherwise. */
#include <Python.h>

#include "cayleydickson/complex.h"

static PyObject *
algpy_norm_to_py(complex n)
{
    PyObject *module, *type, *result;
    if (n.data[1] == 0) {
        return PyFloat_FromDouble(n.data[0]);
    }
    module = PyImport_ImportModule("algpy.algebras.multicomplex");
    if (module == NULL) {
        return NULL;
    }
    type = PyObject_GetAttrString(module, "complex1");
    Py_DECREF(module);
    if (type == NULL) {
        return NULL;
    }
    result = PyObject_CallFunction(type, "dd", n.data[0], n.data[1]);
    Py_DECREF(type);
    return result;
}

#endif /* ALGPY_NORM_H */
