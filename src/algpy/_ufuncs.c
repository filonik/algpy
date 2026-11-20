/*
 * NumPy ufuncs shared by several algpy types. They are created without
 * loops; each extension registers its own loops when it is imported, so a
 * single ufunc dispatches on every dtype that supports it. Currently the
 * Poincare duals ldual and rdual, registered by the exterior lattices and
 * the PGA and CGA algebras, and re-exported by algpy.lattices.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

static int
add_ufunc(PyObject *m, const char *name, const char *doc)
{
    PyObject *ufunc = PyUFunc_FromFuncAndData(
        NULL, NULL, NULL, 0, 1, 1, PyUFunc_None, name, doc, 0);
    if (ufunc == NULL) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, name, ufunc) < 0) {
        Py_DECREF(ufunc);
        return -1;
    }
    Py_DECREF(ufunc);
    return 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_ufuncs", "Ufuncs shared by algpy types",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__ufuncs(void)
{
    PyObject *m;
    import_array();
    import_umath();
    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }
    if (add_ufunc(m, "ldual",
            "ldual(x, /, out=None, *, where=True, ...)\n\n"
            "Left Poincare dual of lattice values: x | ldual(x) == top on\n"
            "basis blades.") < 0 ||
        add_ufunc(m, "rdual",
            "rdual(x, /, out=None, *, where=True, ...)\n\n"
            "Right Poincare dual of lattice values: rdual(x) | x == top on\n"
            "basis blades.") < 0) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
