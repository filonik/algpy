#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *
hello_world(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyUnicode_FromString("Hello, world from algpy!");
}

static PyMethodDef core_methods[] = {
    {"hello_world", hello_world, METH_NOARGS,
     "Return a greeting from C."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef core_module = {
    PyModuleDef_HEAD_INIT, "_core", NULL, -1, core_methods,
};

PyMODINIT_FUNC
PyInit__core(void)
{
    return PyModule_Create(&core_module);
}
