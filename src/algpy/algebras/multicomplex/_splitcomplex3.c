/*
 * splitcomplex3 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in splitcomplex3.h.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "../../_conversions.h"

#include "splitcomplex3.h"

typedef struct {
    PyObject_HEAD
    splitcomplex3 v;
} PySplitcomplex3;

static PyObject *PySplitcomplex3_Type = NULL;
static int npy_splitcomplex3 = -1;

static PyObject *
PySplitcomplex3_FromValue(PyTypeObject *type, splitcomplex3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySplitcomplex3 *p = (PySplitcomplex3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SPLITCOMPLEX3(v) \
    PySplitcomplex3_FromValue((PyTypeObject *)PySplitcomplex3_Type, (v))

static int
splitcomplex3_from_object(PyObject *o, splitcomplex3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySplitcomplex3_Type)) {
        *out = ((PySplitcomplex3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = splitcomplex3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pysplitcomplex3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PySplitcomplex3_FromValue(type, splitcomplex3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]));
}

static PyObject *
pysplitcomplex3_repr(PyObject *self)
{
    splitcomplex3 v = ((PySplitcomplex3 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("splitcomplex3(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pysplitcomplex3_hash(PyObject *self)
{
    splitcomplex3 v = ((PySplitcomplex3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (splitcomplex3_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pysplitcomplex3_richcompare(PyObject *a, PyObject *b, int op)
{
    splitcomplex3 x, y;
    int r;
    if ((r = splitcomplex3_from_object(a, &x)) <= 0 ||
            (r = splitcomplex3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(splitcomplex3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!splitcomplex3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SPLITCOMPLEX3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysplitcomplex3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        splitcomplex3 x, y;                                                       \
        int r;                                                               \
        if ((r = splitcomplex3_from_object(a, &x)) <= 0 ||                        \
                (r = splitcomplex3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SPLITCOMPLEX3(expr(x, y));                                     \
    }

SPLITCOMPLEX3_BINARY(add, splitcomplex3_add)
SPLITCOMPLEX3_BINARY(sub, splitcomplex3_sub)
SPLITCOMPLEX3_BINARY(mul, splitcomplex3_mul)
SPLITCOMPLEX3_BINARY(div, splitcomplex3_div)

static PyObject *
pysplitcomplex3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    splitcomplex3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "splitcomplex3 does not support modular power");
        return NULL;
    }
    if ((r = splitcomplex3_from_object(a, &x)) <= 0 ||
            (r = splitcomplex3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SPLITCOMPLEX3(splitcomplex3_pow(x, y));
}

#define SPLITCOMPLEX3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysplitcomplex3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SPLITCOMPLEX3(expr(((PySplitcomplex3 *)self)->v));                  \
    }

SPLITCOMPLEX3_UNARY(pos, splitcomplex3_pos)
SPLITCOMPLEX3_UNARY(neg, splitcomplex3_neg)

static PyObject *
pysplitcomplex3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue((PyTypeObject *)cls, splitcomplex3_zero());
}

static PyObject *
pysplitcomplex3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue((PyTypeObject *)cls, splitcomplex3_one());
}

static int
pysplitcomplex3_bool(PyObject *self)
{
    splitcomplex3 v = ((PySplitcomplex3 *)self)->v;
    return vec8_any_nonzero(v);
}

#define SPLITCOMPLEX3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysplitcomplex3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SPLITCOMPLEX3(op(((PySplitcomplex3 *)self)->v));                   \
    }

SPLITCOMPLEX3_METHOD(conjugate, splitcomplex3_conjugate)
SPLITCOMPLEX3_METHOD(exp, splitcomplex3_exp)
SPLITCOMPLEX3_METHOD(log, splitcomplex3_log)
SPLITCOMPLEX3_METHOD(sqrt, splitcomplex3_sqrt)

static PyObject *
pysplitcomplex3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = splitcomplex3_norm(((PySplitcomplex3 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pysplitcomplex3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_SPLITCOMPLEX3(splitcomplex3_inv(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex3_lmat(((PySplitcomplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex3_rmat(((PySplitcomplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {9};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex3_minpoly(((PySplitcomplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[8];
    npy_intp dims[1] = {8};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex3_eivals(((PySplitcomplex3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pysplitcomplex3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex3);
    splitcomplex3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex3 *)PyArray_DATA(arr);
    splitcomplex3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_splitcomplex3);
    splitcomplex3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (splitcomplex3 *)PyArray_DATA(arr);
    splitcomplex3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex3 v = ((PySplitcomplex3 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pysplitcomplex3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_square(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_expm1(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_log1p(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_sin(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_cos(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_tan(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_asin(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_acos(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_atan(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_sinh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_cosh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_tanh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_asinh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_acosh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_atanh(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_re(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex3_FromValue(Py_TYPE(self), splitcomplex3_im(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_real(PyObject *self, void *closure)
{
    return pysplitcomplex3_re(self, NULL);
}

static PyObject *
pysplitcomplex3_imag(PyObject *self, void *closure)
{
    return pysplitcomplex3_im(self, NULL);
}

static PyGetSetDef pysplitcomplex3_getset[] = {
    {"real", pysplitcomplex3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pysplitcomplex3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pysplitcomplex3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(splitcomplex3_is_real(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(splitcomplex3_to_re(((PySplitcomplex3 *)self)->v));
}

static PyObject *
pysplitcomplex3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PySplitcomplex3_FromValue((PyTypeObject *)cls, splitcomplex3_from_re(a));
}

static PyObject *
pysplitcomplex3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    splitcomplex3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = splitcomplex3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySplitcomplex3_FromValue(Py_TYPE(self),
        splitcomplex3_mix(((PySplitcomplex3 *)self)->v, b, t));
}

static PyObject *
pysplitcomplex3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_splitcomplex3, 8);
}

static PyObject *
pysplitcomplex3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_splitcomplex3, 8);
}

static PyMethodDef pysplitcomplex3_methods[] = {
    {"to_float", (PyCFunction)pysplitcomplex3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pysplitcomplex3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pysplitcomplex3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pysplitcomplex3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pysplitcomplex3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pysplitcomplex3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pysplitcomplex3_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pysplitcomplex3_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pysplitcomplex3_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sin", pysplitcomplex3_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pysplitcomplex3_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pysplitcomplex3_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pysplitcomplex3_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pysplitcomplex3_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pysplitcomplex3_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pysplitcomplex3_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pysplitcomplex3_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pysplitcomplex3_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pysplitcomplex3_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pysplitcomplex3_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pysplitcomplex3_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pysplitcomplex3_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pysplitcomplex3_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pysplitcomplex3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pysplitcomplex3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pysplitcomplex3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pysplitcomplex3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pysplitcomplex3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pysplitcomplex3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pysplitcomplex3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysplitcomplex3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysplitcomplex3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pysplitcomplex3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pysplitcomplex3_norm, METH_NOARGS, "Norm as a Python complex."},
    {"inv", pysplitcomplex3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pysplitcomplex3_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pysplitcomplex3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pysplitcomplex3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pysplitcomplex3_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pysplitcomplex3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysplitcomplex3_exp, METH_NOARGS, "exp."},
    {"log", pysplitcomplex3_log, METH_NOARGS, "log."},
    {"sqrt", pysplitcomplex3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pysplitcomplex3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysplitcomplex3_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pysplitcomplex3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySplitcomplex3 *)self)->v.data[i]);
}

static PyObject *
pysplitcomplex3_subscript(PyObject *self, PyObject *key)
{
    Py_ssize_t i;
    if (!PyIndex_Check(key)) {
        PyErr_SetString(PyExc_TypeError, "indices must be integers");
        return NULL;
    }
    i = PyNumber_AsSsize_t(key, PyExc_IndexError);
    if (i == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (i < 0) {
        i += 8;
    }
    return pysplitcomplex3_item(self, i);
}

static int
pysplitcomplex3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysplitcomplex3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysplitcomplex3_slots[] = {
    {Py_tp_dealloc, pysplitcomplex3_dealloc},
    {Py_tp_traverse, pysplitcomplex3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysplitcomplex3_repr},
    {Py_tp_str, pysplitcomplex3_repr},
    {Py_tp_hash, pysplitcomplex3_hash},
    {Py_tp_doc, (void *)"splitcomplex3 number with float64 components"},
    {Py_tp_richcompare, pysplitcomplex3_richcompare},
    {Py_tp_methods, pysplitcomplex3_methods},
    {Py_tp_getset, pysplitcomplex3_getset},
    {Py_sq_length, pysplitcomplex3_length},
    {Py_sq_item, pysplitcomplex3_item},
    {Py_mp_subscript, pysplitcomplex3_subscript},
    {Py_tp_new, pysplitcomplex3_new},
    {Py_nb_add, pysplitcomplex3_add},
    {Py_nb_subtract, pysplitcomplex3_sub},
    {Py_nb_multiply, pysplitcomplex3_mul},
    {Py_nb_true_divide, pysplitcomplex3_div},
    {Py_nb_power, (void *)pysplitcomplex3_pow},
    {Py_nb_positive, pysplitcomplex3_pos},
    {Py_nb_negative, pysplitcomplex3_neg},
    {Py_nb_bool, pysplitcomplex3_bool},
    {0, NULL}
};

static PyType_Spec pysplitcomplex3_spec = {
    .name = "algpy.algebras.multicomplex.splitcomplex3",
    .basicsize = sizeof(PySplitcomplex3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysplitcomplex3_slots,
};

static PyObject *
npysplitcomplex3_getitem(void *data, void *arr)
{
    splitcomplex3 v;
    memcpy(&v, data, sizeof(v));
    return NEW_SPLITCOMPLEX3(v);
}

static int
npysplitcomplex3_setitem(PyObject *item, void *data, void *arr)
{
    splitcomplex3 v;
    int r = splitcomplex3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to splitcomplex3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "splitcomplex3 needs exactly 8 components");
            return -1;
        }
        for (int i = 0; i < 8; i++) {
            v.data[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        }
        Py_DECREF(seq);
        if (PyErr_Occurred()) {
            return -1;
        }
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pysplitcomplex3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(splitcomplex3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 8; j++) {
                char *d = p + 8 * j;
                int b;
                for (b = 0; b < 4; b++) {
                    char t = d[b];
                    d[b] = d[7 - b];
                    d[7 - b] = t;
                }
            }
        }
    }
}

static void
pysplitcomplex3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysplitcomplex3_copyswapn(dst, sizeof(splitcomplex3), src, sizeof(splitcomplex3), 1,
                         swap, arr);
}

static npy_bool
pysplitcomplex3_nonzero(void *data, void *arr)
{
    splitcomplex3 v;
    memcpy(&v, data, sizeof(v));
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysplitcomplex3_arrfuncs;
typedef struct { char c; splitcomplex3 v; } pysplitcomplex3_align_test;

static PyArray_DescrProto pysplitcomplex3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(splitcomplex3), offsetof(pysplitcomplex3_align_test, v), 0, 0, 0,
    &pysplitcomplex3_arrfuncs,
};

#define SPLITCOMPLEX3_CAST_FROM(name, type)                                       \
    static void                                                              \
    pysplitcomplex3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        splitcomplex3 *dst = (splitcomplex3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = splitcomplex3_from_re((double)src[i]);                      \
        }                                                                    \
    }

#define SPLITCOMPLEX3_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

SPLITCOMPLEX3_REAL_TYPES(SPLITCOMPLEX3_CAST_FROM)

static int
pysplitcomplex3_register_casts(void)
{
#define REGISTER_SPLITCOMPLEX3_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_splitcomplex3,                    \
                                     pysplitcomplex3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_splitcomplex3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SPLITCOMPLEX3_REAL_TYPES(REGISTER_SPLITCOMPLEX3_CAST)
#undef REGISTER_SPLITCOMPLEX3_CAST
    return 0;
}

typedef splitcomplex3 (*splitcomplex3_unary_op)(splitcomplex3);
typedef splitcomplex3 (*splitcomplex3_binary_op)(splitcomplex3, splitcomplex3);
typedef int (*splitcomplex3_predicate)(splitcomplex3, splitcomplex3);

static void
loop_splitcomplex3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    splitcomplex3_unary_op op = (splitcomplex3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_splitcomplex3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    splitcomplex3_binary_op op = (splitcomplex3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex3 a, b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = op(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_splitcomplex3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    splitcomplex3_predicate op = (splitcomplex3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysplitcomplex3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_splitcomplex3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_splitcomplex3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysplitcomplex3_register_ufuncs(PyObject *numpy)
{
    return
        pysplitcomplex3_register_loop(numpy, "add", loop_splitcomplex3_binary, 3, 0, splitcomplex3_add) ||
        pysplitcomplex3_register_loop(numpy, "subtract", loop_splitcomplex3_binary, 3, 0, splitcomplex3_sub) ||
        pysplitcomplex3_register_loop(numpy, "multiply", loop_splitcomplex3_binary, 3, 0, splitcomplex3_mul) ||
        pysplitcomplex3_register_loop(numpy, "divide", loop_splitcomplex3_binary, 3, 0, splitcomplex3_div) ||
        pysplitcomplex3_register_loop(numpy, "power", loop_splitcomplex3_binary, 3, 0, splitcomplex3_pow) ||
        pysplitcomplex3_register_loop(numpy, "positive", loop_splitcomplex3_unary, 2, 0, splitcomplex3_pos) ||
        pysplitcomplex3_register_loop(numpy, "negative", loop_splitcomplex3_unary, 2, 0, splitcomplex3_neg) ||
        pysplitcomplex3_register_loop(numpy, "conjugate", loop_splitcomplex3_unary, 2, 0, splitcomplex3_conjugate) ||
        pysplitcomplex3_register_loop(numpy, "reciprocal", loop_splitcomplex3_unary, 2, 0, splitcomplex3_inv) ||
        pysplitcomplex3_register_loop(numpy, "square", loop_splitcomplex3_unary, 2, 0, splitcomplex3_square) ||
        pysplitcomplex3_register_loop(numpy, "exp", loop_splitcomplex3_unary, 2, 0, splitcomplex3_exp) ||
        pysplitcomplex3_register_loop(numpy, "expm1", loop_splitcomplex3_unary, 2, 0, splitcomplex3_expm1) ||
        pysplitcomplex3_register_loop(numpy, "log", loop_splitcomplex3_unary, 2, 0, splitcomplex3_log) ||
        pysplitcomplex3_register_loop(numpy, "log1p", loop_splitcomplex3_unary, 2, 0, splitcomplex3_log1p) ||
        pysplitcomplex3_register_loop(numpy, "sqrt", loop_splitcomplex3_unary, 2, 0, splitcomplex3_sqrt) ||
        pysplitcomplex3_register_loop(numpy, "sin", loop_splitcomplex3_unary, 2, 0, splitcomplex3_sin) ||
        pysplitcomplex3_register_loop(numpy, "cos", loop_splitcomplex3_unary, 2, 0, splitcomplex3_cos) ||
        pysplitcomplex3_register_loop(numpy, "tan", loop_splitcomplex3_unary, 2, 0, splitcomplex3_tan) ||
        pysplitcomplex3_register_loop(numpy, "arcsin", loop_splitcomplex3_unary, 2, 0, splitcomplex3_asin) ||
        pysplitcomplex3_register_loop(numpy, "arccos", loop_splitcomplex3_unary, 2, 0, splitcomplex3_acos) ||
        pysplitcomplex3_register_loop(numpy, "arctan", loop_splitcomplex3_unary, 2, 0, splitcomplex3_atan) ||
        pysplitcomplex3_register_loop(numpy, "sinh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_sinh) ||
        pysplitcomplex3_register_loop(numpy, "cosh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_cosh) ||
        pysplitcomplex3_register_loop(numpy, "tanh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_tanh) ||
        pysplitcomplex3_register_loop(numpy, "arcsinh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_asinh) ||
        pysplitcomplex3_register_loop(numpy, "arccosh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_acosh) ||
        pysplitcomplex3_register_loop(numpy, "arctanh", loop_splitcomplex3_unary, 2, 0, splitcomplex3_atanh) ||
        pysplitcomplex3_register_loop(numpy, "equal", loop_splitcomplex3_predicate, 3, 1, splitcomplex3_equal) ||
        pysplitcomplex3_register_loop(numpy, "not_equal", loop_splitcomplex3_predicate, 3, 1, splitcomplex3_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_splitcomplex3", "splitcomplex3 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__splitcomplex3(void)
{
    PyObject *m, *numpy;
    PyArray_Descr *descr;
    import_array();
    import_umath();
    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }
    numpy = PyImport_ImportModule("numpy");
    if (numpy == NULL) {
        goto fail;
    }
    PySplitcomplex3_Type = PyType_FromModuleAndSpec(
        m, &pysplitcomplex3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySplitcomplex3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysplitcomplex3_arrfuncs);
    pysplitcomplex3_arrfuncs.getitem = npysplitcomplex3_getitem;
    pysplitcomplex3_arrfuncs.setitem = npysplitcomplex3_setitem;
    pysplitcomplex3_arrfuncs.copyswapn = pysplitcomplex3_copyswapn;
    pysplitcomplex3_arrfuncs.copyswap = pysplitcomplex3_copyswap;
    pysplitcomplex3_arrfuncs.nonzero = pysplitcomplex3_nonzero;
    pysplitcomplex3_descr_proto.typeobj = (PyTypeObject *)PySplitcomplex3_Type;
    Py_SET_TYPE((PyObject *)&pysplitcomplex3_descr_proto, &PyArrayDescr_Type);
    npy_splitcomplex3 = PyArray_RegisterDataType(&pysplitcomplex3_descr_proto);
    if (npy_splitcomplex3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_splitcomplex3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySplitcomplex3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PySplitcomplex3_Type, 8) < 0) {
        goto fail;
    }
    if (pysplitcomplex3_register_casts() < 0 ||
            pysplitcomplex3_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "splitcomplex3", PySplitcomplex3_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
