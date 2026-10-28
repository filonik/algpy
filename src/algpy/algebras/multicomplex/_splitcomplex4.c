/*
 * splitcomplex4 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in splitcomplex4.h.
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

#include "splitcomplex4.h"

typedef struct {
    PyObject_HEAD
    splitcomplex4 v;
} PySplitcomplex4;

static PyObject *PySplitcomplex4_Type = NULL;
static int npy_splitcomplex4 = -1;

static PyObject *
PySplitcomplex4_FromValue(PyTypeObject *type, splitcomplex4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySplitcomplex4 *p = (PySplitcomplex4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SPLITCOMPLEX4(v) \
    PySplitcomplex4_FromValue((PyTypeObject *)PySplitcomplex4_Type, (v))

static int
splitcomplex4_from_object(PyObject *o, splitcomplex4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySplitcomplex4_Type)) {
        *out = ((PySplitcomplex4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = splitcomplex4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pysplitcomplex4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PySplitcomplex4_FromValue(type, splitcomplex4_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]));
}

static PyObject *
pysplitcomplex4_repr(PyObject *self)
{
    splitcomplex4 v = ((PySplitcomplex4 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("splitcomplex4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pysplitcomplex4_hash(PyObject *self)
{
    splitcomplex4 v = ((PySplitcomplex4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (splitcomplex4_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pysplitcomplex4_richcompare(PyObject *a, PyObject *b, int op)
{
    splitcomplex4 x, y;
    int r;
    if ((r = splitcomplex4_from_object(a, &x)) <= 0 ||
            (r = splitcomplex4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(splitcomplex4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!splitcomplex4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SPLITCOMPLEX4_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysplitcomplex4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        splitcomplex4 x, y;                                                       \
        int r;                                                               \
        if ((r = splitcomplex4_from_object(a, &x)) <= 0 ||                        \
                (r = splitcomplex4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SPLITCOMPLEX4(expr(x, y));                                     \
    }

SPLITCOMPLEX4_BINARY(add, splitcomplex4_add)
SPLITCOMPLEX4_BINARY(sub, splitcomplex4_sub)
SPLITCOMPLEX4_BINARY(mul, splitcomplex4_mul)
SPLITCOMPLEX4_BINARY(div, splitcomplex4_div)

static PyObject *
pysplitcomplex4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    splitcomplex4 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "splitcomplex4 does not support modular power");
        return NULL;
    }
    if ((r = splitcomplex4_from_object(a, &x)) <= 0 ||
            (r = splitcomplex4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SPLITCOMPLEX4(splitcomplex4_pow(x, y));
}

#define SPLITCOMPLEX4_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysplitcomplex4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SPLITCOMPLEX4(expr(((PySplitcomplex4 *)self)->v));                  \
    }

SPLITCOMPLEX4_UNARY(pos, splitcomplex4_pos)
SPLITCOMPLEX4_UNARY(neg, splitcomplex4_neg)

static PyObject *
pysplitcomplex4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue((PyTypeObject *)cls, splitcomplex4_zero());
}

static PyObject *
pysplitcomplex4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue((PyTypeObject *)cls, splitcomplex4_one());
}

static int
pysplitcomplex4_bool(PyObject *self)
{
    splitcomplex4 v = ((PySplitcomplex4 *)self)->v;
    return vec16_any_nonzero(v);
}

#define SPLITCOMPLEX4_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysplitcomplex4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SPLITCOMPLEX4(op(((PySplitcomplex4 *)self)->v));                   \
    }

SPLITCOMPLEX4_METHOD(conjugate, splitcomplex4_conjugate)
SPLITCOMPLEX4_METHOD(exp, splitcomplex4_exp)
SPLITCOMPLEX4_METHOD(log, splitcomplex4_log)
SPLITCOMPLEX4_METHOD(sqrt, splitcomplex4_sqrt)

static PyObject *
pysplitcomplex4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = splitcomplex4_norm(((PySplitcomplex4 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pysplitcomplex4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_SPLITCOMPLEX4(splitcomplex4_inv(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex4_lmat(((PySplitcomplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex4_rmat(((PySplitcomplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {17};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex4_minpoly(((PySplitcomplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[16];
    npy_intp dims[1] = {16};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex4_eivals(((PySplitcomplex4 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 16; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pysplitcomplex4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex4);
    splitcomplex4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex4 *)PyArray_DATA(arr);
    splitcomplex4_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex4_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_splitcomplex4);
    splitcomplex4 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (splitcomplex4 *)PyArray_DATA(arr);
    splitcomplex4_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex4 v = ((PySplitcomplex4 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pysplitcomplex4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_square(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_expm1(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_log1p(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_sin(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_cos(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_tan(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_asin(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_acos(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_atan(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_sinh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_cosh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_tanh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_asinh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_acosh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_atanh(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_re(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex4_FromValue(Py_TYPE(self), splitcomplex4_im(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_real(PyObject *self, void *closure)
{
    return pysplitcomplex4_re(self, NULL);
}

static PyObject *
pysplitcomplex4_imag(PyObject *self, void *closure)
{
    return pysplitcomplex4_im(self, NULL);
}

static PyGetSetDef pysplitcomplex4_getset[] = {
    {"real", pysplitcomplex4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pysplitcomplex4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pysplitcomplex4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(splitcomplex4_is_real(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(splitcomplex4_to_re(((PySplitcomplex4 *)self)->v));
}

static PyObject *
pysplitcomplex4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PySplitcomplex4_FromValue((PyTypeObject *)cls, splitcomplex4_from_re(a));
}

static PyObject *
pysplitcomplex4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    splitcomplex4 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = splitcomplex4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySplitcomplex4_FromValue(Py_TYPE(self),
        splitcomplex4_mix(((PySplitcomplex4 *)self)->v, b, t));
}

static PyObject *
pysplitcomplex4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_splitcomplex4, 16);
}

static PyObject *
pysplitcomplex4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_splitcomplex4, 16);
}

static PyMethodDef pysplitcomplex4_methods[] = {
    {"to_float", (PyCFunction)pysplitcomplex4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pysplitcomplex4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pysplitcomplex4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pysplitcomplex4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pysplitcomplex4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pysplitcomplex4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pysplitcomplex4_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pysplitcomplex4_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pysplitcomplex4_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sin", pysplitcomplex4_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pysplitcomplex4_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pysplitcomplex4_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pysplitcomplex4_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pysplitcomplex4_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pysplitcomplex4_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pysplitcomplex4_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pysplitcomplex4_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pysplitcomplex4_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pysplitcomplex4_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pysplitcomplex4_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pysplitcomplex4_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pysplitcomplex4_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pysplitcomplex4_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pysplitcomplex4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pysplitcomplex4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pysplitcomplex4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pysplitcomplex4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pysplitcomplex4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pysplitcomplex4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pysplitcomplex4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysplitcomplex4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysplitcomplex4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pysplitcomplex4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pysplitcomplex4_norm, METH_NOARGS, "Norm as a Python complex."},
    {"inv", pysplitcomplex4_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pysplitcomplex4_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pysplitcomplex4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pysplitcomplex4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pysplitcomplex4_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pysplitcomplex4_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysplitcomplex4_exp, METH_NOARGS, "exp."},
    {"log", pysplitcomplex4_log, METH_NOARGS, "log."},
    {"sqrt", pysplitcomplex4_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pysplitcomplex4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysplitcomplex4_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pysplitcomplex4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySplitcomplex4 *)self)->v.data[i]);
}

static PyObject *
pysplitcomplex4_subscript(PyObject *self, PyObject *key)
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
        i += 16;
    }
    return pysplitcomplex4_item(self, i);
}

static int
pysplitcomplex4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysplitcomplex4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysplitcomplex4_slots[] = {
    {Py_tp_dealloc, pysplitcomplex4_dealloc},
    {Py_tp_traverse, pysplitcomplex4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysplitcomplex4_repr},
    {Py_tp_str, pysplitcomplex4_repr},
    {Py_tp_hash, pysplitcomplex4_hash},
    {Py_tp_doc, (void *)"splitcomplex4 number with float64 components"},
    {Py_tp_richcompare, pysplitcomplex4_richcompare},
    {Py_tp_methods, pysplitcomplex4_methods},
    {Py_tp_getset, pysplitcomplex4_getset},
    {Py_sq_length, pysplitcomplex4_length},
    {Py_sq_item, pysplitcomplex4_item},
    {Py_mp_subscript, pysplitcomplex4_subscript},
    {Py_tp_new, pysplitcomplex4_new},
    {Py_nb_add, pysplitcomplex4_add},
    {Py_nb_subtract, pysplitcomplex4_sub},
    {Py_nb_multiply, pysplitcomplex4_mul},
    {Py_nb_true_divide, pysplitcomplex4_div},
    {Py_nb_power, (void *)pysplitcomplex4_pow},
    {Py_nb_positive, pysplitcomplex4_pos},
    {Py_nb_negative, pysplitcomplex4_neg},
    {Py_nb_bool, pysplitcomplex4_bool},
    {0, NULL}
};

static PyType_Spec pysplitcomplex4_spec = {
    .name = "algpy.algebras.multicomplex.splitcomplex4",
    .basicsize = sizeof(PySplitcomplex4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysplitcomplex4_slots,
};

static PyObject *
npysplitcomplex4_getitem(void *data, void *arr)
{
    splitcomplex4 v;
    memcpy(&v, data, sizeof(v));
    return NEW_SPLITCOMPLEX4(v);
}

static int
npysplitcomplex4_setitem(PyObject *item, void *data, void *arr)
{
    splitcomplex4 v;
    int r = splitcomplex4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to splitcomplex4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "splitcomplex4 needs exactly 16 components");
            return -1;
        }
        for (int i = 0; i < 16; i++) {
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
pysplitcomplex4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(splitcomplex4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 16; j++) {
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
pysplitcomplex4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysplitcomplex4_copyswapn(dst, sizeof(splitcomplex4), src, sizeof(splitcomplex4), 1,
                         swap, arr);
}

static npy_bool
pysplitcomplex4_nonzero(void *data, void *arr)
{
    splitcomplex4 v;
    memcpy(&v, data, sizeof(v));
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysplitcomplex4_arrfuncs;
typedef struct { char c; splitcomplex4 v; } pysplitcomplex4_align_test;

static PyArray_DescrProto pysplitcomplex4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(splitcomplex4), offsetof(pysplitcomplex4_align_test, v), 0, 0, 0,
    &pysplitcomplex4_arrfuncs,
};

#define SPLITCOMPLEX4_CAST_FROM(name, type)                                       \
    static void                                                              \
    pysplitcomplex4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        splitcomplex4 *dst = (splitcomplex4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = splitcomplex4_from_re((double)src[i]);                      \
        }                                                                    \
    }

#define SPLITCOMPLEX4_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

SPLITCOMPLEX4_REAL_TYPES(SPLITCOMPLEX4_CAST_FROM)

static int
pysplitcomplex4_register_casts(void)
{
#define REGISTER_SPLITCOMPLEX4_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_splitcomplex4,                    \
                                     pysplitcomplex4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_splitcomplex4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SPLITCOMPLEX4_REAL_TYPES(REGISTER_SPLITCOMPLEX4_CAST)
#undef REGISTER_SPLITCOMPLEX4_CAST
    return 0;
}

typedef splitcomplex4 (*splitcomplex4_unary_op)(splitcomplex4);
typedef splitcomplex4 (*splitcomplex4_binary_op)(splitcomplex4, splitcomplex4);
typedef int (*splitcomplex4_predicate)(splitcomplex4, splitcomplex4);

static void
loop_splitcomplex4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    splitcomplex4_unary_op op = (splitcomplex4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_splitcomplex4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    splitcomplex4_binary_op op = (splitcomplex4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex4 a, b, r;
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
loop_splitcomplex4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    splitcomplex4_predicate op = (splitcomplex4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysplitcomplex4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_splitcomplex4;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_splitcomplex4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysplitcomplex4_register_ufuncs(PyObject *numpy)
{
    return
        pysplitcomplex4_register_loop(numpy, "add", loop_splitcomplex4_binary, 3, 0, splitcomplex4_add) ||
        pysplitcomplex4_register_loop(numpy, "subtract", loop_splitcomplex4_binary, 3, 0, splitcomplex4_sub) ||
        pysplitcomplex4_register_loop(numpy, "multiply", loop_splitcomplex4_binary, 3, 0, splitcomplex4_mul) ||
        pysplitcomplex4_register_loop(numpy, "divide", loop_splitcomplex4_binary, 3, 0, splitcomplex4_div) ||
        pysplitcomplex4_register_loop(numpy, "power", loop_splitcomplex4_binary, 3, 0, splitcomplex4_pow) ||
        pysplitcomplex4_register_loop(numpy, "positive", loop_splitcomplex4_unary, 2, 0, splitcomplex4_pos) ||
        pysplitcomplex4_register_loop(numpy, "negative", loop_splitcomplex4_unary, 2, 0, splitcomplex4_neg) ||
        pysplitcomplex4_register_loop(numpy, "conjugate", loop_splitcomplex4_unary, 2, 0, splitcomplex4_conjugate) ||
        pysplitcomplex4_register_loop(numpy, "reciprocal", loop_splitcomplex4_unary, 2, 0, splitcomplex4_inv) ||
        pysplitcomplex4_register_loop(numpy, "square", loop_splitcomplex4_unary, 2, 0, splitcomplex4_square) ||
        pysplitcomplex4_register_loop(numpy, "exp", loop_splitcomplex4_unary, 2, 0, splitcomplex4_exp) ||
        pysplitcomplex4_register_loop(numpy, "expm1", loop_splitcomplex4_unary, 2, 0, splitcomplex4_expm1) ||
        pysplitcomplex4_register_loop(numpy, "log", loop_splitcomplex4_unary, 2, 0, splitcomplex4_log) ||
        pysplitcomplex4_register_loop(numpy, "log1p", loop_splitcomplex4_unary, 2, 0, splitcomplex4_log1p) ||
        pysplitcomplex4_register_loop(numpy, "sqrt", loop_splitcomplex4_unary, 2, 0, splitcomplex4_sqrt) ||
        pysplitcomplex4_register_loop(numpy, "sin", loop_splitcomplex4_unary, 2, 0, splitcomplex4_sin) ||
        pysplitcomplex4_register_loop(numpy, "cos", loop_splitcomplex4_unary, 2, 0, splitcomplex4_cos) ||
        pysplitcomplex4_register_loop(numpy, "tan", loop_splitcomplex4_unary, 2, 0, splitcomplex4_tan) ||
        pysplitcomplex4_register_loop(numpy, "arcsin", loop_splitcomplex4_unary, 2, 0, splitcomplex4_asin) ||
        pysplitcomplex4_register_loop(numpy, "arccos", loop_splitcomplex4_unary, 2, 0, splitcomplex4_acos) ||
        pysplitcomplex4_register_loop(numpy, "arctan", loop_splitcomplex4_unary, 2, 0, splitcomplex4_atan) ||
        pysplitcomplex4_register_loop(numpy, "sinh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_sinh) ||
        pysplitcomplex4_register_loop(numpy, "cosh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_cosh) ||
        pysplitcomplex4_register_loop(numpy, "tanh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_tanh) ||
        pysplitcomplex4_register_loop(numpy, "arcsinh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_asinh) ||
        pysplitcomplex4_register_loop(numpy, "arccosh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_acosh) ||
        pysplitcomplex4_register_loop(numpy, "arctanh", loop_splitcomplex4_unary, 2, 0, splitcomplex4_atanh) ||
        pysplitcomplex4_register_loop(numpy, "equal", loop_splitcomplex4_predicate, 3, 1, splitcomplex4_equal) ||
        pysplitcomplex4_register_loop(numpy, "not_equal", loop_splitcomplex4_predicate, 3, 1, splitcomplex4_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_splitcomplex4", "splitcomplex4 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__splitcomplex4(void)
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
    PySplitcomplex4_Type = PyType_FromModuleAndSpec(
        m, &pysplitcomplex4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySplitcomplex4_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysplitcomplex4_arrfuncs);
    pysplitcomplex4_arrfuncs.getitem = npysplitcomplex4_getitem;
    pysplitcomplex4_arrfuncs.setitem = npysplitcomplex4_setitem;
    pysplitcomplex4_arrfuncs.copyswapn = pysplitcomplex4_copyswapn;
    pysplitcomplex4_arrfuncs.copyswap = pysplitcomplex4_copyswap;
    pysplitcomplex4_arrfuncs.nonzero = pysplitcomplex4_nonzero;
    pysplitcomplex4_descr_proto.typeobj = (PyTypeObject *)PySplitcomplex4_Type;
    Py_SET_TYPE((PyObject *)&pysplitcomplex4_descr_proto, &PyArrayDescr_Type);
    npy_splitcomplex4 = PyArray_RegisterDataType(&pysplitcomplex4_descr_proto);
    if (npy_splitcomplex4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_splitcomplex4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySplitcomplex4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PySplitcomplex4_Type, 16) < 0) {
        goto fail;
    }
    if (pysplitcomplex4_register_casts() < 0 ||
            pysplitcomplex4_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "splitcomplex4", PySplitcomplex4_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
