/*
 * splitcomplex5 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in splitcomplex5.h.
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
#include "../../_grades.h"

#include "splitcomplex5.h"

typedef struct {
    PyObject_HEAD
    splitcomplex5 v;
} PySplitcomplex5;

static PyObject *PySplitcomplex5_Type = NULL;
static int npy_splitcomplex5 = -1;

static PyObject *
PySplitcomplex5_FromValue(PyTypeObject *type, splitcomplex5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySplitcomplex5 *p = (PySplitcomplex5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SPLITCOMPLEX5(v) \
    PySplitcomplex5_FromValue((PyTypeObject *)PySplitcomplex5_Type, (v))

static int
splitcomplex5_from_object(PyObject *o, splitcomplex5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySplitcomplex5_Type)) {
        *out = ((PySplitcomplex5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = splitcomplex5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pysplitcomplex5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PySplitcomplex5_FromValue(type, splitcomplex5_make(d));
}

static PyObject *
pysplitcomplex5_repr(PyObject *self)
{
    splitcomplex5 v = ((PySplitcomplex5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("splitcomplex5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pysplitcomplex5_hash(PyObject *self)
{
    splitcomplex5 v = ((PySplitcomplex5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_SPLITCOMPLEX5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (splitcomplex5_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddddddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pysplitcomplex5_richcompare(PyObject *a, PyObject *b, int op)
{
    splitcomplex5 x, y;
    int r;
    if ((r = splitcomplex5_from_object(a, &x)) <= 0 ||
            (r = splitcomplex5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(splitcomplex5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!splitcomplex5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SPLITCOMPLEX5_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysplitcomplex5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        splitcomplex5 x, y;                                                       \
        int r;                                                               \
        if ((r = splitcomplex5_from_object(a, &x)) <= 0 ||                        \
                (r = splitcomplex5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SPLITCOMPLEX5(expr(x, y));                                     \
    }

SPLITCOMPLEX5_BINARY(add, splitcomplex5_add)
SPLITCOMPLEX5_BINARY(sub, splitcomplex5_sub)
SPLITCOMPLEX5_BINARY(mul, splitcomplex5_mul)
SPLITCOMPLEX5_BINARY(div, splitcomplex5_div)

static PyObject *
pysplitcomplex5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    splitcomplex5 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "splitcomplex5 does not support modular power");
        return NULL;
    }
    if ((r = splitcomplex5_from_object(a, &x)) <= 0 ||
            (r = splitcomplex5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SPLITCOMPLEX5(splitcomplex5_rpow(x, y));
}

#define SPLITCOMPLEX5_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysplitcomplex5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SPLITCOMPLEX5(expr(((PySplitcomplex5 *)self)->v));                  \
    }

SPLITCOMPLEX5_UNARY(pos, splitcomplex5_pos)
SPLITCOMPLEX5_UNARY(neg, splitcomplex5_neg)

static PyObject *
pysplitcomplex5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue((PyTypeObject *)cls, splitcomplex5_zero());
}

static PyObject *
pysplitcomplex5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue((PyTypeObject *)cls, splitcomplex5_one());
}

static int
pysplitcomplex5_bool(PyObject *self)
{
    splitcomplex5 v = ((PySplitcomplex5 *)self)->v;
    return vec32_any_nonzero(v);
}

#define SPLITCOMPLEX5_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysplitcomplex5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SPLITCOMPLEX5(op(((PySplitcomplex5 *)self)->v));                   \
    }

SPLITCOMPLEX5_METHOD(conjugate, splitcomplex5_conjugate)
SPLITCOMPLEX5_METHOD(exp, splitcomplex5_exp)
SPLITCOMPLEX5_METHOD(log, splitcomplex5_log)
SPLITCOMPLEX5_METHOD(sqrt, splitcomplex5_sqrt)

static PyObject *
pysplitcomplex5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex value = splitcomplex5_norm(((PySplitcomplex5 *)self)->v);
    return PyComplex_FromDoubles(value.data[0], value.data[1]);
}

static PyObject *
pysplitcomplex5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_SPLITCOMPLEX5(splitcomplex5_inv(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex5_lmat(((PySplitcomplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex5_rmat(((PySplitcomplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {33};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex5_minpoly(((PySplitcomplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[32];
    npy_intp dims[1] = {32};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex5_eivals(((PySplitcomplex5 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 32; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pysplitcomplex5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex5);
    splitcomplex5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex5 *)PyArray_DATA(arr);
    splitcomplex5_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex5);
    splitcomplex5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex5 *)PyArray_DATA(arr);
    splitcomplex5_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex5_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_splitcomplex5);
    splitcomplex5 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (splitcomplex5 *)PyArray_DATA(arr);
    splitcomplex5_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex5 v = ((PySplitcomplex5 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pysplitcomplex5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_square(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_expm1(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_log1p(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_sin(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_cos(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_tan(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_asin(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_acos(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_atan(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_sinh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_cosh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_tanh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_asinh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_acosh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_atanh(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_re(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_im(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_real(PyObject *self, void *closure)
{
    return pysplitcomplex5_re(self, NULL);
}

static PyObject *
pysplitcomplex5_imag(PyObject *self, void *closure)
{
    return pysplitcomplex5_im(self, NULL);
}

static PyObject *pysplitcomplex5_grades(PyObject *self, void *closure);

static PyGetSetDef pysplitcomplex5_getset[] = {
    {"grade", pysplitcomplex5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pysplitcomplex5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pysplitcomplex5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pysplitcomplex5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(splitcomplex5_is_real(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(splitcomplex5_to_re(((PySplitcomplex5 *)self)->v));
}

static PyObject *
pysplitcomplex5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PySplitcomplex5_FromValue((PyTypeObject *)cls, splitcomplex5_from_re(a));
}

static PyObject *
pysplitcomplex5_lpow(PyObject *self, PyObject *other)
{
    splitcomplex5 b;
    int r;
    r = splitcomplex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySplitcomplex5_FromValue(Py_TYPE(self),
        splitcomplex5_lpow(((PySplitcomplex5 *)self)->v, b));
}

static PyObject *
pysplitcomplex5_rpow(PyObject *self, PyObject *other)
{
    splitcomplex5 b;
    int r;
    r = splitcomplex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySplitcomplex5_FromValue(Py_TYPE(self),
        splitcomplex5_rpow(((PySplitcomplex5 *)self)->v, b));
}

static PyObject *
pysplitcomplex5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    splitcomplex5 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = splitcomplex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySplitcomplex5_FromValue(Py_TYPE(self),
        splitcomplex5_mix(((PySplitcomplex5 *)self)->v, b, t));
}

static PyObject *
pysplitcomplex5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_splitcomplex5, 32);
}

static PyObject *
pysplitcomplex5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_splitcomplex5, 32);
}

static PyObject *
pysplitcomplex5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_SPLITCOMPLEX5(splitcomplex5_ipow(((PySplitcomplex5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pysplitcomplex5_grade(PyObject *self, Py_ssize_t k)
{
    splitcomplex5 v = ((PySplitcomplex5 *)self)->v;
    switch (k) {
        case 0:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade0(v));
        case 1:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade1(v));
        case 2:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade2(v));
        case 3:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade3(v));
        case 4:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade4(v));
        case 5:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_grade5(v));
        default:
            return PySplitcomplex5_FromValue(Py_TYPE(self), splitcomplex5_zero());
    }
}

static PyObject *
pysplitcomplex5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pysplitcomplex5_grade, 6);
}

static PyMethodDef pysplitcomplex5_methods[] = {
    {"to_float", (PyCFunction)pysplitcomplex5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pysplitcomplex5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pysplitcomplex5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pysplitcomplex5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pysplitcomplex5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pysplitcomplex5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pysplitcomplex5_square, METH_NOARGS, "square."},
    {"expm1", pysplitcomplex5_expm1, METH_NOARGS, "expm1."},
    {"log1p", pysplitcomplex5_log1p, METH_NOARGS, "log1p."},
    {"sin", pysplitcomplex5_sin, METH_NOARGS, "sin."},
    {"cos", pysplitcomplex5_cos, METH_NOARGS, "cos."},
    {"tan", pysplitcomplex5_tan, METH_NOARGS, "tan."},
    {"asin", pysplitcomplex5_asin, METH_NOARGS, "asin."},
    {"acos", pysplitcomplex5_acos, METH_NOARGS, "acos."},
    {"atan", pysplitcomplex5_atan, METH_NOARGS, "atan."},
    {"sinh", pysplitcomplex5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pysplitcomplex5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pysplitcomplex5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pysplitcomplex5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pysplitcomplex5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pysplitcomplex5_atanh, METH_NOARGS, "atanh."},
    {"re", pysplitcomplex5_re, METH_NOARGS, "re."},
    {"im", pysplitcomplex5_im, METH_NOARGS, "im."},
    {"arcsin", pysplitcomplex5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pysplitcomplex5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pysplitcomplex5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pysplitcomplex5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pysplitcomplex5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pysplitcomplex5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pysplitcomplex5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysplitcomplex5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysplitcomplex5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pysplitcomplex5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pysplitcomplex5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pysplitcomplex5_norm, METH_NOARGS, "Norm as a Python complex."},
    {"lpow", pysplitcomplex5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pysplitcomplex5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pysplitcomplex5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pysplitcomplex5_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pysplitcomplex5_conjugate, METH_NOARGS, "Conjugation of the outer split-complex generator."},
    {"lmat", pysplitcomplex5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pysplitcomplex5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pysplitcomplex5_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pysplitcomplex5_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysplitcomplex5_exp, METH_NOARGS, "exp."},
    {"log", pysplitcomplex5_log, METH_NOARGS, "log."},
    {"sqrt", pysplitcomplex5_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pysplitcomplex5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysplitcomplex5_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pysplitcomplex5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySplitcomplex5 *)self)->v.data[i]);
}

static PyObject *
pysplitcomplex5_subscript(PyObject *self, PyObject *key)
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
        i += 32;
    }
    return pysplitcomplex5_item(self, i);
}

static int
pysplitcomplex5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysplitcomplex5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysplitcomplex5_slots[] = {
    {Py_tp_dealloc, pysplitcomplex5_dealloc},
    {Py_tp_traverse, pysplitcomplex5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysplitcomplex5_repr},
    {Py_tp_str, pysplitcomplex5_repr},
    {Py_tp_hash, pysplitcomplex5_hash},
    {Py_tp_doc, (void *)"Fivefold split-complex number with 32 float64 coefficients"},
    {Py_tp_richcompare, pysplitcomplex5_richcompare},
    {Py_tp_methods, pysplitcomplex5_methods},
    {Py_tp_getset, pysplitcomplex5_getset},
    {Py_sq_length, pysplitcomplex5_length},
    {Py_sq_item, pysplitcomplex5_item},
    {Py_mp_subscript, pysplitcomplex5_subscript},
    {Py_tp_new, pysplitcomplex5_new},
    {Py_nb_add, pysplitcomplex5_add},
    {Py_nb_subtract, pysplitcomplex5_sub},
    {Py_nb_multiply, pysplitcomplex5_mul},
    {Py_nb_true_divide, pysplitcomplex5_div},
    {Py_nb_power, (void *)pysplitcomplex5_pow},
    {Py_nb_positive, pysplitcomplex5_pos},
    {Py_nb_negative, pysplitcomplex5_neg},
    {Py_nb_bool, pysplitcomplex5_bool},
    {0, NULL}
};

static PyType_Spec pysplitcomplex5_spec = {
    .name = "algpy.algebras.multicomplex.splitcomplex5",
    .basicsize = sizeof(PySplitcomplex5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysplitcomplex5_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npysplitcomplex5_getitem(void *data, void *arr)
{
    splitcomplex5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_SPLITCOMPLEX5(v);
}

static int
npysplitcomplex5_setitem(PyObject *item, void *data, void *arr)
{
    splitcomplex5 v;
    int r = splitcomplex5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to splitcomplex5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "splitcomplex5 needs exactly 32 components");
            return -1;
        }
        for (int i = 0; i < 32; i++) {
            v.data[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        }
        Py_DECREF(seq);
        if (PyErr_Occurred()) {
            return -1;
        }
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pysplitcomplex5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(splitcomplex5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pysplitcomplex5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysplitcomplex5_copyswapn(dst, sizeof(splitcomplex5), src, sizeof(splitcomplex5), 1,
                         swap, arr);
}

static npy_bool
pysplitcomplex5_nonzero(void *data, void *arr)
{
    splitcomplex5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysplitcomplex5_arrfuncs;
typedef struct { char c; splitcomplex5 v; } pysplitcomplex5_align_test;

static PyArray_DescrProto pysplitcomplex5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(splitcomplex5), offsetof(pysplitcomplex5_align_test, v), 0, 0, 0,
    &pysplitcomplex5_arrfuncs,
};

#define SPLITCOMPLEX5_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pysplitcomplex5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        splitcomplex5 *dst = (splitcomplex5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = splitcomplex5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define SPLITCOMPLEX5_REAL_TYPES(X) \
    X(BOOL, npy_bool, (double)) \
    X(BYTE, npy_byte, (double)) \
    X(UBYTE, npy_ubyte, (double)) \
    X(SHORT, npy_short, (double)) \
    X(USHORT, npy_ushort, (double)) \
    X(INT, npy_int, (double)) \
    X(UINT, npy_uint, (double)) \
    X(LONG, npy_long, (double)) \
    X(ULONG, npy_ulong, (double)) \
    X(LONGLONG, npy_longlong, (double)) \
    X(ULONGLONG, npy_ulonglong, (double)) \
    X(HALF, npy_half, algpy_half_to_double) \
    X(FLOAT, npy_float, (double)) \
    X(DOUBLE, npy_double, (double)) \
    X(LONGDOUBLE, npy_longdouble, (double))

SPLITCOMPLEX5_REAL_TYPES(SPLITCOMPLEX5_CAST_FROM)

static int
pysplitcomplex5_register_casts(void)
{
#define REGISTER_SPLITCOMPLEX5_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_splitcomplex5,                    \
                                     pysplitcomplex5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_splitcomplex5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SPLITCOMPLEX5_REAL_TYPES(REGISTER_SPLITCOMPLEX5_CAST)
#undef REGISTER_SPLITCOMPLEX5_CAST
    return 0;
}

typedef splitcomplex5 (*splitcomplex5_unary_op)(splitcomplex5);
typedef splitcomplex5 (*splitcomplex5_binary_op)(splitcomplex5, splitcomplex5);
typedef int (*splitcomplex5_predicate)(splitcomplex5, splitcomplex5);

static void
loop_splitcomplex5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    splitcomplex5_unary_op op = (splitcomplex5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_splitcomplex5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    splitcomplex5_binary_op op = (splitcomplex5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex5 a, b, r;
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
loop_splitcomplex5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    splitcomplex5_predicate op = (splitcomplex5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysplitcomplex5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_splitcomplex5;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_splitcomplex5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysplitcomplex5_register_ufuncs(PyObject *numpy)
{
    return
        pysplitcomplex5_register_loop(numpy, "add", loop_splitcomplex5_binary, 3, 0, splitcomplex5_add) ||
        pysplitcomplex5_register_loop(numpy, "subtract", loop_splitcomplex5_binary, 3, 0, splitcomplex5_sub) ||
        pysplitcomplex5_register_loop(numpy, "multiply", loop_splitcomplex5_binary, 3, 0, splitcomplex5_mul) ||
        pysplitcomplex5_register_loop(numpy, "divide", loop_splitcomplex5_binary, 3, 0, splitcomplex5_div) ||
        pysplitcomplex5_register_loop(numpy, "power", loop_splitcomplex5_binary, 3, 0, splitcomplex5_rpow) ||
        pysplitcomplex5_register_loop(numpy, "positive", loop_splitcomplex5_unary, 2, 0, splitcomplex5_pos) ||
        pysplitcomplex5_register_loop(numpy, "negative", loop_splitcomplex5_unary, 2, 0, splitcomplex5_neg) ||
        pysplitcomplex5_register_loop(numpy, "conjugate", loop_splitcomplex5_unary, 2, 0, splitcomplex5_conjugate) ||
        pysplitcomplex5_register_loop(numpy, "reciprocal", loop_splitcomplex5_unary, 2, 0, splitcomplex5_inv) ||
        pysplitcomplex5_register_loop(numpy, "square", loop_splitcomplex5_unary, 2, 0, splitcomplex5_square) ||
        pysplitcomplex5_register_loop(numpy, "exp", loop_splitcomplex5_unary, 2, 0, splitcomplex5_exp) ||
        pysplitcomplex5_register_loop(numpy, "expm1", loop_splitcomplex5_unary, 2, 0, splitcomplex5_expm1) ||
        pysplitcomplex5_register_loop(numpy, "log", loop_splitcomplex5_unary, 2, 0, splitcomplex5_log) ||
        pysplitcomplex5_register_loop(numpy, "log1p", loop_splitcomplex5_unary, 2, 0, splitcomplex5_log1p) ||
        pysplitcomplex5_register_loop(numpy, "sqrt", loop_splitcomplex5_unary, 2, 0, splitcomplex5_sqrt) ||
        pysplitcomplex5_register_loop(numpy, "sin", loop_splitcomplex5_unary, 2, 0, splitcomplex5_sin) ||
        pysplitcomplex5_register_loop(numpy, "cos", loop_splitcomplex5_unary, 2, 0, splitcomplex5_cos) ||
        pysplitcomplex5_register_loop(numpy, "tan", loop_splitcomplex5_unary, 2, 0, splitcomplex5_tan) ||
        pysplitcomplex5_register_loop(numpy, "arcsin", loop_splitcomplex5_unary, 2, 0, splitcomplex5_asin) ||
        pysplitcomplex5_register_loop(numpy, "arccos", loop_splitcomplex5_unary, 2, 0, splitcomplex5_acos) ||
        pysplitcomplex5_register_loop(numpy, "arctan", loop_splitcomplex5_unary, 2, 0, splitcomplex5_atan) ||
        pysplitcomplex5_register_loop(numpy, "sinh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_sinh) ||
        pysplitcomplex5_register_loop(numpy, "cosh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_cosh) ||
        pysplitcomplex5_register_loop(numpy, "tanh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_tanh) ||
        pysplitcomplex5_register_loop(numpy, "arcsinh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_asinh) ||
        pysplitcomplex5_register_loop(numpy, "arccosh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_acosh) ||
        pysplitcomplex5_register_loop(numpy, "arctanh", loop_splitcomplex5_unary, 2, 0, splitcomplex5_atanh) ||
        pysplitcomplex5_register_loop(numpy, "equal", loop_splitcomplex5_predicate, 3, 1, splitcomplex5_equal) ||
        pysplitcomplex5_register_loop(numpy, "not_equal", loop_splitcomplex5_predicate, 3, 1, splitcomplex5_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_splitcomplex5", "splitcomplex5 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__splitcomplex5(void)
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
    PySplitcomplex5_Type = PyType_FromModuleAndSpec(
        m, &pysplitcomplex5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySplitcomplex5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysplitcomplex5_arrfuncs);
    pysplitcomplex5_arrfuncs.getitem = npysplitcomplex5_getitem;
    pysplitcomplex5_arrfuncs.setitem = npysplitcomplex5_setitem;
    pysplitcomplex5_arrfuncs.copyswapn = pysplitcomplex5_copyswapn;
    pysplitcomplex5_arrfuncs.copyswap = pysplitcomplex5_copyswap;
    pysplitcomplex5_arrfuncs.nonzero = pysplitcomplex5_nonzero;
    pysplitcomplex5_descr_proto.typeobj = (PyTypeObject *)PySplitcomplex5_Type;
    Py_SET_TYPE((PyObject *)&pysplitcomplex5_descr_proto, &PyArrayDescr_Type);
    npy_splitcomplex5 = PyArray_RegisterDataType(&pysplitcomplex5_descr_proto);
    if (npy_splitcomplex5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_splitcomplex5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySplitcomplex5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PySplitcomplex5_Type, 32) < 0) {
        goto fail;
    }
    if (pysplitcomplex5_register_casts() < 0 ||
            pysplitcomplex5_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "splitcomplex5", PySplitcomplex5_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
