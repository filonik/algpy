/*
 * trigintaduonion numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in trigintaduonion.h.
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

#include "trigintaduonion.h"

typedef struct {
    PyObject_HEAD
    trigintaduonion v;
} PyTrigintaduonion;

static PyObject *PyTrigintaduonion_Type = NULL;
static int npy_trigintaduonion = -1;

static PyObject *
PyTrigintaduonion_FromValue(PyTypeObject *type, trigintaduonion v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyTrigintaduonion *p = (PyTrigintaduonion *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_TRIGINTADUONION(v) \
    PyTrigintaduonion_FromValue((PyTypeObject *)PyTrigintaduonion_Type, (v))

static int
trigintaduonion_from_object(PyObject *o, trigintaduonion *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyTrigintaduonion_Type)) {
        *out = ((PyTrigintaduonion *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = trigintaduonion_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pytrigintaduonion_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyTrigintaduonion_FromValue(type, trigintaduonion_make(d));
}

static PyObject *
pytrigintaduonion_repr(PyObject *self)
{
    trigintaduonion v = ((PyTrigintaduonion *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("trigintaduonion(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pytrigintaduonion_hash(PyObject *self)
{
    trigintaduonion v = ((PyTrigintaduonion *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_TRIGINTADUONION_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (trigintaduonion_is_real(v)) {
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
pytrigintaduonion_richcompare(PyObject *a, PyObject *b, int op)
{
    trigintaduonion x, y;
    int r;
    if ((r = trigintaduonion_from_object(a, &x)) <= 0 ||
            (r = trigintaduonion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(trigintaduonion_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!trigintaduonion_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define TRIGINTADUONION_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pytrigintaduonion_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        trigintaduonion x, y;                                                       \
        int r;                                                               \
        if ((r = trigintaduonion_from_object(a, &x)) <= 0 ||                        \
                (r = trigintaduonion_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_TRIGINTADUONION(expr(x, y));                                     \
    }

TRIGINTADUONION_BINARY(add, trigintaduonion_add)
TRIGINTADUONION_BINARY(sub, trigintaduonion_sub)
TRIGINTADUONION_BINARY(mul, trigintaduonion_mul)
TRIGINTADUONION_BINARY(div, trigintaduonion_div)

static PyObject *
pytrigintaduonion_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    trigintaduonion x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "trigintaduonion does not support modular power");
        return NULL;
    }
    if ((r = trigintaduonion_from_object(a, &x)) <= 0 ||
            (r = trigintaduonion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_TRIGINTADUONION(trigintaduonion_rpow(x, y));
}

#define TRIGINTADUONION_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pytrigintaduonion_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_TRIGINTADUONION(expr(((PyTrigintaduonion *)self)->v));                  \
    }

TRIGINTADUONION_UNARY(pos, trigintaduonion_pos)
TRIGINTADUONION_UNARY(neg, trigintaduonion_neg)

static PyObject *
pytrigintaduonion_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue((PyTypeObject *)cls, trigintaduonion_zero());
}

static PyObject *
pytrigintaduonion_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue((PyTypeObject *)cls, trigintaduonion_one());
}

static int
pytrigintaduonion_bool(PyObject *self)
{
    trigintaduonion v = ((PyTrigintaduonion *)self)->v;
    return vec32_any_nonzero(v);
}

#define TRIGINTADUONION_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pytrigintaduonion_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_TRIGINTADUONION(op(((PyTrigintaduonion *)self)->v));                   \
    }

TRIGINTADUONION_METHOD(conjugate, trigintaduonion_conjugate)
TRIGINTADUONION_METHOD(exp, trigintaduonion_exp)
TRIGINTADUONION_METHOD(log, trigintaduonion_log)
TRIGINTADUONION_METHOD(sqrt, trigintaduonion_sqrt)

static PyObject *
pytrigintaduonion_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(trigintaduonion_norm(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_TRIGINTADUONION(trigintaduonion_inv(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    trigintaduonion_lmat(((PyTrigintaduonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pytrigintaduonion_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    trigintaduonion_rmat(((PyTrigintaduonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pytrigintaduonion_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    trigintaduonion_minpoly(((PyTrigintaduonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pytrigintaduonion_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[2];
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    trigintaduonion_eivals(((PyTrigintaduonion *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pytrigintaduonion_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_trigintaduonion);
    trigintaduonion *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (trigintaduonion *)PyArray_DATA(arr);
    trigintaduonion_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pytrigintaduonion_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_trigintaduonion);
    trigintaduonion *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (trigintaduonion *)PyArray_DATA(arr);
    trigintaduonion_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pytrigintaduonion_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_trigintaduonion);
    trigintaduonion *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (trigintaduonion *)PyArray_DATA(arr);
    trigintaduonion_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pytrigintaduonion_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    trigintaduonion v = ((PyTrigintaduonion *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pytrigintaduonion_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_square(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_expm1(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_log1p(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_sin(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_cos(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_tan(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_asin(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_acos(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_atan(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_sinh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_cosh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_tanh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_asinh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_acosh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_atanh(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_re(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_im(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_real(PyObject *self, void *closure)
{
    return pytrigintaduonion_re(self, NULL);
}

static PyObject *
pytrigintaduonion_imag(PyObject *self, void *closure)
{
    return pytrigintaduonion_im(self, NULL);
}

static PyObject *pytrigintaduonion_grades(PyObject *self, void *closure);

static PyGetSetDef pytrigintaduonion_getset[] = {
    {"grade", pytrigintaduonion_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pytrigintaduonion_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pytrigintaduonion_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pytrigintaduonion_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(trigintaduonion_is_real(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(trigintaduonion_to_re(((PyTrigintaduonion *)self)->v));
}

static PyObject *
pytrigintaduonion_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyTrigintaduonion_FromValue((PyTypeObject *)cls, trigintaduonion_from_re(a));
}

static PyObject *
pytrigintaduonion_lpow(PyObject *self, PyObject *other)
{
    trigintaduonion b;
    int r;
    r = trigintaduonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyTrigintaduonion_FromValue(Py_TYPE(self),
        trigintaduonion_lpow(((PyTrigintaduonion *)self)->v, b));
}

static PyObject *
pytrigintaduonion_rpow(PyObject *self, PyObject *other)
{
    trigintaduonion b;
    int r;
    r = trigintaduonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyTrigintaduonion_FromValue(Py_TYPE(self),
        trigintaduonion_rpow(((PyTrigintaduonion *)self)->v, b));
}

static PyObject *
pytrigintaduonion_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    trigintaduonion b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = trigintaduonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyTrigintaduonion_FromValue(Py_TYPE(self),
        trigintaduonion_mix(((PyTrigintaduonion *)self)->v, b, t));
}

static PyObject *
pytrigintaduonion_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_trigintaduonion, 32);
}

static PyObject *
pytrigintaduonion_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_trigintaduonion, 32);
}

static PyObject *
pytrigintaduonion_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_TRIGINTADUONION(trigintaduonion_ipow(((PyTrigintaduonion *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pytrigintaduonion_grade(PyObject *self, Py_ssize_t k)
{
    trigintaduonion v = ((PyTrigintaduonion *)self)->v;
    switch (k) {
        case 0:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade0(v));
        case 1:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade1(v));
        case 2:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade2(v));
        case 3:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade3(v));
        case 4:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade4(v));
        case 5:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_grade5(v));
        default:
            return PyTrigintaduonion_FromValue(Py_TYPE(self), trigintaduonion_zero());
    }
}

static PyObject *
pytrigintaduonion_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pytrigintaduonion_grade, 6);
}

static PyMethodDef pytrigintaduonion_methods[] = {
    {"to_float", (PyCFunction)pytrigintaduonion_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pytrigintaduonion_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pytrigintaduonion_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pytrigintaduonion_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pytrigintaduonion_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pytrigintaduonion_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pytrigintaduonion_square, METH_NOARGS, "square."},
    {"expm1", pytrigintaduonion_expm1, METH_NOARGS, "expm1."},
    {"log1p", pytrigintaduonion_log1p, METH_NOARGS, "log1p."},
    {"sin", pytrigintaduonion_sin, METH_NOARGS, "sin."},
    {"cos", pytrigintaduonion_cos, METH_NOARGS, "cos."},
    {"tan", pytrigintaduonion_tan, METH_NOARGS, "tan."},
    {"asin", pytrigintaduonion_asin, METH_NOARGS, "asin."},
    {"acos", pytrigintaduonion_acos, METH_NOARGS, "acos."},
    {"atan", pytrigintaduonion_atan, METH_NOARGS, "atan."},
    {"sinh", pytrigintaduonion_sinh, METH_NOARGS, "sinh."},
    {"cosh", pytrigintaduonion_cosh, METH_NOARGS, "cosh."},
    {"tanh", pytrigintaduonion_tanh, METH_NOARGS, "tanh."},
    {"asinh", pytrigintaduonion_asinh, METH_NOARGS, "asinh."},
    {"acosh", pytrigintaduonion_acosh, METH_NOARGS, "acosh."},
    {"atanh", pytrigintaduonion_atanh, METH_NOARGS, "atanh."},
    {"re", pytrigintaduonion_re, METH_NOARGS, "re."},
    {"im", pytrigintaduonion_im, METH_NOARGS, "im."},
    {"arcsin", pytrigintaduonion_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pytrigintaduonion_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pytrigintaduonion_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pytrigintaduonion_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pytrigintaduonion_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pytrigintaduonion_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pytrigintaduonion_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pytrigintaduonion_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pytrigintaduonion_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pytrigintaduonion_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pytrigintaduonion_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pytrigintaduonion_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pytrigintaduonion_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pytrigintaduonion_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pytrigintaduonion_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pytrigintaduonion_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pytrigintaduonion_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pytrigintaduonion_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pytrigintaduonion_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pytrigintaduonion_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pytrigintaduonion_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pytrigintaduonion_exp, METH_NOARGS, "exp."},
    {"log", pytrigintaduonion_log, METH_NOARGS, "log."},
    {"sqrt", pytrigintaduonion_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pytrigintaduonion_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pytrigintaduonion_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pytrigintaduonion_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyTrigintaduonion *)self)->v.data[i]);
}

static PyObject *
pytrigintaduonion_subscript(PyObject *self, PyObject *key)
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
    return pytrigintaduonion_item(self, i);
}

static int
pytrigintaduonion_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pytrigintaduonion_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pytrigintaduonion_slots[] = {
    {Py_tp_dealloc, pytrigintaduonion_dealloc},
    {Py_tp_traverse, pytrigintaduonion_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pytrigintaduonion_repr},
    {Py_tp_str, pytrigintaduonion_repr},
    {Py_tp_hash, pytrigintaduonion_hash},
    {Py_tp_doc, (void *)"Trigintaduonion with 32 float64 coefficients"},
    {Py_tp_richcompare, pytrigintaduonion_richcompare},
    {Py_tp_methods, pytrigintaduonion_methods},
    {Py_tp_getset, pytrigintaduonion_getset},
    {Py_sq_length, pytrigintaduonion_length},
    {Py_sq_item, pytrigintaduonion_item},
    {Py_mp_subscript, pytrigintaduonion_subscript},
    {Py_tp_new, pytrigintaduonion_new},
    {Py_nb_add, pytrigintaduonion_add},
    {Py_nb_subtract, pytrigintaduonion_sub},
    {Py_nb_multiply, pytrigintaduonion_mul},
    {Py_nb_true_divide, pytrigintaduonion_div},
    {Py_nb_power, (void *)pytrigintaduonion_pow},
    {Py_nb_positive, pytrigintaduonion_pos},
    {Py_nb_negative, pytrigintaduonion_neg},
    {Py_nb_bool, pytrigintaduonion_bool},
    {0, NULL}
};

static PyType_Spec pytrigintaduonion_spec = {
    .name = "algpy.algebras.cayleydickson.trigintaduonion",
    .basicsize = sizeof(PyTrigintaduonion),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pytrigintaduonion_slots,
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
npytrigintaduonion_getitem(void *data, void *arr)
{
    trigintaduonion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_TRIGINTADUONION(v);
}

static int
npytrigintaduonion_setitem(PyObject *item, void *data, void *arr)
{
    trigintaduonion v;
    int r = trigintaduonion_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to trigintaduonion");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "trigintaduonion needs exactly 32 components");
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
pytrigintaduonion_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(trigintaduonion));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pytrigintaduonion_copyswap(void *dst, void *src, int swap, void *arr)
{
    pytrigintaduonion_copyswapn(dst, sizeof(trigintaduonion), src, sizeof(trigintaduonion), 1,
                         swap, arr);
}

static npy_bool
pytrigintaduonion_nonzero(void *data, void *arr)
{
    trigintaduonion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pytrigintaduonion_arrfuncs;
typedef struct { char c; trigintaduonion v; } pytrigintaduonion_align_test;

static PyArray_DescrProto pytrigintaduonion_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(trigintaduonion), offsetof(pytrigintaduonion_align_test, v), 0, 0, 0,
    &pytrigintaduonion_arrfuncs,
};

#define TRIGINTADUONION_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pytrigintaduonion_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        trigintaduonion *dst = (trigintaduonion *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = trigintaduonion_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define TRIGINTADUONION_REAL_TYPES(X) \
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

TRIGINTADUONION_REAL_TYPES(TRIGINTADUONION_CAST_FROM)

static int
pytrigintaduonion_register_casts(void)
{
#define REGISTER_TRIGINTADUONION_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_trigintaduonion,                    \
                                     pytrigintaduonion_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_trigintaduonion, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    TRIGINTADUONION_REAL_TYPES(REGISTER_TRIGINTADUONION_CAST)
#undef REGISTER_TRIGINTADUONION_CAST
    return 0;
}

typedef trigintaduonion (*trigintaduonion_unary_op)(trigintaduonion);
typedef trigintaduonion (*trigintaduonion_binary_op)(trigintaduonion, trigintaduonion);
typedef int (*trigintaduonion_predicate)(trigintaduonion, trigintaduonion);

static void
loop_trigintaduonion_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    trigintaduonion_unary_op op = (trigintaduonion_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        trigintaduonion a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_trigintaduonion_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    trigintaduonion_binary_op op = (trigintaduonion_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        trigintaduonion a, b, r;
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
loop_trigintaduonion_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    trigintaduonion_predicate op = (trigintaduonion_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        trigintaduonion a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pytrigintaduonion_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_trigintaduonion;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_trigintaduonion,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pytrigintaduonion_register_ufuncs(PyObject *numpy)
{
    return
        pytrigintaduonion_register_loop(numpy, "add", loop_trigintaduonion_binary, 3, 0, trigintaduonion_add) ||
        pytrigintaduonion_register_loop(numpy, "subtract", loop_trigintaduonion_binary, 3, 0, trigintaduonion_sub) ||
        pytrigintaduonion_register_loop(numpy, "multiply", loop_trigintaduonion_binary, 3, 0, trigintaduonion_mul) ||
        pytrigintaduonion_register_loop(numpy, "divide", loop_trigintaduonion_binary, 3, 0, trigintaduonion_div) ||
        pytrigintaduonion_register_loop(numpy, "power", loop_trigintaduonion_binary, 3, 0, trigintaduonion_rpow) ||
        pytrigintaduonion_register_loop(numpy, "positive", loop_trigintaduonion_unary, 2, 0, trigintaduonion_pos) ||
        pytrigintaduonion_register_loop(numpy, "negative", loop_trigintaduonion_unary, 2, 0, trigintaduonion_neg) ||
        pytrigintaduonion_register_loop(numpy, "conjugate", loop_trigintaduonion_unary, 2, 0, trigintaduonion_conjugate) ||
        pytrigintaduonion_register_loop(numpy, "reciprocal", loop_trigintaduonion_unary, 2, 0, trigintaduonion_inv) ||
        pytrigintaduonion_register_loop(numpy, "square", loop_trigintaduonion_unary, 2, 0, trigintaduonion_square) ||
        pytrigintaduonion_register_loop(numpy, "exp", loop_trigintaduonion_unary, 2, 0, trigintaduonion_exp) ||
        pytrigintaduonion_register_loop(numpy, "expm1", loop_trigintaduonion_unary, 2, 0, trigintaduonion_expm1) ||
        pytrigintaduonion_register_loop(numpy, "log", loop_trigintaduonion_unary, 2, 0, trigintaduonion_log) ||
        pytrigintaduonion_register_loop(numpy, "log1p", loop_trigintaduonion_unary, 2, 0, trigintaduonion_log1p) ||
        pytrigintaduonion_register_loop(numpy, "sqrt", loop_trigintaduonion_unary, 2, 0, trigintaduonion_sqrt) ||
        pytrigintaduonion_register_loop(numpy, "sin", loop_trigintaduonion_unary, 2, 0, trigintaduonion_sin) ||
        pytrigintaduonion_register_loop(numpy, "cos", loop_trigintaduonion_unary, 2, 0, trigintaduonion_cos) ||
        pytrigintaduonion_register_loop(numpy, "tan", loop_trigintaduonion_unary, 2, 0, trigintaduonion_tan) ||
        pytrigintaduonion_register_loop(numpy, "arcsin", loop_trigintaduonion_unary, 2, 0, trigintaduonion_asin) ||
        pytrigintaduonion_register_loop(numpy, "arccos", loop_trigintaduonion_unary, 2, 0, trigintaduonion_acos) ||
        pytrigintaduonion_register_loop(numpy, "arctan", loop_trigintaduonion_unary, 2, 0, trigintaduonion_atan) ||
        pytrigintaduonion_register_loop(numpy, "sinh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_sinh) ||
        pytrigintaduonion_register_loop(numpy, "cosh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_cosh) ||
        pytrigintaduonion_register_loop(numpy, "tanh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_tanh) ||
        pytrigintaduonion_register_loop(numpy, "arcsinh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_asinh) ||
        pytrigintaduonion_register_loop(numpy, "arccosh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_acosh) ||
        pytrigintaduonion_register_loop(numpy, "arctanh", loop_trigintaduonion_unary, 2, 0, trigintaduonion_atanh) ||
        pytrigintaduonion_register_loop(numpy, "equal", loop_trigintaduonion_predicate, 3, 1, trigintaduonion_equal) ||
        pytrigintaduonion_register_loop(numpy, "not_equal", loop_trigintaduonion_predicate, 3, 1, trigintaduonion_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_trigintaduonion", "trigintaduonion number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__trigintaduonion(void)
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
    PyTrigintaduonion_Type = PyType_FromModuleAndSpec(
        m, &pytrigintaduonion_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyTrigintaduonion_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pytrigintaduonion_arrfuncs);
    pytrigintaduonion_arrfuncs.getitem = npytrigintaduonion_getitem;
    pytrigintaduonion_arrfuncs.setitem = npytrigintaduonion_setitem;
    pytrigintaduonion_arrfuncs.copyswapn = pytrigintaduonion_copyswapn;
    pytrigintaduonion_arrfuncs.copyswap = pytrigintaduonion_copyswap;
    pytrigintaduonion_arrfuncs.nonzero = pytrigintaduonion_nonzero;
    pytrigintaduonion_descr_proto.typeobj = (PyTypeObject *)PyTrigintaduonion_Type;
    Py_SET_TYPE((PyObject *)&pytrigintaduonion_descr_proto, &PyArrayDescr_Type);
    npy_trigintaduonion = PyArray_RegisterDataType(&pytrigintaduonion_descr_proto);
    if (npy_trigintaduonion < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_trigintaduonion);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyTrigintaduonion_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyTrigintaduonion_Type, 32) < 0) {
        goto fail;
    }
    if (pytrigintaduonion_register_casts() < 0 ||
            pytrigintaduonion_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "trigintaduonion", PyTrigintaduonion_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
