/*
 * sedenion numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in sedenion.h.
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

#include "sedenion.h"

typedef struct {
    PyObject_HEAD
    sedenion v;
} PySedenion;

static PyObject *PySedenion_Type = NULL;
static int npy_sedenion = -1;

static PyObject *
PySedenion_FromValue(PyTypeObject *type, sedenion v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySedenion *p = (PySedenion *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SEDENION(v) \
    PySedenion_FromValue((PyTypeObject *)PySedenion_Type, (v))

static int
sedenion_from_object(PyObject *o, sedenion *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySedenion_Type)) {
        *out = ((PySedenion *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = sedenion_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pysedenion_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PySedenion_FromValue(type, sedenion_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]));
}

static PyObject *
pysedenion_repr(PyObject *self)
{
    sedenion v = ((PySedenion *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("sedenion(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pysedenion_hash(PyObject *self)
{
    sedenion v = ((PySedenion *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_SEDENION_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (sedenion_is_real(v)) {
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
pysedenion_richcompare(PyObject *a, PyObject *b, int op)
{
    sedenion x, y;
    int r;
    if ((r = sedenion_from_object(a, &x)) <= 0 ||
            (r = sedenion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(sedenion_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!sedenion_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SEDENION_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysedenion_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        sedenion x, y;                                                       \
        int r;                                                               \
        if ((r = sedenion_from_object(a, &x)) <= 0 ||                        \
                (r = sedenion_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SEDENION(expr(x, y));                                     \
    }

SEDENION_BINARY(add, sedenion_add)
SEDENION_BINARY(sub, sedenion_sub)
SEDENION_BINARY(mul, sedenion_mul)
SEDENION_BINARY(div, sedenion_div)

static PyObject *
pysedenion_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    sedenion x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "sedenion does not support modular power");
        return NULL;
    }
    if ((r = sedenion_from_object(a, &x)) <= 0 ||
            (r = sedenion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SEDENION(sedenion_rpow(x, y));
}

#define SEDENION_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysedenion_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SEDENION(expr(((PySedenion *)self)->v));                  \
    }

SEDENION_UNARY(pos, sedenion_pos)
SEDENION_UNARY(neg, sedenion_neg)

static PyObject *
pysedenion_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue((PyTypeObject *)cls, sedenion_zero());
}

static PyObject *
pysedenion_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue((PyTypeObject *)cls, sedenion_one());
}

static int
pysedenion_bool(PyObject *self)
{
    sedenion v = ((PySedenion *)self)->v;
    return vec16_any_nonzero(v);
}

#define SEDENION_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysedenion_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SEDENION(op(((PySedenion *)self)->v));                   \
    }

SEDENION_METHOD(conjugate, sedenion_conjugate)
SEDENION_METHOD(exp, sedenion_exp)
SEDENION_METHOD(log, sedenion_log)
SEDENION_METHOD(sqrt, sedenion_sqrt)

static PyObject *
pysedenion_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(sedenion_norm(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_SEDENION(sedenion_inv(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    sedenion_lmat(((PySedenion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysedenion_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    sedenion_rmat(((PySedenion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysedenion_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    sedenion_minpoly(((PySedenion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysedenion_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[2];
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    sedenion_eivals(((PySedenion *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pysedenion_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_sedenion);
    sedenion *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (sedenion *)PyArray_DATA(arr);
    sedenion_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysedenion_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_sedenion);
    sedenion *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (sedenion *)PyArray_DATA(arr);
    sedenion_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysedenion_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    sedenion v = ((PySedenion *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pysedenion_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_square(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_expm1(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_log1p(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_sin(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_cos(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_tan(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_asin(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_acos(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_atan(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_sinh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_cosh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_tanh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_asinh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_acosh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_atanh(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_re(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PySedenion_FromValue(Py_TYPE(self), sedenion_im(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_real(PyObject *self, void *closure)
{
    return pysedenion_re(self, NULL);
}

static PyObject *
pysedenion_imag(PyObject *self, void *closure)
{
    return pysedenion_im(self, NULL);
}

static PyGetSetDef pysedenion_getset[] = {
    {"real", pysedenion_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pysedenion_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pysedenion_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(sedenion_is_real(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(sedenion_to_re(((PySedenion *)self)->v));
}

static PyObject *
pysedenion_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PySedenion_FromValue((PyTypeObject *)cls, sedenion_from_re(a));
}

static PyObject *
pysedenion_lpow(PyObject *self, PyObject *other)
{
    sedenion b;
    int r;
    r = sedenion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySedenion_FromValue(Py_TYPE(self),
        sedenion_lpow(((PySedenion *)self)->v, b));
}

static PyObject *
pysedenion_rpow(PyObject *self, PyObject *other)
{
    sedenion b;
    int r;
    r = sedenion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySedenion_FromValue(Py_TYPE(self),
        sedenion_rpow(((PySedenion *)self)->v, b));
}

static PyObject *
pysedenion_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    sedenion b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = sedenion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PySedenion_FromValue(Py_TYPE(self),
        sedenion_mix(((PySedenion *)self)->v, b, t));
}

static PyObject *
pysedenion_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_sedenion, 16);
}

static PyObject *
pysedenion_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_sedenion, 16);
}

static PyObject *
pysedenion_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_SEDENION(sedenion_ipow(((PySedenion *)self)->v, (int64_t)exponent));
}

static PyMethodDef pysedenion_methods[] = {
    {"to_float", (PyCFunction)pysedenion_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pysedenion_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pysedenion_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pysedenion_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pysedenion_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pysedenion_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pysedenion_square, METH_NOARGS, "square."},
    {"expm1", pysedenion_expm1, METH_NOARGS, "expm1."},
    {"log1p", pysedenion_log1p, METH_NOARGS, "log1p."},
    {"sin", pysedenion_sin, METH_NOARGS, "sin."},
    {"cos", pysedenion_cos, METH_NOARGS, "cos."},
    {"tan", pysedenion_tan, METH_NOARGS, "tan."},
    {"asin", pysedenion_asin, METH_NOARGS, "asin."},
    {"acos", pysedenion_acos, METH_NOARGS, "acos."},
    {"atan", pysedenion_atan, METH_NOARGS, "atan."},
    {"sinh", pysedenion_sinh, METH_NOARGS, "sinh."},
    {"cosh", pysedenion_cosh, METH_NOARGS, "cosh."},
    {"tanh", pysedenion_tanh, METH_NOARGS, "tanh."},
    {"asinh", pysedenion_asinh, METH_NOARGS, "asinh."},
    {"acosh", pysedenion_acosh, METH_NOARGS, "acosh."},
    {"atanh", pysedenion_atanh, METH_NOARGS, "atanh."},
    {"re", pysedenion_re, METH_NOARGS, "re."},
    {"im", pysedenion_im, METH_NOARGS, "im."},
    {"arcsin", pysedenion_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pysedenion_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pysedenion_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pysedenion_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pysedenion_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pysedenion_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pysedenion_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysedenion_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysedenion_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pysedenion_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pysedenion_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pysedenion_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pysedenion_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pysedenion_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pysedenion_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pysedenion_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pysedenion_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pysedenion_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pysedenion_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pysedenion_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysedenion_exp, METH_NOARGS, "exp."},
    {"log", pysedenion_log, METH_NOARGS, "log."},
    {"sqrt", pysedenion_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pysedenion_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysedenion_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pysedenion_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySedenion *)self)->v.data[i]);
}

static PyObject *
pysedenion_subscript(PyObject *self, PyObject *key)
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
    return pysedenion_item(self, i);
}

static int
pysedenion_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysedenion_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysedenion_slots[] = {
    {Py_tp_dealloc, pysedenion_dealloc},
    {Py_tp_traverse, pysedenion_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysedenion_repr},
    {Py_tp_str, pysedenion_repr},
    {Py_tp_hash, pysedenion_hash},
    {Py_tp_doc, (void *)"sedenion number with float64 components"},
    {Py_tp_richcompare, pysedenion_richcompare},
    {Py_tp_methods, pysedenion_methods},
    {Py_tp_getset, pysedenion_getset},
    {Py_sq_length, pysedenion_length},
    {Py_sq_item, pysedenion_item},
    {Py_mp_subscript, pysedenion_subscript},
    {Py_tp_new, pysedenion_new},
    {Py_nb_add, pysedenion_add},
    {Py_nb_subtract, pysedenion_sub},
    {Py_nb_multiply, pysedenion_mul},
    {Py_nb_true_divide, pysedenion_div},
    {Py_nb_power, (void *)pysedenion_pow},
    {Py_nb_positive, pysedenion_pos},
    {Py_nb_negative, pysedenion_neg},
    {Py_nb_bool, pysedenion_bool},
    {0, NULL}
};

static PyType_Spec pysedenion_spec = {
    .name = "algpy.algebras.cayleydickson.sedenion",
    .basicsize = sizeof(PySedenion),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysedenion_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npysedenion_getitem(void *data, void *arr)
{
    sedenion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_SEDENION(v);
}

static int
npysedenion_setitem(PyObject *item, void *data, void *arr)
{
    sedenion v;
    int r = sedenion_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to sedenion");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "sedenion needs exactly 16 components");
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
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pysedenion_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(sedenion));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pysedenion_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysedenion_copyswapn(dst, sizeof(sedenion), src, sizeof(sedenion), 1,
                         swap, arr);
}

static npy_bool
pysedenion_nonzero(void *data, void *arr)
{
    sedenion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysedenion_arrfuncs;
typedef struct { char c; sedenion v; } pysedenion_align_test;

static PyArray_DescrProto pysedenion_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(sedenion), offsetof(pysedenion_align_test, v), 0, 0, 0,
    &pysedenion_arrfuncs,
};

#define SEDENION_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pysedenion_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        sedenion *dst = (sedenion *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = sedenion_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define SEDENION_REAL_TYPES(X) \
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

SEDENION_REAL_TYPES(SEDENION_CAST_FROM)

static int
pysedenion_register_casts(void)
{
#define REGISTER_SEDENION_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_sedenion,                    \
                                     pysedenion_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_sedenion, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SEDENION_REAL_TYPES(REGISTER_SEDENION_CAST)
#undef REGISTER_SEDENION_CAST
    return 0;
}

typedef sedenion (*sedenion_unary_op)(sedenion);
typedef sedenion (*sedenion_binary_op)(sedenion, sedenion);
typedef int (*sedenion_predicate)(sedenion, sedenion);

static void
loop_sedenion_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    sedenion_unary_op op = (sedenion_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        sedenion a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_sedenion_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    sedenion_binary_op op = (sedenion_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        sedenion a, b, r;
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
loop_sedenion_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    sedenion_predicate op = (sedenion_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        sedenion a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysedenion_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_sedenion;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_sedenion,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysedenion_register_ufuncs(PyObject *numpy)
{
    return
        pysedenion_register_loop(numpy, "add", loop_sedenion_binary, 3, 0, sedenion_add) ||
        pysedenion_register_loop(numpy, "subtract", loop_sedenion_binary, 3, 0, sedenion_sub) ||
        pysedenion_register_loop(numpy, "multiply", loop_sedenion_binary, 3, 0, sedenion_mul) ||
        pysedenion_register_loop(numpy, "divide", loop_sedenion_binary, 3, 0, sedenion_div) ||
        pysedenion_register_loop(numpy, "power", loop_sedenion_binary, 3, 0, sedenion_rpow) ||
        pysedenion_register_loop(numpy, "positive", loop_sedenion_unary, 2, 0, sedenion_pos) ||
        pysedenion_register_loop(numpy, "negative", loop_sedenion_unary, 2, 0, sedenion_neg) ||
        pysedenion_register_loop(numpy, "conjugate", loop_sedenion_unary, 2, 0, sedenion_conjugate) ||
        pysedenion_register_loop(numpy, "reciprocal", loop_sedenion_unary, 2, 0, sedenion_inv) ||
        pysedenion_register_loop(numpy, "square", loop_sedenion_unary, 2, 0, sedenion_square) ||
        pysedenion_register_loop(numpy, "exp", loop_sedenion_unary, 2, 0, sedenion_exp) ||
        pysedenion_register_loop(numpy, "expm1", loop_sedenion_unary, 2, 0, sedenion_expm1) ||
        pysedenion_register_loop(numpy, "log", loop_sedenion_unary, 2, 0, sedenion_log) ||
        pysedenion_register_loop(numpy, "log1p", loop_sedenion_unary, 2, 0, sedenion_log1p) ||
        pysedenion_register_loop(numpy, "sqrt", loop_sedenion_unary, 2, 0, sedenion_sqrt) ||
        pysedenion_register_loop(numpy, "sin", loop_sedenion_unary, 2, 0, sedenion_sin) ||
        pysedenion_register_loop(numpy, "cos", loop_sedenion_unary, 2, 0, sedenion_cos) ||
        pysedenion_register_loop(numpy, "tan", loop_sedenion_unary, 2, 0, sedenion_tan) ||
        pysedenion_register_loop(numpy, "arcsin", loop_sedenion_unary, 2, 0, sedenion_asin) ||
        pysedenion_register_loop(numpy, "arccos", loop_sedenion_unary, 2, 0, sedenion_acos) ||
        pysedenion_register_loop(numpy, "arctan", loop_sedenion_unary, 2, 0, sedenion_atan) ||
        pysedenion_register_loop(numpy, "sinh", loop_sedenion_unary, 2, 0, sedenion_sinh) ||
        pysedenion_register_loop(numpy, "cosh", loop_sedenion_unary, 2, 0, sedenion_cosh) ||
        pysedenion_register_loop(numpy, "tanh", loop_sedenion_unary, 2, 0, sedenion_tanh) ||
        pysedenion_register_loop(numpy, "arcsinh", loop_sedenion_unary, 2, 0, sedenion_asinh) ||
        pysedenion_register_loop(numpy, "arccosh", loop_sedenion_unary, 2, 0, sedenion_acosh) ||
        pysedenion_register_loop(numpy, "arctanh", loop_sedenion_unary, 2, 0, sedenion_atanh) ||
        pysedenion_register_loop(numpy, "equal", loop_sedenion_predicate, 3, 1, sedenion_equal) ||
        pysedenion_register_loop(numpy, "not_equal", loop_sedenion_predicate, 3, 1, sedenion_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_sedenion", "sedenion number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__sedenion(void)
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
    PySedenion_Type = PyType_FromModuleAndSpec(
        m, &pysedenion_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySedenion_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysedenion_arrfuncs);
    pysedenion_arrfuncs.getitem = npysedenion_getitem;
    pysedenion_arrfuncs.setitem = npysedenion_setitem;
    pysedenion_arrfuncs.copyswapn = pysedenion_copyswapn;
    pysedenion_arrfuncs.copyswap = pysedenion_copyswap;
    pysedenion_arrfuncs.nonzero = pysedenion_nonzero;
    pysedenion_descr_proto.typeobj = (PyTypeObject *)PySedenion_Type;
    Py_SET_TYPE((PyObject *)&pysedenion_descr_proto, &PyArrayDescr_Type);
    npy_sedenion = PyArray_RegisterDataType(&pysedenion_descr_proto);
    if (npy_sedenion < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_sedenion);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySedenion_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PySedenion_Type, 16) < 0) {
        goto fail;
    }
    if (pysedenion_register_casts() < 0 ||
            pysedenion_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "sedenion", PySedenion_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
