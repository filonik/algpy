/*
 * complex4 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex4.h.
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

#include "complex4.h"

typedef struct {
    PyObject_HEAD
    complex4 v;
} PyComplex4;

static PyObject *PyComplex4_Type = NULL;
static int npy_complex4 = -1;

static PyObject *
PyComplex4_FromValue(PyTypeObject *type, complex4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplex4 *p = (PyComplex4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX4(v) \
    PyComplex4_FromValue((PyTypeObject *)PyComplex4_Type, (v))

static int
complex4_from_object(PyObject *o, complex4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplex4_Type)) {
        *out = ((PyComplex4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PyComplex4_FromValue(type, complex4_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]));
}

static PyObject *
pycomplex4_repr(PyObject *self)
{
    complex4 v = ((PyComplex4 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex4_hash(PyObject *self)
{
    complex4 v = ((PyComplex4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_COMPLEX4_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (complex4_is_real(v)) {
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
pycomplex4_richcompare(PyObject *a, PyObject *b, int op)
{
    complex4 x, y;
    int r;
    if ((r = complex4_from_object(a, &x)) <= 0 ||
            (r = complex4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(complex4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!complex4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define COMPLEX4_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex4 x, y;                                                       \
        int r;                                                               \
        if ((r = complex4_from_object(a, &x)) <= 0 ||                        \
                (r = complex4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX4(expr(x, y));                                     \
    }

COMPLEX4_BINARY(add, complex4_add)
COMPLEX4_BINARY(sub, complex4_sub)
COMPLEX4_BINARY(mul, complex4_mul)
COMPLEX4_BINARY(div, complex4_div)

static PyObject *
pycomplex4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex4 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex4 does not support modular power");
        return NULL;
    }
    if ((r = complex4_from_object(a, &x)) <= 0 ||
            (r = complex4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX4(complex4_rpow(x, y));
}

#define COMPLEX4_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX4(expr(((PyComplex4 *)self)->v));                  \
    }

COMPLEX4_UNARY(pos, complex4_pos)
COMPLEX4_UNARY(neg, complex4_neg)

static PyObject *
pycomplex4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue((PyTypeObject *)cls, complex4_zero());
}

static PyObject *
pycomplex4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue((PyTypeObject *)cls, complex4_one());
}

static int
pycomplex4_bool(PyObject *self)
{
    complex4 v = ((PyComplex4 *)self)->v;
    return vec16_any_nonzero(v);
}

#define COMPLEX4_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX4(op(((PyComplex4 *)self)->v));                   \
    }

COMPLEX4_METHOD(conjugate, complex4_conjugate)
COMPLEX4_METHOD(exp, complex4_exp)
COMPLEX4_METHOD(log, complex4_log)
COMPLEX4_METHOD(sqrt, complex4_sqrt)

static PyObject *
pycomplex4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex4_norm(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX4(complex4_inv(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex4_lmat(((PyComplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex4_rmat(((PyComplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {17};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex4_minpoly(((PyComplex4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[16];
    npy_intp dims[1] = {16};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex4_eivals(((PyComplex4 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 16; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex4);
    complex4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex4 *)PyArray_DATA(arr);
    complex4_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex4_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex4);
    complex4 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex4 *)PyArray_DATA(arr);
    complex4_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex4 v = ((PyComplex4 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pycomplex4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_square(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_expm1(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_log1p(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_sin(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_cos(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_tan(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_asin(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_acos(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_atan(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_sinh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_cosh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_tanh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_asinh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_acosh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_atanh(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_re(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex4_FromValue(Py_TYPE(self), complex4_im(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_real(PyObject *self, void *closure)
{
    return pycomplex4_re(self, NULL);
}

static PyObject *
pycomplex4_imag(PyObject *self, void *closure)
{
    return pycomplex4_im(self, NULL);
}

static PyGetSetDef pycomplex4_getset[] = {
    {"real", pycomplex4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex4_is_real(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex4_to_re(((PyComplex4 *)self)->v));
}

static PyObject *
pycomplex4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplex4_FromValue((PyTypeObject *)cls, complex4_from_re(a));
}

static PyObject *
pycomplex4_lpow(PyObject *self, PyObject *other)
{
    complex4 b;
    int r;
    r = complex4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex4_FromValue(Py_TYPE(self),
        complex4_lpow(((PyComplex4 *)self)->v, b));
}

static PyObject *
pycomplex4_rpow(PyObject *self, PyObject *other)
{
    complex4 b;
    int r;
    r = complex4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex4_FromValue(Py_TYPE(self),
        complex4_rpow(((PyComplex4 *)self)->v, b));
}

static PyObject *
pycomplex4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex4 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex4_FromValue(Py_TYPE(self),
        complex4_mix(((PyComplex4 *)self)->v, b, t));
}

static PyObject *
pycomplex4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex4, 16);
}

static PyObject *
pycomplex4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex4, 16);
}

static PyObject *
pycomplex4_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_COMPLEX4(complex4_ipow(((PyComplex4 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pycomplex4_methods[] = {
    {"to_float", (PyCFunction)pycomplex4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycomplex4_square, METH_NOARGS, "square."},
    {"expm1", pycomplex4_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycomplex4_log1p, METH_NOARGS, "log1p."},
    {"sin", pycomplex4_sin, METH_NOARGS, "sin."},
    {"cos", pycomplex4_cos, METH_NOARGS, "cos."},
    {"tan", pycomplex4_tan, METH_NOARGS, "tan."},
    {"asin", pycomplex4_asin, METH_NOARGS, "asin."},
    {"acos", pycomplex4_acos, METH_NOARGS, "acos."},
    {"atan", pycomplex4_atan, METH_NOARGS, "atan."},
    {"sinh", pycomplex4_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycomplex4_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycomplex4_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycomplex4_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycomplex4_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycomplex4_atanh, METH_NOARGS, "atanh."},
    {"re", pycomplex4_re, METH_NOARGS, "re."},
    {"im", pycomplex4_im, METH_NOARGS, "im."},
    {"arcsin", pycomplex4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycomplex4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycomplex4_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pycomplex4_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycomplex4_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycomplex4_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycomplex4_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex4_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pycomplex4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex4_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pycomplex4_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex4_exp, METH_NOARGS, "exp."},
    {"log", pycomplex4_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex4_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex4_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pycomplex4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplex4 *)self)->v.data[i]);
}

static PyObject *
pycomplex4_subscript(PyObject *self, PyObject *key)
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
    return pycomplex4_item(self, i);
}

static int
pycomplex4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex4_slots[] = {
    {Py_tp_dealloc, pycomplex4_dealloc},
    {Py_tp_traverse, pycomplex4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex4_repr},
    {Py_tp_str, pycomplex4_repr},
    {Py_tp_hash, pycomplex4_hash},
    {Py_tp_doc, (void *)"complex4 number with float64 components"},
    {Py_tp_richcompare, pycomplex4_richcompare},
    {Py_tp_methods, pycomplex4_methods},
    {Py_tp_getset, pycomplex4_getset},
    {Py_sq_length, pycomplex4_length},
    {Py_sq_item, pycomplex4_item},
    {Py_mp_subscript, pycomplex4_subscript},
    {Py_tp_new, pycomplex4_new},
    {Py_nb_add, pycomplex4_add},
    {Py_nb_subtract, pycomplex4_sub},
    {Py_nb_multiply, pycomplex4_mul},
    {Py_nb_true_divide, pycomplex4_div},
    {Py_nb_power, (void *)pycomplex4_pow},
    {Py_nb_positive, pycomplex4_pos},
    {Py_nb_negative, pycomplex4_neg},
    {Py_nb_bool, pycomplex4_bool},
    {0, NULL}
};

static PyType_Spec pycomplex4_spec = {
    .name = "algpy.algebras.multicomplex.complex4",
    .basicsize = sizeof(PyComplex4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex4_slots,
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
npycomplex4_getitem(void *data, void *arr)
{
    complex4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_COMPLEX4(v);
}

static int
npycomplex4_setitem(PyObject *item, void *data, void *arr)
{
    complex4 v;
    int r = complex4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex4 needs exactly 16 components");
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
pycomplex4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycomplex4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex4_copyswapn(dst, sizeof(complex4), src, sizeof(complex4), 1,
                         swap, arr);
}

static npy_bool
pycomplex4_nonzero(void *data, void *arr)
{
    complex4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex4_arrfuncs;
typedef struct { char c; complex4 v; } pycomplex4_align_test;

static PyArray_DescrProto pycomplex4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex4), offsetof(pycomplex4_align_test, v), 0, 0, 0,
    &pycomplex4_arrfuncs,
};

#define COMPLEX4_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycomplex4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex4 *dst = (complex4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex4_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define COMPLEX4_REAL_TYPES(X) \
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

COMPLEX4_REAL_TYPES(COMPLEX4_CAST_FROM)

static int
pycomplex4_register_casts(void)
{
#define REGISTER_COMPLEX4_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex4,                    \
                                     pycomplex4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX4_REAL_TYPES(REGISTER_COMPLEX4_CAST)
#undef REGISTER_COMPLEX4_CAST
    return 0;
}

typedef complex4 (*complex4_unary_op)(complex4);
typedef complex4 (*complex4_binary_op)(complex4, complex4);
typedef int (*complex4_predicate)(complex4, complex4);

static void
loop_complex4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex4_unary_op op = (complex4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex4_binary_op op = (complex4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex4 a, b, r;
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
loop_complex4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex4_predicate op = (complex4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex4;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex4_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex4_register_loop(numpy, "add", loop_complex4_binary, 3, 0, complex4_add) ||
        pycomplex4_register_loop(numpy, "subtract", loop_complex4_binary, 3, 0, complex4_sub) ||
        pycomplex4_register_loop(numpy, "multiply", loop_complex4_binary, 3, 0, complex4_mul) ||
        pycomplex4_register_loop(numpy, "divide", loop_complex4_binary, 3, 0, complex4_div) ||
        pycomplex4_register_loop(numpy, "power", loop_complex4_binary, 3, 0, complex4_rpow) ||
        pycomplex4_register_loop(numpy, "positive", loop_complex4_unary, 2, 0, complex4_pos) ||
        pycomplex4_register_loop(numpy, "negative", loop_complex4_unary, 2, 0, complex4_neg) ||
        pycomplex4_register_loop(numpy, "conjugate", loop_complex4_unary, 2, 0, complex4_conjugate) ||
        pycomplex4_register_loop(numpy, "reciprocal", loop_complex4_unary, 2, 0, complex4_inv) ||
        pycomplex4_register_loop(numpy, "square", loop_complex4_unary, 2, 0, complex4_square) ||
        pycomplex4_register_loop(numpy, "exp", loop_complex4_unary, 2, 0, complex4_exp) ||
        pycomplex4_register_loop(numpy, "expm1", loop_complex4_unary, 2, 0, complex4_expm1) ||
        pycomplex4_register_loop(numpy, "log", loop_complex4_unary, 2, 0, complex4_log) ||
        pycomplex4_register_loop(numpy, "log1p", loop_complex4_unary, 2, 0, complex4_log1p) ||
        pycomplex4_register_loop(numpy, "sqrt", loop_complex4_unary, 2, 0, complex4_sqrt) ||
        pycomplex4_register_loop(numpy, "sin", loop_complex4_unary, 2, 0, complex4_sin) ||
        pycomplex4_register_loop(numpy, "cos", loop_complex4_unary, 2, 0, complex4_cos) ||
        pycomplex4_register_loop(numpy, "tan", loop_complex4_unary, 2, 0, complex4_tan) ||
        pycomplex4_register_loop(numpy, "arcsin", loop_complex4_unary, 2, 0, complex4_asin) ||
        pycomplex4_register_loop(numpy, "arccos", loop_complex4_unary, 2, 0, complex4_acos) ||
        pycomplex4_register_loop(numpy, "arctan", loop_complex4_unary, 2, 0, complex4_atan) ||
        pycomplex4_register_loop(numpy, "sinh", loop_complex4_unary, 2, 0, complex4_sinh) ||
        pycomplex4_register_loop(numpy, "cosh", loop_complex4_unary, 2, 0, complex4_cosh) ||
        pycomplex4_register_loop(numpy, "tanh", loop_complex4_unary, 2, 0, complex4_tanh) ||
        pycomplex4_register_loop(numpy, "arcsinh", loop_complex4_unary, 2, 0, complex4_asinh) ||
        pycomplex4_register_loop(numpy, "arccosh", loop_complex4_unary, 2, 0, complex4_acosh) ||
        pycomplex4_register_loop(numpy, "arctanh", loop_complex4_unary, 2, 0, complex4_atanh) ||
        pycomplex4_register_loop(numpy, "equal", loop_complex4_predicate, 3, 1, complex4_equal) ||
        pycomplex4_register_loop(numpy, "not_equal", loop_complex4_predicate, 3, 1, complex4_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex4", "complex4 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex4(void)
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
    PyComplex4_Type = PyType_FromModuleAndSpec(
        m, &pycomplex4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplex4_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex4_arrfuncs);
    pycomplex4_arrfuncs.getitem = npycomplex4_getitem;
    pycomplex4_arrfuncs.setitem = npycomplex4_setitem;
    pycomplex4_arrfuncs.copyswapn = pycomplex4_copyswapn;
    pycomplex4_arrfuncs.copyswap = pycomplex4_copyswap;
    pycomplex4_arrfuncs.nonzero = pycomplex4_nonzero;
    pycomplex4_descr_proto.typeobj = (PyTypeObject *)PyComplex4_Type;
    Py_SET_TYPE((PyObject *)&pycomplex4_descr_proto, &PyArrayDescr_Type);
    npy_complex4 = PyArray_RegisterDataType(&pycomplex4_descr_proto);
    if (npy_complex4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplex4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplex4_Type, 16) < 0) {
        goto fail;
    }
    if (pycomplex4_register_casts() < 0 ||
            pycomplex4_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex4", PyComplex4_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
