/*
 * octonion numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in octonion.h.
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

#include "octonion.h"

typedef struct {
    PyObject_HEAD
    octonion v;
} PyOctonion;

static PyObject *PyOctonion_Type = NULL;
static int npy_octonion = -1;

static PyObject *
PyOctonion_FromValue(PyTypeObject *type, octonion v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyOctonion *p = (PyOctonion *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_OCTONION(v) \
    PyOctonion_FromValue((PyTypeObject *)PyOctonion_Type, (v))

static int
octonion_from_object(PyObject *o, octonion *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyOctonion_Type)) {
        *out = ((PyOctonion *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = octonion_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyoctonion_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyOctonion_FromValue(type, octonion_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]));
}

static PyObject *
pyoctonion_repr(PyObject *self)
{
    octonion v = ((PyOctonion *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("octonion(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyoctonion_hash(PyObject *self)
{
    octonion v = ((PyOctonion *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_OCTONION_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (octonion_is_real(v)) {
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
pyoctonion_richcompare(PyObject *a, PyObject *b, int op)
{
    octonion x, y;
    int r;
    if ((r = octonion_from_object(a, &x)) <= 0 ||
            (r = octonion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(octonion_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!octonion_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define OCTONION_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyoctonion_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        octonion x, y;                                                       \
        int r;                                                               \
        if ((r = octonion_from_object(a, &x)) <= 0 ||                        \
                (r = octonion_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_OCTONION(expr(x, y));                                     \
    }

OCTONION_BINARY(add, octonion_add)
OCTONION_BINARY(sub, octonion_sub)
OCTONION_BINARY(mul, octonion_mul)
OCTONION_BINARY(div, octonion_div)

static PyObject *
pyoctonion_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    octonion x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "octonion does not support modular power");
        return NULL;
    }
    if ((r = octonion_from_object(a, &x)) <= 0 ||
            (r = octonion_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_OCTONION(octonion_rpow(x, y));
}

#define OCTONION_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyoctonion_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_OCTONION(expr(((PyOctonion *)self)->v));                  \
    }

OCTONION_UNARY(pos, octonion_pos)
OCTONION_UNARY(neg, octonion_neg)

static PyObject *
pyoctonion_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue((PyTypeObject *)cls, octonion_zero());
}

static PyObject *
pyoctonion_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue((PyTypeObject *)cls, octonion_one());
}

static int
pyoctonion_bool(PyObject *self)
{
    octonion v = ((PyOctonion *)self)->v;
    return vec8_any_nonzero(v);
}

#define OCTONION_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyoctonion_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_OCTONION(op(((PyOctonion *)self)->v));                   \
    }

OCTONION_METHOD(conjugate, octonion_conjugate)
OCTONION_METHOD(exp, octonion_exp)
OCTONION_METHOD(log, octonion_log)
OCTONION_METHOD(sqrt, octonion_sqrt)

static PyObject *
pyoctonion_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(octonion_norm(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_OCTONION(octonion_inv(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    octonion_lmat(((PyOctonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyoctonion_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    octonion_rmat(((PyOctonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyoctonion_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    octonion_minpoly(((PyOctonion *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyoctonion_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[2];
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    octonion_eivals(((PyOctonion *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyoctonion_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_octonion);
    octonion *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (octonion *)PyArray_DATA(arr);
    octonion_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyoctonion_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_octonion);
    octonion *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (octonion *)PyArray_DATA(arr);
    octonion_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyoctonion_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    octonion v = ((PyOctonion *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pyoctonion_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_square(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_expm1(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_log1p(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_sin(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_cos(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_tan(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_asin(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_acos(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_atan(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_sinh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_cosh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_tanh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_asinh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_acosh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_atanh(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_re(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyOctonion_FromValue(Py_TYPE(self), octonion_im(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_real(PyObject *self, void *closure)
{
    return pyoctonion_re(self, NULL);
}

static PyObject *
pyoctonion_imag(PyObject *self, void *closure)
{
    return pyoctonion_im(self, NULL);
}

static PyGetSetDef pyoctonion_getset[] = {
    {"real", pyoctonion_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyoctonion_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyoctonion_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(octonion_is_real(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(octonion_to_re(((PyOctonion *)self)->v));
}

static PyObject *
pyoctonion_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyOctonion_FromValue((PyTypeObject *)cls, octonion_from_re(a));
}

static PyObject *
pyoctonion_lpow(PyObject *self, PyObject *other)
{
    octonion b;
    int r;
    r = octonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyOctonion_FromValue(Py_TYPE(self),
        octonion_lpow(((PyOctonion *)self)->v, b));
}

static PyObject *
pyoctonion_rpow(PyObject *self, PyObject *other)
{
    octonion b;
    int r;
    r = octonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyOctonion_FromValue(Py_TYPE(self),
        octonion_rpow(((PyOctonion *)self)->v, b));
}

static PyObject *
pyoctonion_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    octonion b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = octonion_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyOctonion_FromValue(Py_TYPE(self),
        octonion_mix(((PyOctonion *)self)->v, b, t));
}

static PyObject *
pyoctonion_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_octonion, 8);
}

static PyObject *
pyoctonion_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_octonion, 8);
}

static PyObject *
pyoctonion_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_OCTONION(octonion_ipow(((PyOctonion *)self)->v, (int64_t)exponent));
}

static PyMethodDef pyoctonion_methods[] = {
    {"to_float", (PyCFunction)pyoctonion_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyoctonion_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyoctonion_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyoctonion_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyoctonion_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyoctonion_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyoctonion_square, METH_NOARGS, "square."},
    {"expm1", pyoctonion_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyoctonion_log1p, METH_NOARGS, "log1p."},
    {"sin", pyoctonion_sin, METH_NOARGS, "sin."},
    {"cos", pyoctonion_cos, METH_NOARGS, "cos."},
    {"tan", pyoctonion_tan, METH_NOARGS, "tan."},
    {"asin", pyoctonion_asin, METH_NOARGS, "asin."},
    {"acos", pyoctonion_acos, METH_NOARGS, "acos."},
    {"atan", pyoctonion_atan, METH_NOARGS, "atan."},
    {"sinh", pyoctonion_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyoctonion_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyoctonion_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyoctonion_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyoctonion_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyoctonion_atanh, METH_NOARGS, "atanh."},
    {"re", pyoctonion_re, METH_NOARGS, "re."},
    {"im", pyoctonion_im, METH_NOARGS, "im."},
    {"arcsin", pyoctonion_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyoctonion_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyoctonion_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyoctonion_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyoctonion_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyoctonion_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyoctonion_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyoctonion_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyoctonion_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pyoctonion_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyoctonion_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyoctonion_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyoctonion_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyoctonion_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyoctonion_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyoctonion_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pyoctonion_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyoctonion_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyoctonion_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyoctonion_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyoctonion_exp, METH_NOARGS, "exp."},
    {"log", pyoctonion_log, METH_NOARGS, "log."},
    {"sqrt", pyoctonion_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyoctonion_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyoctonion_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pyoctonion_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyOctonion *)self)->v.data[i]);
}

static PyObject *
pyoctonion_subscript(PyObject *self, PyObject *key)
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
    return pyoctonion_item(self, i);
}

static int
pyoctonion_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyoctonion_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyoctonion_slots[] = {
    {Py_tp_dealloc, pyoctonion_dealloc},
    {Py_tp_traverse, pyoctonion_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyoctonion_repr},
    {Py_tp_str, pyoctonion_repr},
    {Py_tp_hash, pyoctonion_hash},
    {Py_tp_doc, (void *)"octonion number with float64 components"},
    {Py_tp_richcompare, pyoctonion_richcompare},
    {Py_tp_methods, pyoctonion_methods},
    {Py_tp_getset, pyoctonion_getset},
    {Py_sq_length, pyoctonion_length},
    {Py_sq_item, pyoctonion_item},
    {Py_mp_subscript, pyoctonion_subscript},
    {Py_tp_new, pyoctonion_new},
    {Py_nb_add, pyoctonion_add},
    {Py_nb_subtract, pyoctonion_sub},
    {Py_nb_multiply, pyoctonion_mul},
    {Py_nb_true_divide, pyoctonion_div},
    {Py_nb_power, (void *)pyoctonion_pow},
    {Py_nb_positive, pyoctonion_pos},
    {Py_nb_negative, pyoctonion_neg},
    {Py_nb_bool, pyoctonion_bool},
    {0, NULL}
};

static PyType_Spec pyoctonion_spec = {
    .name = "algpy.algebras.cayleydickson.octonion",
    .basicsize = sizeof(PyOctonion),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyoctonion_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyoctonion_getitem(void *data, void *arr)
{
    octonion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_OCTONION(v);
}

static int
npyoctonion_setitem(PyObject *item, void *data, void *arr)
{
    octonion v;
    int r = octonion_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to octonion");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "octonion needs exactly 8 components");
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
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pyoctonion_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(octonion));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyoctonion_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyoctonion_copyswapn(dst, sizeof(octonion), src, sizeof(octonion), 1,
                         swap, arr);
}

static npy_bool
pyoctonion_nonzero(void *data, void *arr)
{
    octonion v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyoctonion_arrfuncs;
typedef struct { char c; octonion v; } pyoctonion_align_test;

static PyArray_DescrProto pyoctonion_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(octonion), offsetof(pyoctonion_align_test, v), 0, 0, 0,
    &pyoctonion_arrfuncs,
};

#define OCTONION_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyoctonion_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        octonion *dst = (octonion *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = octonion_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define OCTONION_REAL_TYPES(X) \
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

OCTONION_REAL_TYPES(OCTONION_CAST_FROM)

static int
pyoctonion_register_casts(void)
{
#define REGISTER_OCTONION_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_octonion,                    \
                                     pyoctonion_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_octonion, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    OCTONION_REAL_TYPES(REGISTER_OCTONION_CAST)
#undef REGISTER_OCTONION_CAST
    return 0;
}

typedef octonion (*octonion_unary_op)(octonion);
typedef octonion (*octonion_binary_op)(octonion, octonion);
typedef int (*octonion_predicate)(octonion, octonion);

static void
loop_octonion_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    octonion_unary_op op = (octonion_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        octonion a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_octonion_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    octonion_binary_op op = (octonion_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        octonion a, b, r;
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
loop_octonion_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    octonion_predicate op = (octonion_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        octonion a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyoctonion_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_octonion;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_octonion,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyoctonion_register_ufuncs(PyObject *numpy)
{
    return
        pyoctonion_register_loop(numpy, "add", loop_octonion_binary, 3, 0, octonion_add) ||
        pyoctonion_register_loop(numpy, "subtract", loop_octonion_binary, 3, 0, octonion_sub) ||
        pyoctonion_register_loop(numpy, "multiply", loop_octonion_binary, 3, 0, octonion_mul) ||
        pyoctonion_register_loop(numpy, "divide", loop_octonion_binary, 3, 0, octonion_div) ||
        pyoctonion_register_loop(numpy, "power", loop_octonion_binary, 3, 0, octonion_rpow) ||
        pyoctonion_register_loop(numpy, "positive", loop_octonion_unary, 2, 0, octonion_pos) ||
        pyoctonion_register_loop(numpy, "negative", loop_octonion_unary, 2, 0, octonion_neg) ||
        pyoctonion_register_loop(numpy, "conjugate", loop_octonion_unary, 2, 0, octonion_conjugate) ||
        pyoctonion_register_loop(numpy, "reciprocal", loop_octonion_unary, 2, 0, octonion_inv) ||
        pyoctonion_register_loop(numpy, "square", loop_octonion_unary, 2, 0, octonion_square) ||
        pyoctonion_register_loop(numpy, "exp", loop_octonion_unary, 2, 0, octonion_exp) ||
        pyoctonion_register_loop(numpy, "expm1", loop_octonion_unary, 2, 0, octonion_expm1) ||
        pyoctonion_register_loop(numpy, "log", loop_octonion_unary, 2, 0, octonion_log) ||
        pyoctonion_register_loop(numpy, "log1p", loop_octonion_unary, 2, 0, octonion_log1p) ||
        pyoctonion_register_loop(numpy, "sqrt", loop_octonion_unary, 2, 0, octonion_sqrt) ||
        pyoctonion_register_loop(numpy, "sin", loop_octonion_unary, 2, 0, octonion_sin) ||
        pyoctonion_register_loop(numpy, "cos", loop_octonion_unary, 2, 0, octonion_cos) ||
        pyoctonion_register_loop(numpy, "tan", loop_octonion_unary, 2, 0, octonion_tan) ||
        pyoctonion_register_loop(numpy, "arcsin", loop_octonion_unary, 2, 0, octonion_asin) ||
        pyoctonion_register_loop(numpy, "arccos", loop_octonion_unary, 2, 0, octonion_acos) ||
        pyoctonion_register_loop(numpy, "arctan", loop_octonion_unary, 2, 0, octonion_atan) ||
        pyoctonion_register_loop(numpy, "sinh", loop_octonion_unary, 2, 0, octonion_sinh) ||
        pyoctonion_register_loop(numpy, "cosh", loop_octonion_unary, 2, 0, octonion_cosh) ||
        pyoctonion_register_loop(numpy, "tanh", loop_octonion_unary, 2, 0, octonion_tanh) ||
        pyoctonion_register_loop(numpy, "arcsinh", loop_octonion_unary, 2, 0, octonion_asinh) ||
        pyoctonion_register_loop(numpy, "arccosh", loop_octonion_unary, 2, 0, octonion_acosh) ||
        pyoctonion_register_loop(numpy, "arctanh", loop_octonion_unary, 2, 0, octonion_atanh) ||
        pyoctonion_register_loop(numpy, "equal", loop_octonion_predicate, 3, 1, octonion_equal) ||
        pyoctonion_register_loop(numpy, "not_equal", loop_octonion_predicate, 3, 1, octonion_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_octonion", "octonion number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__octonion(void)
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
    PyOctonion_Type = PyType_FromModuleAndSpec(
        m, &pyoctonion_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyOctonion_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyoctonion_arrfuncs);
    pyoctonion_arrfuncs.getitem = npyoctonion_getitem;
    pyoctonion_arrfuncs.setitem = npyoctonion_setitem;
    pyoctonion_arrfuncs.copyswapn = pyoctonion_copyswapn;
    pyoctonion_arrfuncs.copyswap = pyoctonion_copyswap;
    pyoctonion_arrfuncs.nonzero = pyoctonion_nonzero;
    pyoctonion_descr_proto.typeobj = (PyTypeObject *)PyOctonion_Type;
    Py_SET_TYPE((PyObject *)&pyoctonion_descr_proto, &PyArrayDescr_Type);
    npy_octonion = PyArray_RegisterDataType(&pyoctonion_descr_proto);
    if (npy_octonion < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_octonion);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyOctonion_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyOctonion_Type, 8) < 0) {
        goto fail;
    }
    if (pyoctonion_register_casts() < 0 ||
            pyoctonion_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "octonion", PyOctonion_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
