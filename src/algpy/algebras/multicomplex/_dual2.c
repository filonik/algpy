/*
 * Bi-dual numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in dual2.h.
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

#include "dual2.h"

typedef struct {
    PyObject_HEAD
    dual2 v;
} PyDual2;

static PyObject *PyDual2_Type = NULL;
static int npy_dual2 = -1;

static PyObject *
PyDual2_FromValue(PyTypeObject *type, dual2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyDual2 *p = (PyDual2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_DUAL2(v) \
    PyDual2_FromValue((PyTypeObject *)PyDual2_Type, (v))

static int
dual2_from_object(PyObject *o, dual2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyDual2_Type)) {
        *out = ((PyDual2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = dual2_make(d, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pydual2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"w", "x", "y", "z", NULL};
    double w = 0, x = 0, y = 0, z = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                     &w, &x, &y, &z)) {
        return NULL;
    }
    return PyDual2_FromValue(type, dual2_make(w, x, y, z));
}

static PyObject *
pydual2_repr(PyObject *self)
{
    dual2 v = ((PyDual2 *)self)->v;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("dual2(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pydual2_hash(PyObject *self)
{
    dual2 v = ((PyDual2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_DUAL2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (v.data[1] == 0 && v.data[2] == 0 && v.data[3] == 0) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddd)", v.data[0], v.data[1], v.data[2], v.data[3]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pydual2_richcompare(PyObject *a, PyObject *b, int op)
{
    dual2 x, y;
    int r;
    if ((r = dual2_from_object(a, &x)) <= 0 ||
            (r = dual2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(dual2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!dual2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define DUAL2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pydual2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        dual2 x, y;                                                       \
        int r;                                                               \
        if ((r = dual2_from_object(a, &x)) <= 0 ||                        \
                (r = dual2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_DUAL2(expr(x, y));                                     \
    }

DUAL2_BINARY(add, dual2_add)
DUAL2_BINARY(sub, dual2_sub)
DUAL2_BINARY(mul, dual2_mul)
DUAL2_BINARY(div, dual2_div)

static PyObject *
pydual2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    dual2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "dual2 does not support modular power");
        return NULL;
    }
    if ((r = dual2_from_object(a, &x)) <= 0 ||
            (r = dual2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_DUAL2(dual2_rpow(x, y));
}

#define DUAL2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pydual2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_DUAL2(expr(((PyDual2 *)self)->v));                  \
    }

DUAL2_UNARY(pos, dual2_pos)
DUAL2_UNARY(neg, dual2_neg)

static PyObject *
pydual2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue((PyTypeObject *)cls, dual2_zero());
}

static PyObject *
pydual2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue((PyTypeObject *)cls, dual2_one());
}

static int
pydual2_bool(PyObject *self)
{
    dual2 v = ((PyDual2 *)self)->v;
    return vec4_any_nonzero(v);
}

#define DUAL2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pydual2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_DUAL2(op(((PyDual2 *)self)->v));                   \
    }

DUAL2_METHOD(conjugate, dual2_conjugate)
DUAL2_METHOD(exp, dual2_exp)
DUAL2_METHOD(log, dual2_log)
DUAL2_METHOD(sqrt, dual2_sqrt)

static PyObject *
pydual2_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 n = dual2_norm2(((PyDual2 *)self)->v);
    return Py_BuildValue("(dd)", n.data[0], n.data[1]);
}

static PyObject *
pydual2_norm4(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual2_norm4(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual2_norm(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_DUAL2(dual2_inv(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual2_lmat(((PyDual2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual2_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual2_rmat(((PyDual2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual2_minpoly(((PyDual2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    dual2_eivals(((PyDual2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pydual2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual2);
    dual2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual2 *)PyArray_DATA(arr);
    dual2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_dual2);
    dual2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (dual2 *)PyArray_DATA(arr);
    dual2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pydual2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual2 v = ((PyDual2 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyObject *
pydual2_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_square(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_expm1(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_log1p(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_sin(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_cos(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_tan(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_asin(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_acos(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_atan(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_sinh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_cosh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_tanh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_asinh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_acosh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_atanh(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_re(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual2_FromValue(Py_TYPE(self), dual2_im(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_real(PyObject *self, void *closure)
{
    return pydual2_re(self, NULL);
}

static PyObject *
pydual2_imag(PyObject *self, void *closure)
{
    return pydual2_im(self, NULL);
}

static PyGetSetDef pydual2_getset[] = {
    {"real", pydual2_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pydual2_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pydual2_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(dual2_is_real(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual2_to_re(((PyDual2 *)self)->v));
}

static PyObject *
pydual2_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual2_FromValue((PyTypeObject *)cls, dual2_from_re(a));
}

static PyObject *
pydual2_lpow(PyObject *self, PyObject *other)
{
    dual2 b;
    int r;
    r = dual2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual2_FromValue(Py_TYPE(self),
        dual2_lpow(((PyDual2 *)self)->v, b));
}

static PyObject *
pydual2_rpow(PyObject *self, PyObject *other)
{
    dual2 b;
    int r;
    r = dual2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual2_FromValue(Py_TYPE(self),
        dual2_rpow(((PyDual2 *)self)->v, b));
}

static PyObject *
pydual2_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    dual2 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = dual2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual2_FromValue(Py_TYPE(self),
        dual2_mix(((PyDual2 *)self)->v, b, t));
}

static PyObject *
pydual2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_dual2, 4);
}

static PyObject *
pydual2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_dual2, 4);
}

static PyObject *
pydual2_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_DUAL2(dual2_ipow(((PyDual2 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pydual2_methods[] = {
    {"to_float", (PyCFunction)pydual2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pydual2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pydual2_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pydual2_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pydual2_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pydual2_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pydual2_square, METH_NOARGS, "square."},
    {"expm1", pydual2_expm1, METH_NOARGS, "expm1."},
    {"log1p", pydual2_log1p, METH_NOARGS, "log1p."},
    {"sin", pydual2_sin, METH_NOARGS, "sin."},
    {"cos", pydual2_cos, METH_NOARGS, "cos."},
    {"tan", pydual2_tan, METH_NOARGS, "tan."},
    {"asin", pydual2_asin, METH_NOARGS, "asin."},
    {"acos", pydual2_acos, METH_NOARGS, "acos."},
    {"atan", pydual2_atan, METH_NOARGS, "atan."},
    {"sinh", pydual2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pydual2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pydual2_tanh, METH_NOARGS, "tanh."},
    {"asinh", pydual2_asinh, METH_NOARGS, "asinh."},
    {"acosh", pydual2_acosh, METH_NOARGS, "acosh."},
    {"atanh", pydual2_atanh, METH_NOARGS, "atanh."},
    {"re", pydual2_re, METH_NOARGS, "re."},
    {"im", pydual2_im, METH_NOARGS, "im."},
    {"arcsin", pydual2_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pydual2_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pydual2_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pydual2_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pydual2_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pydual2_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pydual2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pydual2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pydual2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pydual2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pydual2_norm2, METH_NOARGS, "Squared norm (a dual pair)."},
    {"norm4", pydual2_norm4, METH_NOARGS, "Fourth power of the norm."},
    {"norm", pydual2_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pydual2_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pydual2_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pydual2_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pydual2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pydual2_conjugate, METH_NOARGS, "Conjugate of the e1 part."},
    {"lmat", pydual2_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pydual2_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pydual2_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pydual2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pydual2_exp, METH_NOARGS, "exp."},
    {"log", pydual2_log, METH_NOARGS, "log."},
    {"sqrt", pydual2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pydual2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pydual2_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pydual2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyDual2 *)self)->v.data[i]);
}

static PyObject *
pydual2_subscript(PyObject *self, PyObject *key)
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
        i += 4;
    }
    return pydual2_item(self, i);
}

static int
pydual2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pydual2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pydual2_slots[] = {
    {Py_tp_dealloc, pydual2_dealloc},
    {Py_tp_traverse, pydual2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pydual2_repr},
    {Py_tp_str, pydual2_repr},
    {Py_tp_hash, pydual2_hash},
    {Py_tp_doc, (void *)"Bi-dual number with float64 components"},
    {Py_tp_richcompare, pydual2_richcompare},
    {Py_tp_methods, pydual2_methods},
    {Py_tp_getset, pydual2_getset},
    {Py_sq_length, pydual2_length},
    {Py_sq_item, pydual2_item},
    {Py_mp_subscript, pydual2_subscript},
    {Py_tp_new, pydual2_new},
    {Py_nb_add, pydual2_add},
    {Py_nb_subtract, pydual2_sub},
    {Py_nb_multiply, pydual2_mul},
    {Py_nb_true_divide, pydual2_div},
    {Py_nb_power, (void *)pydual2_pow},
    {Py_nb_positive, pydual2_pos},
    {Py_nb_negative, pydual2_neg},
    {Py_nb_bool, pydual2_bool},
    {0, NULL}
};

static PyType_Spec pydual2_spec = {
    .name = "algpy.algebras.multicomplex.dual2",
    .basicsize = sizeof(PyDual2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pydual2_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npydual2_getitem(void *data, void *arr)
{
    dual2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_DUAL2(v);
}

static int
npydual2_setitem(PyObject *item, void *data, void *arr)
{
    dual2 v;
    int r = dual2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to dual2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "dual2 needs exactly 4 components");
            return -1;
        }
        for (int i = 0; i < 4; i++) {
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
pydual2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(dual2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pydual2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pydual2_copyswapn(dst, sizeof(dual2), src, sizeof(dual2), 1,
                         swap, arr);
}

static npy_bool
pydual2_nonzero(void *data, void *arr)
{
    dual2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pydual2_arrfuncs;
typedef struct { char c; dual2 v; } pydual2_align_test;

static PyArray_DescrProto pydual2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(dual2), offsetof(pydual2_align_test, v), 0, 0, 0,
    &pydual2_arrfuncs,
};

#define DUAL2_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pydual2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        dual2 *dst = (dual2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = dual2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define DUAL2_REAL_TYPES(X) \
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

DUAL2_REAL_TYPES(DUAL2_CAST_FROM)

static int
pydual2_register_casts(void)
{
#define REGISTER_DUAL2_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_dual2,                    \
                                     pydual2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_dual2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    DUAL2_REAL_TYPES(REGISTER_DUAL2_CAST)
#undef REGISTER_DUAL2_CAST
    return 0;
}

typedef dual2 (*dual2_unary_op)(dual2);
typedef dual2 (*dual2_binary_op)(dual2, dual2);
typedef int (*dual2_predicate)(dual2, dual2);

static void
loop_dual2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    dual2_unary_op op = (dual2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_dual2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    dual2_binary_op op = (dual2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual2 a, b, r;
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
loop_dual2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    dual2_predicate op = (dual2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pydual2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_dual2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_dual2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pydual2_register_ufuncs(PyObject *numpy)
{
    return
        pydual2_register_loop(numpy, "add", loop_dual2_binary, 3, 0, dual2_add) ||
        pydual2_register_loop(numpy, "subtract", loop_dual2_binary, 3, 0, dual2_sub) ||
        pydual2_register_loop(numpy, "multiply", loop_dual2_binary, 3, 0, dual2_mul) ||
        pydual2_register_loop(numpy, "divide", loop_dual2_binary, 3, 0, dual2_div) ||
        pydual2_register_loop(numpy, "power", loop_dual2_binary, 3, 0, dual2_rpow) ||
        pydual2_register_loop(numpy, "positive", loop_dual2_unary, 2, 0, dual2_pos) ||
        pydual2_register_loop(numpy, "negative", loop_dual2_unary, 2, 0, dual2_neg) ||
        pydual2_register_loop(numpy, "conjugate", loop_dual2_unary, 2, 0, dual2_conjugate) ||
        pydual2_register_loop(numpy, "reciprocal", loop_dual2_unary, 2, 0, dual2_inv) ||
        pydual2_register_loop(numpy, "square", loop_dual2_unary, 2, 0, dual2_square) ||
        pydual2_register_loop(numpy, "exp", loop_dual2_unary, 2, 0, dual2_exp) ||
        pydual2_register_loop(numpy, "expm1", loop_dual2_unary, 2, 0, dual2_expm1) ||
        pydual2_register_loop(numpy, "log", loop_dual2_unary, 2, 0, dual2_log) ||
        pydual2_register_loop(numpy, "log1p", loop_dual2_unary, 2, 0, dual2_log1p) ||
        pydual2_register_loop(numpy, "sqrt", loop_dual2_unary, 2, 0, dual2_sqrt) ||
        pydual2_register_loop(numpy, "sin", loop_dual2_unary, 2, 0, dual2_sin) ||
        pydual2_register_loop(numpy, "cos", loop_dual2_unary, 2, 0, dual2_cos) ||
        pydual2_register_loop(numpy, "tan", loop_dual2_unary, 2, 0, dual2_tan) ||
        pydual2_register_loop(numpy, "arcsin", loop_dual2_unary, 2, 0, dual2_asin) ||
        pydual2_register_loop(numpy, "arccos", loop_dual2_unary, 2, 0, dual2_acos) ||
        pydual2_register_loop(numpy, "arctan", loop_dual2_unary, 2, 0, dual2_atan) ||
        pydual2_register_loop(numpy, "sinh", loop_dual2_unary, 2, 0, dual2_sinh) ||
        pydual2_register_loop(numpy, "cosh", loop_dual2_unary, 2, 0, dual2_cosh) ||
        pydual2_register_loop(numpy, "tanh", loop_dual2_unary, 2, 0, dual2_tanh) ||
        pydual2_register_loop(numpy, "arcsinh", loop_dual2_unary, 2, 0, dual2_asinh) ||
        pydual2_register_loop(numpy, "arccosh", loop_dual2_unary, 2, 0, dual2_acosh) ||
        pydual2_register_loop(numpy, "arctanh", loop_dual2_unary, 2, 0, dual2_atanh) ||
        pydual2_register_loop(numpy, "equal", loop_dual2_predicate, 3, 1, dual2_equal) ||
        pydual2_register_loop(numpy, "not_equal", loop_dual2_predicate, 3, 1, dual2_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_dual2", "Bi-dual number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__dual2(void)
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
    PyDual2_Type = PyType_FromModuleAndSpec(
        m, &pydual2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyDual2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pydual2_arrfuncs);
    pydual2_arrfuncs.getitem = npydual2_getitem;
    pydual2_arrfuncs.setitem = npydual2_setitem;
    pydual2_arrfuncs.copyswapn = pydual2_copyswapn;
    pydual2_arrfuncs.copyswap = pydual2_copyswap;
    pydual2_arrfuncs.nonzero = pydual2_nonzero;
    pydual2_descr_proto.typeobj = (PyTypeObject *)PyDual2_Type;
    Py_SET_TYPE((PyObject *)&pydual2_descr_proto, &PyArrayDescr_Type);
    npy_dual2 = PyArray_RegisterDataType(&pydual2_descr_proto);
    if (npy_dual2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_dual2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyDual2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyDual2_Type, 4) < 0) {
        goto fail;
    }
    if (pydual2_register_casts() < 0 ||
            pydual2_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "dual2", PyDual2_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
