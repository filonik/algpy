/*
 * dual3 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in dual3.h.
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

#include "dual3.h"

typedef struct {
    PyObject_HEAD
    dual3 v;
} PyDual3;

static PyObject *PyDual3_Type = NULL;
static int npy_dual3 = -1;

static PyObject *
PyDual3_FromValue(PyTypeObject *type, dual3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyDual3 *p = (PyDual3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_DUAL3(v) \
    PyDual3_FromValue((PyTypeObject *)PyDual3_Type, (v))

static int
dual3_from_object(PyObject *o, dual3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyDual3_Type)) {
        *out = ((PyDual3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = dual3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pydual3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyDual3_FromValue(type, dual3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]));
}

static PyObject *
pydual3_repr(PyObject *self)
{
    dual3 v = ((PyDual3 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("dual3(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pydual3_hash(PyObject *self)
{
    dual3 v = ((PyDual3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_DUAL3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (dual3_is_real(v)) {
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
pydual3_richcompare(PyObject *a, PyObject *b, int op)
{
    dual3 x, y;
    int r;
    if ((r = dual3_from_object(a, &x)) <= 0 ||
            (r = dual3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(dual3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!dual3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define DUAL3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pydual3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        dual3 x, y;                                                       \
        int r;                                                               \
        if ((r = dual3_from_object(a, &x)) <= 0 ||                        \
                (r = dual3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_DUAL3(expr(x, y));                                     \
    }

DUAL3_BINARY(add, dual3_add)
DUAL3_BINARY(sub, dual3_sub)
DUAL3_BINARY(mul, dual3_mul)
DUAL3_BINARY(div, dual3_div)

static PyObject *
pydual3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    dual3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "dual3 does not support modular power");
        return NULL;
    }
    if ((r = dual3_from_object(a, &x)) <= 0 ||
            (r = dual3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_DUAL3(dual3_rpow(x, y));
}

#define DUAL3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pydual3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_DUAL3(expr(((PyDual3 *)self)->v));                  \
    }

DUAL3_UNARY(pos, dual3_pos)
DUAL3_UNARY(neg, dual3_neg)

static PyObject *
pydual3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue((PyTypeObject *)cls, dual3_zero());
}

static PyObject *
pydual3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue((PyTypeObject *)cls, dual3_one());
}

static int
pydual3_bool(PyObject *self)
{
    dual3 v = ((PyDual3 *)self)->v;
    return vec8_any_nonzero(v);
}

#define DUAL3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pydual3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_DUAL3(op(((PyDual3 *)self)->v));                   \
    }

DUAL3_METHOD(conjugate, dual3_conjugate)
DUAL3_METHOD(exp, dual3_exp)
DUAL3_METHOD(log, dual3_log)
DUAL3_METHOD(sqrt, dual3_sqrt)

static PyObject *
pydual3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual3_norm(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_DUAL3(dual3_inv(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual3_lmat(((PyDual3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual3_rmat(((PyDual3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {9};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual3_minpoly(((PyDual3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[8];
    npy_intp dims[1] = {8};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    dual3_eivals(((PyDual3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pydual3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual3);
    dual3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual3 *)PyArray_DATA(arr);
    dual3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_dual3);
    dual3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (dual3 *)PyArray_DATA(arr);
    dual3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pydual3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual3 v = ((PyDual3 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pydual3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_square(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_expm1(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_log1p(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_sin(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_cos(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_tan(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_asin(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_acos(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_atan(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_sinh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_cosh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_tanh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_asinh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_acosh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_atanh(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_re(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual3_FromValue(Py_TYPE(self), dual3_im(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_real(PyObject *self, void *closure)
{
    return pydual3_re(self, NULL);
}

static PyObject *
pydual3_imag(PyObject *self, void *closure)
{
    return pydual3_im(self, NULL);
}

static PyGetSetDef pydual3_getset[] = {
    {"real", pydual3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pydual3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pydual3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(dual3_is_real(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual3_to_re(((PyDual3 *)self)->v));
}

static PyObject *
pydual3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual3_FromValue((PyTypeObject *)cls, dual3_from_re(a));
}

static PyObject *
pydual3_lpow(PyObject *self, PyObject *other)
{
    dual3 b;
    int r;
    r = dual3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual3_FromValue(Py_TYPE(self),
        dual3_lpow(((PyDual3 *)self)->v, b));
}

static PyObject *
pydual3_rpow(PyObject *self, PyObject *other)
{
    dual3 b;
    int r;
    r = dual3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual3_FromValue(Py_TYPE(self),
        dual3_rpow(((PyDual3 *)self)->v, b));
}

static PyObject *
pydual3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    dual3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = dual3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual3_FromValue(Py_TYPE(self),
        dual3_mix(((PyDual3 *)self)->v, b, t));
}

static PyObject *
pydual3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_dual3, 8);
}

static PyObject *
pydual3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_dual3, 8);
}

static PyObject *
pydual3_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_DUAL3(dual3_ipow(((PyDual3 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pydual3_methods[] = {
    {"to_float", (PyCFunction)pydual3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pydual3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pydual3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pydual3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pydual3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pydual3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pydual3_square, METH_NOARGS, "square."},
    {"expm1", pydual3_expm1, METH_NOARGS, "expm1."},
    {"log1p", pydual3_log1p, METH_NOARGS, "log1p."},
    {"sin", pydual3_sin, METH_NOARGS, "sin."},
    {"cos", pydual3_cos, METH_NOARGS, "cos."},
    {"tan", pydual3_tan, METH_NOARGS, "tan."},
    {"asin", pydual3_asin, METH_NOARGS, "asin."},
    {"acos", pydual3_acos, METH_NOARGS, "acos."},
    {"atan", pydual3_atan, METH_NOARGS, "atan."},
    {"sinh", pydual3_sinh, METH_NOARGS, "sinh."},
    {"cosh", pydual3_cosh, METH_NOARGS, "cosh."},
    {"tanh", pydual3_tanh, METH_NOARGS, "tanh."},
    {"asinh", pydual3_asinh, METH_NOARGS, "asinh."},
    {"acosh", pydual3_acosh, METH_NOARGS, "acosh."},
    {"atanh", pydual3_atanh, METH_NOARGS, "atanh."},
    {"re", pydual3_re, METH_NOARGS, "re."},
    {"im", pydual3_im, METH_NOARGS, "im."},
    {"arcsin", pydual3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pydual3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pydual3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pydual3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pydual3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pydual3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pydual3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pydual3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pydual3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pydual3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pydual3_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pydual3_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pydual3_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pydual3_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pydual3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pydual3_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pydual3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pydual3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pydual3_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pydual3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pydual3_exp, METH_NOARGS, "exp."},
    {"log", pydual3_log, METH_NOARGS, "log."},
    {"sqrt", pydual3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pydual3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pydual3_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pydual3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyDual3 *)self)->v.data[i]);
}

static PyObject *
pydual3_subscript(PyObject *self, PyObject *key)
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
    return pydual3_item(self, i);
}

static int
pydual3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pydual3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pydual3_slots[] = {
    {Py_tp_dealloc, pydual3_dealloc},
    {Py_tp_traverse, pydual3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pydual3_repr},
    {Py_tp_str, pydual3_repr},
    {Py_tp_hash, pydual3_hash},
    {Py_tp_doc, (void *)"dual3 number with float64 components"},
    {Py_tp_richcompare, pydual3_richcompare},
    {Py_tp_methods, pydual3_methods},
    {Py_tp_getset, pydual3_getset},
    {Py_sq_length, pydual3_length},
    {Py_sq_item, pydual3_item},
    {Py_mp_subscript, pydual3_subscript},
    {Py_tp_new, pydual3_new},
    {Py_nb_add, pydual3_add},
    {Py_nb_subtract, pydual3_sub},
    {Py_nb_multiply, pydual3_mul},
    {Py_nb_true_divide, pydual3_div},
    {Py_nb_power, (void *)pydual3_pow},
    {Py_nb_positive, pydual3_pos},
    {Py_nb_negative, pydual3_neg},
    {Py_nb_bool, pydual3_bool},
    {0, NULL}
};

static PyType_Spec pydual3_spec = {
    .name = "algpy.algebras.multicomplex.dual3",
    .basicsize = sizeof(PyDual3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pydual3_slots,
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
npydual3_getitem(void *data, void *arr)
{
    dual3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_DUAL3(v);
}

static int
npydual3_setitem(PyObject *item, void *data, void *arr)
{
    dual3 v;
    int r = dual3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to dual3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "dual3 needs exactly 8 components");
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
pydual3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(dual3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pydual3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pydual3_copyswapn(dst, sizeof(dual3), src, sizeof(dual3), 1,
                         swap, arr);
}

static npy_bool
pydual3_nonzero(void *data, void *arr)
{
    dual3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pydual3_arrfuncs;
typedef struct { char c; dual3 v; } pydual3_align_test;

static PyArray_DescrProto pydual3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(dual3), offsetof(pydual3_align_test, v), 0, 0, 0,
    &pydual3_arrfuncs,
};

#define DUAL3_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pydual3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        dual3 *dst = (dual3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = dual3_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define DUAL3_REAL_TYPES(X) \
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

DUAL3_REAL_TYPES(DUAL3_CAST_FROM)

static int
pydual3_register_casts(void)
{
#define REGISTER_DUAL3_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_dual3,                    \
                                     pydual3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_dual3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    DUAL3_REAL_TYPES(REGISTER_DUAL3_CAST)
#undef REGISTER_DUAL3_CAST
    return 0;
}

typedef dual3 (*dual3_unary_op)(dual3);
typedef dual3 (*dual3_binary_op)(dual3, dual3);
typedef int (*dual3_predicate)(dual3, dual3);

static void
loop_dual3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    dual3_unary_op op = (dual3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_dual3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    dual3_binary_op op = (dual3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual3 a, b, r;
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
loop_dual3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    dual3_predicate op = (dual3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pydual3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_dual3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_dual3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pydual3_register_ufuncs(PyObject *numpy)
{
    return
        pydual3_register_loop(numpy, "add", loop_dual3_binary, 3, 0, dual3_add) ||
        pydual3_register_loop(numpy, "subtract", loop_dual3_binary, 3, 0, dual3_sub) ||
        pydual3_register_loop(numpy, "multiply", loop_dual3_binary, 3, 0, dual3_mul) ||
        pydual3_register_loop(numpy, "divide", loop_dual3_binary, 3, 0, dual3_div) ||
        pydual3_register_loop(numpy, "power", loop_dual3_binary, 3, 0, dual3_rpow) ||
        pydual3_register_loop(numpy, "positive", loop_dual3_unary, 2, 0, dual3_pos) ||
        pydual3_register_loop(numpy, "negative", loop_dual3_unary, 2, 0, dual3_neg) ||
        pydual3_register_loop(numpy, "conjugate", loop_dual3_unary, 2, 0, dual3_conjugate) ||
        pydual3_register_loop(numpy, "reciprocal", loop_dual3_unary, 2, 0, dual3_inv) ||
        pydual3_register_loop(numpy, "square", loop_dual3_unary, 2, 0, dual3_square) ||
        pydual3_register_loop(numpy, "exp", loop_dual3_unary, 2, 0, dual3_exp) ||
        pydual3_register_loop(numpy, "expm1", loop_dual3_unary, 2, 0, dual3_expm1) ||
        pydual3_register_loop(numpy, "log", loop_dual3_unary, 2, 0, dual3_log) ||
        pydual3_register_loop(numpy, "log1p", loop_dual3_unary, 2, 0, dual3_log1p) ||
        pydual3_register_loop(numpy, "sqrt", loop_dual3_unary, 2, 0, dual3_sqrt) ||
        pydual3_register_loop(numpy, "sin", loop_dual3_unary, 2, 0, dual3_sin) ||
        pydual3_register_loop(numpy, "cos", loop_dual3_unary, 2, 0, dual3_cos) ||
        pydual3_register_loop(numpy, "tan", loop_dual3_unary, 2, 0, dual3_tan) ||
        pydual3_register_loop(numpy, "arcsin", loop_dual3_unary, 2, 0, dual3_asin) ||
        pydual3_register_loop(numpy, "arccos", loop_dual3_unary, 2, 0, dual3_acos) ||
        pydual3_register_loop(numpy, "arctan", loop_dual3_unary, 2, 0, dual3_atan) ||
        pydual3_register_loop(numpy, "sinh", loop_dual3_unary, 2, 0, dual3_sinh) ||
        pydual3_register_loop(numpy, "cosh", loop_dual3_unary, 2, 0, dual3_cosh) ||
        pydual3_register_loop(numpy, "tanh", loop_dual3_unary, 2, 0, dual3_tanh) ||
        pydual3_register_loop(numpy, "arcsinh", loop_dual3_unary, 2, 0, dual3_asinh) ||
        pydual3_register_loop(numpy, "arccosh", loop_dual3_unary, 2, 0, dual3_acosh) ||
        pydual3_register_loop(numpy, "arctanh", loop_dual3_unary, 2, 0, dual3_atanh) ||
        pydual3_register_loop(numpy, "equal", loop_dual3_predicate, 3, 1, dual3_equal) ||
        pydual3_register_loop(numpy, "not_equal", loop_dual3_predicate, 3, 1, dual3_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_dual3", "dual3 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__dual3(void)
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
    PyDual3_Type = PyType_FromModuleAndSpec(
        m, &pydual3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyDual3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pydual3_arrfuncs);
    pydual3_arrfuncs.getitem = npydual3_getitem;
    pydual3_arrfuncs.setitem = npydual3_setitem;
    pydual3_arrfuncs.copyswapn = pydual3_copyswapn;
    pydual3_arrfuncs.copyswap = pydual3_copyswap;
    pydual3_arrfuncs.nonzero = pydual3_nonzero;
    pydual3_descr_proto.typeobj = (PyTypeObject *)PyDual3_Type;
    Py_SET_TYPE((PyObject *)&pydual3_descr_proto, &PyArrayDescr_Type);
    npy_dual3 = PyArray_RegisterDataType(&pydual3_descr_proto);
    if (npy_dual3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_dual3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyDual3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyDual3_Type, 8) < 0) {
        goto fail;
    }
    if (pydual3_register_casts() < 0 ||
            pydual3_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "dual3", PyDual3_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
