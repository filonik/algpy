/*
 * Dual numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in dual1.h.
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

#include "dual1.h"

typedef struct {
    PyObject_HEAD
    dual1 v;
} PyDual1;

static PyObject *PyDual1_Type = NULL;
static int npy_dual1 = -1;

static PyObject *
PyDual1_FromValue(PyTypeObject *type, dual1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyDual1 *p = (PyDual1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_DUAL1(v) \
    PyDual1_FromValue((PyTypeObject *)PyDual1_Type, (v))

static int
dual1_from_object(PyObject *o, dual1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyDual1_Type)) {
        *out = ((PyDual1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = dual1_make(d, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pydual1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"real", "imag", NULL};
    double re = 0, im = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist, &re, &im)) {
        return NULL;
    }
    return PyDual1_FromValue(type, dual1_make(re, im));
}

static PyObject *
pydual1_repr(PyObject *self)
{
    dual1 v = ((PyDual1 *)self)->v;
    char *re = PyOS_double_to_string(v.data[0], 'g', 6, 0, NULL);
    char *im = PyOS_double_to_string(v.data[1], 'g', 6, 0, NULL);
    PyObject *result;
    if (re == NULL || im == NULL) {
        PyMem_Free(re);
        PyMem_Free(im);
        return NULL;
    }
    result = PyUnicode_FromFormat("dual1(%s, %s)", re, im);
    PyMem_Free(re);
    PyMem_Free(im);
    return result;
}

static Py_hash_t
pydual1_hash(PyObject *self)
{
    dual1 v = ((PyDual1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_DUAL1_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (v.data[1] == 0) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dd)", v.data[0], v.data[1]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pydual1_richcompare(PyObject *a, PyObject *b, int op)
{
    dual1 x, y;
    int r;
    if ((r = dual1_from_object(a, &x)) <= 0 ||
            (r = dual1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(x.data[0] == y.data[0] &&
                                   x.data[1] == y.data[1]);
        case Py_NE:
            return PyBool_FromLong(x.data[0] != y.data[0] ||
                                   x.data[1] != y.data[1]);
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define DUAL1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pydual1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        dual1 x, y;                                                       \
        int r;                                                               \
        if ((r = dual1_from_object(a, &x)) <= 0 ||                        \
                (r = dual1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_DUAL1(expr(x, y));                                     \
    }

DUAL1_BINARY(add, dual1_add)
DUAL1_BINARY(sub, dual1_sub)
DUAL1_BINARY(mul, dual1_mul)
DUAL1_BINARY(div, dual1_div)

static PyObject *
pydual1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    dual1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "dual1 does not support modular power");
        return NULL;
    }
    if ((r = dual1_from_object(a, &x)) <= 0 ||
            (r = dual1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_DUAL1(dual1_rpow(x, y));
}

#define DUAL1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pydual1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_DUAL1(expr(((PyDual1 *)self)->v));                  \
    }

DUAL1_UNARY(pos, dual1_pos)
DUAL1_UNARY(neg, dual1_neg)

static PyObject *
pydual1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue((PyTypeObject *)cls, dual1_zero());
}

static PyObject *
pydual1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue((PyTypeObject *)cls, dual1_one());
}

static int
pydual1_bool(PyObject *self)
{
    dual1 v = ((PyDual1 *)self)->v;
    return v.data[0] != 0 || v.data[1] != 0;
}

#define DUAL1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pydual1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_DUAL1(op(((PyDual1 *)self)->v));                   \
    }

DUAL1_METHOD(exp, dual1_exp)
DUAL1_METHOD(log, dual1_log)
DUAL1_METHOD(sin, dual1_sin)
DUAL1_METHOD(cos, dual1_cos)
DUAL1_METHOD(tan, dual1_tan)
DUAL1_METHOD(sinh, dual1_sinh)
DUAL1_METHOD(cosh, dual1_cosh)
DUAL1_METHOD(tanh, dual1_tanh)
DUAL1_METHOD(sqrt, dual1_sqrt)
DUAL1_METHOD(conjugate, dual1_conjugate)

static PyObject *
pydual1_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual1_norm2(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 v = ((PyDual1 *)self)->v;
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *d;
    if (arr == NULL) {
        return NULL;
    }
    d = (double *)PyArray_DATA((PyArrayObject *)arr);
    dual1_lmat(v, d);
    return arr;
}

static PyObject *
pydual1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 v = ((PyDual1 *)self)->v;
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double *d;
    if (arr == NULL) {
        return NULL;
    }
    d = (double *)PyArray_DATA((PyArrayObject *)arr);
    dual1_rmat(v, d);
    return arr;
}

static PyObject *
pydual1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 v = ((PyDual1 *)self)->v;
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double *d;
    if (arr == NULL) {
        return NULL;
    }
    d = (double *)PyArray_DATA((PyArrayObject *)arr);
    d[0] = dual1_norm2(v);
    d[1] = -2 * v.data[0];
    d[2] = 1;
    return arr;
}

static PyObject *
pydual1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 v = ((PyDual1 *)self)->v;
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    npy_csetreal(&d[0], v.data[0]);
    npy_csetimag(&d[0], 0);
    npy_csetreal(&d[1], v.data[0]);
    npy_csetimag(&d[1], 0);
    return arr;
}

static PyObject *
pydual1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual1);
    dual1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual1 *)PyArray_DATA(arr);
    dual1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_dual1);
    dual1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (dual1 *)PyArray_DATA(arr);
    dual1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pydual1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 v = ((PyDual1 *)self)->v;
    return Py_BuildValue("(O(dd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1]);
}

static PyObject *
pydual1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual1_norm(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_inv(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_square(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_expm1(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_log1p(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_asin(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_acos(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_atan(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_asinh(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_acosh(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_atanh(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_re(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual1_FromValue(Py_TYPE(self), dual1_im(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_real(PyObject *self, void *closure)
{
    return pydual1_re(self, NULL);
}

static PyObject *
pydual1_imag(PyObject *self, void *closure)
{
    return pydual1_im(self, NULL);
}

static PyGetSetDef pydual1_getset[] = {
    {"real", pydual1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pydual1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pydual1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(dual1_is_real(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual1_to_re(((PyDual1 *)self)->v));
}

static PyObject *
pydual1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual1_FromValue((PyTypeObject *)cls, dual1_from_re(a));
}

static PyObject *
pydual1_lpow(PyObject *self, PyObject *other)
{
    dual1 b;
    int r;
    r = dual1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual1_FromValue(Py_TYPE(self),
        dual1_lpow(((PyDual1 *)self)->v, b));
}

static PyObject *
pydual1_rpow(PyObject *self, PyObject *other)
{
    dual1 b;
    int r;
    r = dual1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual1_FromValue(Py_TYPE(self),
        dual1_rpow(((PyDual1 *)self)->v, b));
}

static PyObject *
pydual1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    dual1 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = dual1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual1_FromValue(Py_TYPE(self),
        dual1_mix(((PyDual1 *)self)->v, b, t));
}

static PyObject *
pydual1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_dual1, 2);
}

static PyObject *
pydual1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_dual1, 2);
}

static PyObject *
pydual1_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual1_FromValue(Py_TYPE(self), dual1_ipow(((PyDual1 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pydual1_methods[] = {
    {"to_float", (PyCFunction)pydual1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pydual1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pydual1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pydual1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pydual1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pydual1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pydual1_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pydual1_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pydual1_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pydual1_inv, METH_NOARGS, "inv."},
    {"square", pydual1_square, METH_NOARGS, "square."},
    {"expm1", pydual1_expm1, METH_NOARGS, "expm1."},
    {"log1p", pydual1_log1p, METH_NOARGS, "log1p."},
    {"asin", pydual1_asin, METH_NOARGS, "asin."},
    {"acos", pydual1_acos, METH_NOARGS, "acos."},
    {"atan", pydual1_atan, METH_NOARGS, "atan."},
    {"asinh", pydual1_asinh, METH_NOARGS, "asinh."},
    {"acosh", pydual1_acosh, METH_NOARGS, "acosh."},
    {"atanh", pydual1_atanh, METH_NOARGS, "atanh."},
    {"re", pydual1_re, METH_NOARGS, "re."},
    {"im", pydual1_im, METH_NOARGS, "im."},
    {"arcsin", pydual1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pydual1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pydual1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pydual1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pydual1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pydual1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pydual1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pydual1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pydual1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pydual1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pydual1_norm2, METH_NOARGS, "Squared norm."},
    {"norm", pydual1_norm, METH_NOARGS, "Norm as a Python float."},
    {"conjugate", pydual1_conjugate, METH_NOARGS, "Dual conjugate."},
    {"lmat", pydual1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pydual1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pydual1_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pydual1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pydual1_exp, METH_NOARGS, "Dual exponential."},
    {"log", pydual1_log, METH_NOARGS, "Principal dual logarithm."},
    {"sin", pydual1_sin, METH_NOARGS, "Dual sine."},
    {"cos", pydual1_cos, METH_NOARGS, "Dual cosine."},
    {"tan", pydual1_tan, METH_NOARGS, "Dual tangent."},
    {"sinh", pydual1_sinh, METH_NOARGS, "Dual hyperbolic sine."},
    {"cosh", pydual1_cosh, METH_NOARGS, "Dual hyperbolic cosine."},
    {"tanh", pydual1_tanh, METH_NOARGS, "Dual hyperbolic tangent."},
    {"sqrt", pydual1_sqrt, METH_NOARGS, "Principal dual square root."},
    {"__reduce__", pydual1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pydual1_length(PyObject *Py_UNUSED(self))
{
    return 2;
}

static PyObject *
pydual1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyDual1 *)self)->v.data[i]);
}

static PyObject *
pydual1_subscript(PyObject *self, PyObject *key)
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
        i += 2;
    }
    return pydual1_item(self, i);
}

static int
pydual1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pydual1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pydual1_slots[] = {
    {Py_tp_dealloc, pydual1_dealloc},
    {Py_tp_traverse, pydual1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pydual1_repr},
    {Py_tp_str, pydual1_repr},
    {Py_tp_hash, pydual1_hash},
    {Py_tp_doc, (void *)"Dual number with float64 components"},
    {Py_tp_richcompare, pydual1_richcompare},
    {Py_tp_methods, pydual1_methods},
    {Py_tp_getset, pydual1_getset},
    {Py_sq_length, pydual1_length},
    {Py_sq_item, pydual1_item},
    {Py_mp_subscript, pydual1_subscript},
    {Py_tp_new, pydual1_new},
    {Py_nb_add, pydual1_add},
    {Py_nb_subtract, pydual1_sub},
    {Py_nb_multiply, pydual1_mul},
    {Py_nb_true_divide, pydual1_div},
    {Py_nb_power, (void *)pydual1_pow},
    {Py_nb_positive, pydual1_pos},
    {Py_nb_negative, pydual1_neg},
    {Py_nb_bool, pydual1_bool},
    {0, NULL}
};

static PyType_Spec pydual1_spec = {
    .name = "algpy.algebras.multicomplex.dual1",
    .basicsize = sizeof(PyDual1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pydual1_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npydual1_getitem(void *data, void *arr)
{
    dual1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_DUAL1(v);
}

static int
npydual1_setitem(PyObject *item, void *data, void *arr)
{
    dual1 v;
    int r = dual1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to dual1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "dual1 needs exactly 2 components");
            return -1;
        }
        v.data[0] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 0));
        v.data[1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 1));
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
pydual1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(dual1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pydual1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pydual1_copyswapn(dst, sizeof(dual1), src, sizeof(dual1), 1,
                         swap, arr);
}

static npy_bool
pydual1_nonzero(void *data, void *arr)
{
    dual1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return (v.data[0] != 0 || v.data[1] != 0) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pydual1_arrfuncs;
typedef struct { char c; dual1 v; } pydual1_align_test;

static PyArray_DescrProto pydual1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(dual1), offsetof(pydual1_align_test, v), 0, 0, 0,
    &pydual1_arrfuncs,
};

#define DUAL1_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pydual1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        dual1 *dst = (dual1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = dual1_make(convert(src[i]), 0);                      \
        }                                                                    \
    }

#define DUAL1_REAL_TYPES(X) \
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

DUAL1_REAL_TYPES(DUAL1_CAST_FROM)

static int
pydual1_register_casts(void)
{
#define REGISTER_DUAL1_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_dual1,                    \
                                     pydual1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_dual1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    DUAL1_REAL_TYPES(REGISTER_DUAL1_CAST)
#undef REGISTER_DUAL1_CAST
    return 0;
}

typedef dual1 (*dual1_unary_op)(dual1);
typedef dual1 (*dual1_binary_op)(dual1, dual1);
typedef int (*dual1_predicate)(dual1, dual1);

static void
loop_dual1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    dual1_unary_op op = (dual1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_dual1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    dual1_binary_op op = (dual1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual1 a, b, r;
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
loop_dual1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    dual1_predicate op = (dual1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pydual1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_dual1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_dual1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pydual1_register_ufuncs(PyObject *numpy)
{
    return
        pydual1_register_loop(numpy, "add", loop_dual1_binary, 3, 0, dual1_add) ||
        pydual1_register_loop(numpy, "subtract", loop_dual1_binary, 3, 0, dual1_sub) ||
        pydual1_register_loop(numpy, "multiply", loop_dual1_binary, 3, 0, dual1_mul) ||
        pydual1_register_loop(numpy, "divide", loop_dual1_binary, 3, 0, dual1_div) ||
        pydual1_register_loop(numpy, "power", loop_dual1_binary, 3, 0, dual1_rpow) ||
        pydual1_register_loop(numpy, "positive", loop_dual1_unary, 2, 0, dual1_pos) ||
        pydual1_register_loop(numpy, "negative", loop_dual1_unary, 2, 0, dual1_neg) ||
        pydual1_register_loop(numpy, "conjugate", loop_dual1_unary, 2, 0, dual1_conjugate) ||
        pydual1_register_loop(numpy, "reciprocal", loop_dual1_unary, 2, 0, dual1_inv) ||
        pydual1_register_loop(numpy, "square", loop_dual1_unary, 2, 0, dual1_square) ||
        pydual1_register_loop(numpy, "exp", loop_dual1_unary, 2, 0, dual1_exp) ||
        pydual1_register_loop(numpy, "expm1", loop_dual1_unary, 2, 0, dual1_expm1) ||
        pydual1_register_loop(numpy, "log", loop_dual1_unary, 2, 0, dual1_log) ||
        pydual1_register_loop(numpy, "log1p", loop_dual1_unary, 2, 0, dual1_log1p) ||
        pydual1_register_loop(numpy, "sqrt", loop_dual1_unary, 2, 0, dual1_sqrt) ||
        pydual1_register_loop(numpy, "sin", loop_dual1_unary, 2, 0, dual1_sin) ||
        pydual1_register_loop(numpy, "cos", loop_dual1_unary, 2, 0, dual1_cos) ||
        pydual1_register_loop(numpy, "tan", loop_dual1_unary, 2, 0, dual1_tan) ||
        pydual1_register_loop(numpy, "arcsin", loop_dual1_unary, 2, 0, dual1_asin) ||
        pydual1_register_loop(numpy, "arccos", loop_dual1_unary, 2, 0, dual1_acos) ||
        pydual1_register_loop(numpy, "arctan", loop_dual1_unary, 2, 0, dual1_atan) ||
        pydual1_register_loop(numpy, "sinh", loop_dual1_unary, 2, 0, dual1_sinh) ||
        pydual1_register_loop(numpy, "cosh", loop_dual1_unary, 2, 0, dual1_cosh) ||
        pydual1_register_loop(numpy, "tanh", loop_dual1_unary, 2, 0, dual1_tanh) ||
        pydual1_register_loop(numpy, "arcsinh", loop_dual1_unary, 2, 0, dual1_asinh) ||
        pydual1_register_loop(numpy, "arccosh", loop_dual1_unary, 2, 0, dual1_acosh) ||
        pydual1_register_loop(numpy, "arctanh", loop_dual1_unary, 2, 0, dual1_atanh) ||
        pydual1_register_loop(numpy, "equal", loop_dual1_predicate, 3, 1, dual1_equal) ||
        pydual1_register_loop(numpy, "not_equal", loop_dual1_predicate, 3, 1, dual1_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_dual1", "Dual number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__dual1(void)
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
    PyDual1_Type = PyType_FromModuleAndSpec(
        m, &pydual1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyDual1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pydual1_arrfuncs);
    pydual1_arrfuncs.getitem = npydual1_getitem;
    pydual1_arrfuncs.setitem = npydual1_setitem;
    pydual1_arrfuncs.copyswapn = pydual1_copyswapn;
    pydual1_arrfuncs.copyswap = pydual1_copyswap;
    pydual1_arrfuncs.nonzero = pydual1_nonzero;
    pydual1_descr_proto.typeobj = (PyTypeObject *)PyDual1_Type;
    Py_SET_TYPE((PyObject *)&pydual1_descr_proto, &PyArrayDescr_Type);
    npy_dual1 = PyArray_RegisterDataType(&pydual1_descr_proto);
    if (npy_dual1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_dual1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyDual1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyDual1_Type, 2) < 0) {
        goto fail;
    }
    if (pydual1_register_casts() < 0 ||
            pydual1_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "dual1", PyDual1_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
