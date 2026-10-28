/*
 * dual4 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in dual4.h.
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

#include "dual4.h"

typedef struct {
    PyObject_HEAD
    dual4 v;
} PyDual4;

static PyObject *PyDual4_Type = NULL;
static int npy_dual4 = -1;

static PyObject *
PyDual4_FromValue(PyTypeObject *type, dual4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyDual4 *p = (PyDual4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_DUAL4(v) \
    PyDual4_FromValue((PyTypeObject *)PyDual4_Type, (v))

static int
dual4_from_object(PyObject *o, dual4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyDual4_Type)) {
        *out = ((PyDual4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = dual4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pydual4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PyDual4_FromValue(type, dual4_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]));
}

static PyObject *
pydual4_repr(PyObject *self)
{
    dual4 v = ((PyDual4 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("dual4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pydual4_hash(PyObject *self)
{
    dual4 v = ((PyDual4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (dual4_is_real(v)) {
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
pydual4_richcompare(PyObject *a, PyObject *b, int op)
{
    dual4 x, y;
    int r;
    if ((r = dual4_from_object(a, &x)) <= 0 ||
            (r = dual4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(dual4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!dual4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define DUAL4_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pydual4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        dual4 x, y;                                                       \
        int r;                                                               \
        if ((r = dual4_from_object(a, &x)) <= 0 ||                        \
                (r = dual4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_DUAL4(expr(x, y));                                     \
    }

DUAL4_BINARY(add, dual4_add)
DUAL4_BINARY(sub, dual4_sub)
DUAL4_BINARY(mul, dual4_mul)
DUAL4_BINARY(div, dual4_div)

static PyObject *
pydual4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    dual4 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "dual4 does not support modular power");
        return NULL;
    }
    if ((r = dual4_from_object(a, &x)) <= 0 ||
            (r = dual4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_DUAL4(dual4_pow(x, y));
}

#define DUAL4_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pydual4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_DUAL4(expr(((PyDual4 *)self)->v));                  \
    }

DUAL4_UNARY(pos, dual4_pos)
DUAL4_UNARY(neg, dual4_neg)

static PyObject *
pydual4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue((PyTypeObject *)cls, dual4_zero());
}

static PyObject *
pydual4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue((PyTypeObject *)cls, dual4_one());
}

static int
pydual4_bool(PyObject *self)
{
    dual4 v = ((PyDual4 *)self)->v;
    return vec16_any_nonzero(v);
}

#define DUAL4_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pydual4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_DUAL4(op(((PyDual4 *)self)->v));                   \
    }

DUAL4_METHOD(conjugate, dual4_conjugate)
DUAL4_METHOD(exp, dual4_exp)
DUAL4_METHOD(log, dual4_log)
DUAL4_METHOD(sqrt, dual4_sqrt)

static PyObject *
pydual4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual4_norm(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_DUAL4(dual4_inv(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual4_lmat(((PyDual4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual4_rmat(((PyDual4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {17};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual4_minpoly(((PyDual4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[16];
    npy_intp dims[1] = {16};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    dual4_eivals(((PyDual4 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 16; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pydual4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual4);
    dual4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual4 *)PyArray_DATA(arr);
    dual4_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual4_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_dual4);
    dual4 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (dual4 *)PyArray_DATA(arr);
    dual4_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pydual4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual4 v = ((PyDual4 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pydual4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_square(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_expm1(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_log1p(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_sin(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_cos(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_tan(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_asin(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_acos(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_atan(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_sinh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_cosh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_tanh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_asinh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_acosh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_atanh(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_re(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual4_FromValue(Py_TYPE(self), dual4_im(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_real(PyObject *self, void *closure)
{
    return pydual4_re(self, NULL);
}

static PyObject *
pydual4_imag(PyObject *self, void *closure)
{
    return pydual4_im(self, NULL);
}

static PyGetSetDef pydual4_getset[] = {
    {"real", pydual4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pydual4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pydual4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(dual4_is_real(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual4_to_re(((PyDual4 *)self)->v));
}

static PyObject *
pydual4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual4_FromValue((PyTypeObject *)cls, dual4_from_re(a));
}

static PyObject *
pydual4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    dual4 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = dual4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual4_FromValue(Py_TYPE(self),
        dual4_mix(((PyDual4 *)self)->v, b, t));
}

static PyObject *
pydual4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_dual4, 16);
}

static PyObject *
pydual4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_dual4, 16);
}

static PyMethodDef pydual4_methods[] = {
    {"to_float", (PyCFunction)pydual4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pydual4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pydual4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pydual4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pydual4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pydual4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pydual4_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pydual4_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pydual4_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sin", pydual4_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pydual4_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pydual4_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pydual4_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pydual4_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pydual4_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pydual4_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pydual4_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pydual4_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pydual4_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pydual4_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pydual4_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pydual4_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pydual4_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pydual4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pydual4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pydual4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pydual4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pydual4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pydual4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pydual4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pydual4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pydual4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pydual4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pydual4_norm, METH_NOARGS, "Norm as a Python float."},
    {"inv", pydual4_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pydual4_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pydual4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pydual4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pydual4_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pydual4_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pydual4_exp, METH_NOARGS, "exp."},
    {"log", pydual4_log, METH_NOARGS, "log."},
    {"sqrt", pydual4_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pydual4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pydual4_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pydual4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyDual4 *)self)->v.data[i]);
}

static PyObject *
pydual4_subscript(PyObject *self, PyObject *key)
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
    return pydual4_item(self, i);
}

static int
pydual4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pydual4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pydual4_slots[] = {
    {Py_tp_dealloc, pydual4_dealloc},
    {Py_tp_traverse, pydual4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pydual4_repr},
    {Py_tp_str, pydual4_repr},
    {Py_tp_hash, pydual4_hash},
    {Py_tp_doc, (void *)"dual4 number with float64 components"},
    {Py_tp_richcompare, pydual4_richcompare},
    {Py_tp_methods, pydual4_methods},
    {Py_tp_getset, pydual4_getset},
    {Py_sq_length, pydual4_length},
    {Py_sq_item, pydual4_item},
    {Py_mp_subscript, pydual4_subscript},
    {Py_tp_new, pydual4_new},
    {Py_nb_add, pydual4_add},
    {Py_nb_subtract, pydual4_sub},
    {Py_nb_multiply, pydual4_mul},
    {Py_nb_true_divide, pydual4_div},
    {Py_nb_power, (void *)pydual4_pow},
    {Py_nb_positive, pydual4_pos},
    {Py_nb_negative, pydual4_neg},
    {Py_nb_bool, pydual4_bool},
    {0, NULL}
};

static PyType_Spec pydual4_spec = {
    .name = "algpy.algebras.multicomplex.dual4",
    .basicsize = sizeof(PyDual4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pydual4_slots,
};

static PyObject *
npydual4_getitem(void *data, void *arr)
{
    dual4 v;
    memcpy(&v, data, sizeof(v));
    return NEW_DUAL4(v);
}

static int
npydual4_setitem(PyObject *item, void *data, void *arr)
{
    dual4 v;
    int r = dual4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to dual4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "dual4 needs exactly 16 components");
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
pydual4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(dual4));
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
pydual4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pydual4_copyswapn(dst, sizeof(dual4), src, sizeof(dual4), 1,
                         swap, arr);
}

static npy_bool
pydual4_nonzero(void *data, void *arr)
{
    dual4 v;
    memcpy(&v, data, sizeof(v));
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pydual4_arrfuncs;
typedef struct { char c; dual4 v; } pydual4_align_test;

static PyArray_DescrProto pydual4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(dual4), offsetof(pydual4_align_test, v), 0, 0, 0,
    &pydual4_arrfuncs,
};

#define DUAL4_CAST_FROM(name, type)                                       \
    static void                                                              \
    pydual4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        dual4 *dst = (dual4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = dual4_from_re((double)src[i]);                      \
        }                                                                    \
    }

#define DUAL4_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

DUAL4_REAL_TYPES(DUAL4_CAST_FROM)

static int
pydual4_register_casts(void)
{
#define REGISTER_DUAL4_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_dual4,                    \
                                     pydual4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_dual4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    DUAL4_REAL_TYPES(REGISTER_DUAL4_CAST)
#undef REGISTER_DUAL4_CAST
    return 0;
}

typedef dual4 (*dual4_unary_op)(dual4);
typedef dual4 (*dual4_binary_op)(dual4, dual4);
typedef int (*dual4_predicate)(dual4, dual4);

static void
loop_dual4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    dual4_unary_op op = (dual4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_dual4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    dual4_binary_op op = (dual4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual4 a, b, r;
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
loop_dual4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    dual4_predicate op = (dual4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pydual4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_dual4;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_dual4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pydual4_register_ufuncs(PyObject *numpy)
{
    return
        pydual4_register_loop(numpy, "add", loop_dual4_binary, 3, 0, dual4_add) ||
        pydual4_register_loop(numpy, "subtract", loop_dual4_binary, 3, 0, dual4_sub) ||
        pydual4_register_loop(numpy, "multiply", loop_dual4_binary, 3, 0, dual4_mul) ||
        pydual4_register_loop(numpy, "divide", loop_dual4_binary, 3, 0, dual4_div) ||
        pydual4_register_loop(numpy, "power", loop_dual4_binary, 3, 0, dual4_pow) ||
        pydual4_register_loop(numpy, "positive", loop_dual4_unary, 2, 0, dual4_pos) ||
        pydual4_register_loop(numpy, "negative", loop_dual4_unary, 2, 0, dual4_neg) ||
        pydual4_register_loop(numpy, "conjugate", loop_dual4_unary, 2, 0, dual4_conjugate) ||
        pydual4_register_loop(numpy, "reciprocal", loop_dual4_unary, 2, 0, dual4_inv) ||
        pydual4_register_loop(numpy, "square", loop_dual4_unary, 2, 0, dual4_square) ||
        pydual4_register_loop(numpy, "exp", loop_dual4_unary, 2, 0, dual4_exp) ||
        pydual4_register_loop(numpy, "expm1", loop_dual4_unary, 2, 0, dual4_expm1) ||
        pydual4_register_loop(numpy, "log", loop_dual4_unary, 2, 0, dual4_log) ||
        pydual4_register_loop(numpy, "log1p", loop_dual4_unary, 2, 0, dual4_log1p) ||
        pydual4_register_loop(numpy, "sqrt", loop_dual4_unary, 2, 0, dual4_sqrt) ||
        pydual4_register_loop(numpy, "sin", loop_dual4_unary, 2, 0, dual4_sin) ||
        pydual4_register_loop(numpy, "cos", loop_dual4_unary, 2, 0, dual4_cos) ||
        pydual4_register_loop(numpy, "tan", loop_dual4_unary, 2, 0, dual4_tan) ||
        pydual4_register_loop(numpy, "arcsin", loop_dual4_unary, 2, 0, dual4_asin) ||
        pydual4_register_loop(numpy, "arccos", loop_dual4_unary, 2, 0, dual4_acos) ||
        pydual4_register_loop(numpy, "arctan", loop_dual4_unary, 2, 0, dual4_atan) ||
        pydual4_register_loop(numpy, "sinh", loop_dual4_unary, 2, 0, dual4_sinh) ||
        pydual4_register_loop(numpy, "cosh", loop_dual4_unary, 2, 0, dual4_cosh) ||
        pydual4_register_loop(numpy, "tanh", loop_dual4_unary, 2, 0, dual4_tanh) ||
        pydual4_register_loop(numpy, "arcsinh", loop_dual4_unary, 2, 0, dual4_asinh) ||
        pydual4_register_loop(numpy, "arccosh", loop_dual4_unary, 2, 0, dual4_acosh) ||
        pydual4_register_loop(numpy, "arctanh", loop_dual4_unary, 2, 0, dual4_atanh) ||
        pydual4_register_loop(numpy, "equal", loop_dual4_predicate, 3, 1, dual4_equal) ||
        pydual4_register_loop(numpy, "not_equal", loop_dual4_predicate, 3, 1, dual4_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_dual4", "dual4 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__dual4(void)
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
    PyDual4_Type = PyType_FromModuleAndSpec(
        m, &pydual4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyDual4_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pydual4_arrfuncs);
    pydual4_arrfuncs.getitem = npydual4_getitem;
    pydual4_arrfuncs.setitem = npydual4_setitem;
    pydual4_arrfuncs.copyswapn = pydual4_copyswapn;
    pydual4_arrfuncs.copyswap = pydual4_copyswap;
    pydual4_arrfuncs.nonzero = pydual4_nonzero;
    pydual4_descr_proto.typeobj = (PyTypeObject *)PyDual4_Type;
    Py_SET_TYPE((PyObject *)&pydual4_descr_proto, &PyArrayDescr_Type);
    npy_dual4 = PyArray_RegisterDataType(&pydual4_descr_proto);
    if (npy_dual4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_dual4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyDual4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyDual4_Type, 16) < 0) {
        goto fail;
    }
    if (pydual4_register_casts() < 0 ||
            pydual4_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "dual4", PyDual4_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
