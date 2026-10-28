/*
 * real numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in real.h.
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

#include "real.h"

typedef struct {
    PyObject_HEAD
    real v;
} PyReal;

static PyObject *PyReal_Type = NULL;
static int npy_real = -1;

static PyObject *
PyReal_FromValue(PyTypeObject *type, real v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyReal *p = (PyReal *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_REAL(v) \
    PyReal_FromValue((PyTypeObject *)PyReal_Type, (v))

static int
real_from_object(PyObject *o, real *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyReal_Type)) {
        *out = ((PyReal *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = real_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyreal_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[1] = {0};
    static char *kwlist[] = {"x0", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|d", kwlist,
            &d[0])) {
        return NULL;
    }
    return PyReal_FromValue(type, real_make(d[0]));
}

static PyObject *
pyreal_repr(PyObject *self)
{
    real v = ((PyReal *)self)->v;
    char *s[1] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 1; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("real(%s)",
        s[0]);
done:
    for (i = 0; i < 1; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyreal_hash(PyObject *self)
{
    real v = ((PyReal *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (real_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(d)", v.data[0]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyreal_richcompare(PyObject *a, PyObject *b, int op)
{
    real x, y;
    int r;
    if ((r = real_from_object(a, &x)) <= 0 ||
            (r = real_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(real_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!real_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define REAL_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyreal_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        real x, y;                                                       \
        int r;                                                               \
        if ((r = real_from_object(a, &x)) <= 0 ||                        \
                (r = real_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_REAL(expr(x, y));                                     \
    }

REAL_BINARY(add, real_add)
REAL_BINARY(sub, real_sub)
REAL_BINARY(mul, real_mul)
REAL_BINARY(div, real_div)

static PyObject *
pyreal_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    real x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "real does not support modular power");
        return NULL;
    }
    if ((r = real_from_object(a, &x)) <= 0 ||
            (r = real_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_REAL(real_pow(x, y));
}

#define REAL_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyreal_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_REAL(expr(((PyReal *)self)->v));                  \
    }

REAL_UNARY(pos, real_pos)
REAL_UNARY(neg, real_neg)

static PyObject *
pyreal_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue((PyTypeObject *)cls, real_zero());
}

static PyObject *
pyreal_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue((PyTypeObject *)cls, real_one());
}

static int
pyreal_bool(PyObject *self)
{
    real v = ((PyReal *)self)->v;
    return vec1_any_nonzero(v);
}

#define REAL_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyreal_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_REAL(op(((PyReal *)self)->v));                   \
    }

REAL_METHOD(conjugate, real_conjugate)
REAL_METHOD(exp, real_exp)
REAL_METHOD(log, real_log)
REAL_METHOD(sqrt, real_sqrt)

static PyObject *
pyreal_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(real_norm(((PyReal *)self)->v));
}

static PyObject *
pyreal_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_REAL(real_inv(((PyReal *)self)->v));
}

static PyObject *
pyreal_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {1, 1};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    real_lmat(((PyReal *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyreal_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {1, 1};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    real_rmat(((PyReal *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyreal_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    real_minpoly(((PyReal *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyreal_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[1];
    npy_intp dims[1] = {1};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    real_eivals(((PyReal *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 1; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyreal_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_real);
    real *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (real *)PyArray_DATA(arr);
    real_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyreal_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {1, 1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_real);
    real *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (real *)PyArray_DATA(arr);
    real_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyreal_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    real v = ((PyReal *)self)->v;
    return Py_BuildValue("(O(d))", (PyObject *)Py_TYPE(self),
        v.data[0]);
}

static PyObject *
pyreal_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_square(((PyReal *)self)->v));
}

static PyObject *
pyreal_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_expm1(((PyReal *)self)->v));
}

static PyObject *
pyreal_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_log1p(((PyReal *)self)->v));
}

static PyObject *
pyreal_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_sin(((PyReal *)self)->v));
}

static PyObject *
pyreal_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_cos(((PyReal *)self)->v));
}

static PyObject *
pyreal_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_tan(((PyReal *)self)->v));
}

static PyObject *
pyreal_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_asin(((PyReal *)self)->v));
}

static PyObject *
pyreal_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_acos(((PyReal *)self)->v));
}

static PyObject *
pyreal_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_atan(((PyReal *)self)->v));
}

static PyObject *
pyreal_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_sinh(((PyReal *)self)->v));
}

static PyObject *
pyreal_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_cosh(((PyReal *)self)->v));
}

static PyObject *
pyreal_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_tanh(((PyReal *)self)->v));
}

static PyObject *
pyreal_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_asinh(((PyReal *)self)->v));
}

static PyObject *
pyreal_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_acosh(((PyReal *)self)->v));
}

static PyObject *
pyreal_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_atanh(((PyReal *)self)->v));
}

static PyObject *
pyreal_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_re(((PyReal *)self)->v));
}

static PyObject *
pyreal_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyReal_FromValue(Py_TYPE(self), real_im(((PyReal *)self)->v));
}

static PyObject *
pyreal_real(PyObject *self, void *closure)
{
    return pyreal_re(self, NULL);
}

static PyObject *
pyreal_imag(PyObject *self, void *closure)
{
    return pyreal_im(self, NULL);
}

static PyGetSetDef pyreal_getset[] = {
    {"real", pyreal_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyreal_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyreal_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(real_is_real(((PyReal *)self)->v));
}

static PyObject *
pyreal_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(real_to_re(((PyReal *)self)->v));
}

static PyObject *
pyreal_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyReal_FromValue((PyTypeObject *)cls, real_from_re(a));
}

static PyObject *
pyreal_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    real b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = real_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyReal_FromValue(Py_TYPE(self),
        real_mix(((PyReal *)self)->v, b, t));
}

static PyObject *
pyreal_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_real, 1);
}

static PyObject *
pyreal_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_real, 1);
}

static PyMethodDef pyreal_methods[] = {
    {"to_float", (PyCFunction)pyreal_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyreal_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyreal_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyreal_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyreal_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyreal_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyreal_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pyreal_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pyreal_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sin", pyreal_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pyreal_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pyreal_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pyreal_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pyreal_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pyreal_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pyreal_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pyreal_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pyreal_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pyreal_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pyreal_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pyreal_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pyreal_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pyreal_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pyreal_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyreal_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyreal_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyreal_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyreal_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyreal_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyreal_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyreal_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyreal_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pyreal_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyreal_norm, METH_NOARGS, "Norm as a Python float."},
    {"inv", pyreal_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyreal_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pyreal_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyreal_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyreal_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pyreal_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyreal_exp, METH_NOARGS, "exp."},
    {"log", pyreal_log, METH_NOARGS, "log."},
    {"sqrt", pyreal_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyreal_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyreal_length(PyObject *Py_UNUSED(self))
{
    return 1;
}

static PyObject *
pyreal_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 1) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyReal *)self)->v.data[i]);
}

static PyObject *
pyreal_subscript(PyObject *self, PyObject *key)
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
        i += 1;
    }
    return pyreal_item(self, i);
}

static int
pyreal_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyreal_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyreal_slots[] = {
    {Py_tp_dealloc, pyreal_dealloc},
    {Py_tp_traverse, pyreal_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyreal_repr},
    {Py_tp_str, pyreal_repr},
    {Py_tp_hash, pyreal_hash},
    {Py_tp_doc, (void *)"real number with float64 components"},
    {Py_tp_richcompare, pyreal_richcompare},
    {Py_tp_methods, pyreal_methods},
    {Py_tp_getset, pyreal_getset},
    {Py_sq_length, pyreal_length},
    {Py_sq_item, pyreal_item},
    {Py_mp_subscript, pyreal_subscript},
    {Py_tp_new, pyreal_new},
    {Py_nb_add, pyreal_add},
    {Py_nb_subtract, pyreal_sub},
    {Py_nb_multiply, pyreal_mul},
    {Py_nb_true_divide, pyreal_div},
    {Py_nb_power, (void *)pyreal_pow},
    {Py_nb_positive, pyreal_pos},
    {Py_nb_negative, pyreal_neg},
    {Py_nb_bool, pyreal_bool},
    {0, NULL}
};

static PyType_Spec pyreal_spec = {
    .name = "algpy.algebras.cayleydickson.real",
    .basicsize = sizeof(PyReal),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyreal_slots,
};

static PyObject *
npyreal_getitem(void *data, void *arr)
{
    real v;
    memcpy(&v, data, sizeof(v));
    return NEW_REAL(v);
}

static int
npyreal_setitem(PyObject *item, void *data, void *arr)
{
    real v;
    int r = real_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to real");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 1) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "real needs exactly 1 components");
            return -1;
        }
        for (int i = 0; i < 1; i++) {
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
pyreal_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(real));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 1; j++) {
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
pyreal_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyreal_copyswapn(dst, sizeof(real), src, sizeof(real), 1,
                         swap, arr);
}

static npy_bool
pyreal_nonzero(void *data, void *arr)
{
    real v;
    memcpy(&v, data, sizeof(v));
    return vec1_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyreal_arrfuncs;
typedef struct { char c; real v; } pyreal_align_test;

static PyArray_DescrProto pyreal_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(real), offsetof(pyreal_align_test, v), 0, 0, 0,
    &pyreal_arrfuncs,
};

#define REAL_CAST_FROM(name, type)                                       \
    static void                                                              \
    pyreal_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        real *dst = (real *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = real_from_re((double)src[i]);                      \
        }                                                                    \
    }

#define REAL_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

REAL_REAL_TYPES(REAL_CAST_FROM)

static int
pyreal_register_casts(void)
{
#define REGISTER_REAL_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_real,                    \
                                     pyreal_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_real, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    REAL_REAL_TYPES(REGISTER_REAL_CAST)
#undef REGISTER_REAL_CAST
    return 0;
}

typedef real (*real_unary_op)(real);
typedef real (*real_binary_op)(real, real);
typedef int (*real_predicate)(real, real);

static void
loop_real_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    real_unary_op op = (real_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        real a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_real_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    real_binary_op op = (real_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        real a, b, r;
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
loop_real_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    real_predicate op = (real_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        real a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyreal_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_real;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_real,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyreal_register_ufuncs(PyObject *numpy)
{
    return
        pyreal_register_loop(numpy, "add", loop_real_binary, 3, 0, real_add) ||
        pyreal_register_loop(numpy, "subtract", loop_real_binary, 3, 0, real_sub) ||
        pyreal_register_loop(numpy, "multiply", loop_real_binary, 3, 0, real_mul) ||
        pyreal_register_loop(numpy, "divide", loop_real_binary, 3, 0, real_div) ||
        pyreal_register_loop(numpy, "power", loop_real_binary, 3, 0, real_pow) ||
        pyreal_register_loop(numpy, "positive", loop_real_unary, 2, 0, real_pos) ||
        pyreal_register_loop(numpy, "negative", loop_real_unary, 2, 0, real_neg) ||
        pyreal_register_loop(numpy, "conjugate", loop_real_unary, 2, 0, real_conjugate) ||
        pyreal_register_loop(numpy, "reciprocal", loop_real_unary, 2, 0, real_inv) ||
        pyreal_register_loop(numpy, "square", loop_real_unary, 2, 0, real_square) ||
        pyreal_register_loop(numpy, "exp", loop_real_unary, 2, 0, real_exp) ||
        pyreal_register_loop(numpy, "expm1", loop_real_unary, 2, 0, real_expm1) ||
        pyreal_register_loop(numpy, "log", loop_real_unary, 2, 0, real_log) ||
        pyreal_register_loop(numpy, "log1p", loop_real_unary, 2, 0, real_log1p) ||
        pyreal_register_loop(numpy, "sqrt", loop_real_unary, 2, 0, real_sqrt) ||
        pyreal_register_loop(numpy, "sin", loop_real_unary, 2, 0, real_sin) ||
        pyreal_register_loop(numpy, "cos", loop_real_unary, 2, 0, real_cos) ||
        pyreal_register_loop(numpy, "tan", loop_real_unary, 2, 0, real_tan) ||
        pyreal_register_loop(numpy, "arcsin", loop_real_unary, 2, 0, real_asin) ||
        pyreal_register_loop(numpy, "arccos", loop_real_unary, 2, 0, real_acos) ||
        pyreal_register_loop(numpy, "arctan", loop_real_unary, 2, 0, real_atan) ||
        pyreal_register_loop(numpy, "sinh", loop_real_unary, 2, 0, real_sinh) ||
        pyreal_register_loop(numpy, "cosh", loop_real_unary, 2, 0, real_cosh) ||
        pyreal_register_loop(numpy, "tanh", loop_real_unary, 2, 0, real_tanh) ||
        pyreal_register_loop(numpy, "arcsinh", loop_real_unary, 2, 0, real_asinh) ||
        pyreal_register_loop(numpy, "arccosh", loop_real_unary, 2, 0, real_acosh) ||
        pyreal_register_loop(numpy, "arctanh", loop_real_unary, 2, 0, real_atanh) ||
        pyreal_register_loop(numpy, "equal", loop_real_predicate, 3, 1, real_equal) ||
        pyreal_register_loop(numpy, "not_equal", loop_real_predicate, 3, 1, real_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_real", "real number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__real(void)
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
    PyReal_Type = PyType_FromModuleAndSpec(
        m, &pyreal_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyReal_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyreal_arrfuncs);
    pyreal_arrfuncs.getitem = npyreal_getitem;
    pyreal_arrfuncs.setitem = npyreal_setitem;
    pyreal_arrfuncs.copyswapn = pyreal_copyswapn;
    pyreal_arrfuncs.copyswap = pyreal_copyswap;
    pyreal_arrfuncs.nonzero = pyreal_nonzero;
    pyreal_descr_proto.typeobj = (PyTypeObject *)PyReal_Type;
    Py_SET_TYPE((PyObject *)&pyreal_descr_proto, &PyArrayDescr_Type);
    npy_real = PyArray_RegisterDataType(&pyreal_descr_proto);
    if (npy_real < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_real);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyReal_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyReal_Type, 1) < 0) {
        goto fail;
    }
    if (pyreal_register_casts() < 0 ||
            pyreal_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "real", PyReal_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
