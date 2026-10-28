/*
 * Complex numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex1.h.
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

#include "complex1.h"

typedef struct {
    PyObject_HEAD
    complex1 v;
} PyComplex1;

static PyObject *PyComplex1_Type = NULL;
static int npy_complex1 = -1;

static PyObject *
PyComplex1_FromValue(PyTypeObject *type, complex1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplex1 *p = (PyComplex1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX1(v) \
    PyComplex1_FromValue((PyTypeObject *)PyComplex1_Type, (v))

static int
complex1_from_object(PyObject *o, complex1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplex1_Type)) {
        *out = ((PyComplex1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex1_make(d, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"real", "imag", NULL};
    double re = 0, im = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist, &re, &im)) {
        return NULL;
    }
    return PyComplex1_FromValue(type, complex1_make(re, im));
}

static PyObject *
pycomplex1_repr(PyObject *self)
{
    complex1 v = ((PyComplex1 *)self)->v;
    char *re = PyOS_double_to_string(v.data[0], 'g', 6, 0, NULL);
    char *im = PyOS_double_to_string(v.data[1], 'g', 6, 0, NULL);
    PyObject *result;
    if (re == NULL || im == NULL) {
        PyMem_Free(re);
        PyMem_Free(im);
        return NULL;
    }
    result = PyUnicode_FromFormat("complex1(%s, %s)", re, im);
    PyMem_Free(re);
    PyMem_Free(im);
    return result;
}

static Py_hash_t
pycomplex1_hash(PyObject *self)
{
    complex1 v = ((PyComplex1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
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
pycomplex1_richcompare(PyObject *a, PyObject *b, int op)
{
    complex1 x, y;
    int r;
    if ((r = complex1_from_object(a, &x)) <= 0 ||
            (r = complex1_from_object(b, &y)) <= 0) {
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

#define COMPLEX1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex1 x, y;                                                       \
        int r;                                                               \
        if ((r = complex1_from_object(a, &x)) <= 0 ||                        \
                (r = complex1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX1(expr(x, y));                                     \
    }

COMPLEX1_BINARY(add, complex1_add)
COMPLEX1_BINARY(sub, complex1_sub)
COMPLEX1_BINARY(mul, complex1_mul)
COMPLEX1_BINARY(div, complex1_div)

static PyObject *
pycomplex1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex1 does not support modular power");
        return NULL;
    }
    if ((r = complex1_from_object(a, &x)) <= 0 ||
            (r = complex1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX1(complex1_pow(x, y));
}

#define COMPLEX1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX1(expr(((PyComplex1 *)self)->v));                  \
    }

COMPLEX1_UNARY(pos, complex1_pos)
COMPLEX1_UNARY(neg, complex1_neg)

static PyObject *
pycomplex1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue((PyTypeObject *)cls, complex1_zero());
}

static PyObject *
pycomplex1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue((PyTypeObject *)cls, complex1_one());
}

static int
pycomplex1_bool(PyObject *self)
{
    complex1 v = ((PyComplex1 *)self)->v;
    return v.data[0] != 0 || v.data[1] != 0;
}

#define COMPLEX1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX1(op(((PyComplex1 *)self)->v));                   \
    }

COMPLEX1_METHOD(exp, complex1_exp)
COMPLEX1_METHOD(log, complex1_log)
COMPLEX1_METHOD(sin, complex1_sin)
COMPLEX1_METHOD(cos, complex1_cos)
COMPLEX1_METHOD(tan, complex1_tan)
COMPLEX1_METHOD(sinh, complex1_sinh)
COMPLEX1_METHOD(cosh, complex1_cosh)
COMPLEX1_METHOD(tanh, complex1_tanh)
COMPLEX1_METHOD(sqrt, complex1_sqrt)
COMPLEX1_METHOD(conjugate, complex1_conjugate)

static PyObject *
pycomplex1_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex1_norm2(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 v = ((PyComplex1 *)self)->v;
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_lmat(v, (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 v = ((PyComplex1 *)self)->v;
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_rmat(v, (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 v = ((PyComplex1 *)self)->v;
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_minpoly(v, (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 v = ((PyComplex1 *)self)->v;
    npy_intp dims[1] = {2};
    complex ev[2];
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_eivals(v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    npy_csetreal(&d[0], ev[0].data[0]);
    npy_csetimag(&d[0], ev[0].data[1]);
    npy_csetreal(&d[1], ev[1].data[0]);
    npy_csetimag(&d[1], ev[1].data[1]);
    return arr;
}

static PyObject *
pycomplex1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex1);
    complex1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex1 *)PyArray_DATA(arr);
    complex1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex1);
    complex1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex1 *)PyArray_DATA(arr);
    complex1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 v = ((PyComplex1 *)self)->v;
    return Py_BuildValue("(O(dd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1]);
}

static PyObject *
pycomplex1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex1_norm(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_inv(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_square(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_expm1(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_log1p(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_asin(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_acos(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_atan(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_asinh(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_acosh(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_atanh(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_re(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex1_FromValue(Py_TYPE(self), complex1_im(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_real(PyObject *self, void *closure)
{
    return pycomplex1_re(self, NULL);
}

static PyObject *
pycomplex1_imag(PyObject *self, void *closure)
{
    return pycomplex1_im(self, NULL);
}

static PyGetSetDef pycomplex1_getset[] = {
    {"real", pycomplex1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex1_is_real(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex1_to_re(((PyComplex1 *)self)->v));
}

static PyObject *
pycomplex1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplex1_FromValue((PyTypeObject *)cls, complex1_from_re(a));
}

static PyObject *
pycomplex1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex1 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex1_FromValue(Py_TYPE(self),
        complex1_mix(((PyComplex1 *)self)->v, b, t));
}

static PyObject *
pycomplex1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex1, 2);
}

static PyObject *
pycomplex1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex1, 2);
}

static PyMethodDef pycomplex1_methods[] = {
    {"to_float", (PyCFunction)pycomplex1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"inv", pycomplex1_inv, METH_NOARGS, "inv. See the C header for TODO status."},
    {"square", pycomplex1_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pycomplex1_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pycomplex1_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"asin", pycomplex1_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pycomplex1_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pycomplex1_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"asinh", pycomplex1_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pycomplex1_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pycomplex1_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pycomplex1_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pycomplex1_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pycomplex1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycomplex1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pycomplex1_norm2, METH_NOARGS, "Squared norm."},
    {"norm", pycomplex1_norm, METH_NOARGS, "Norm as a Python float."},
    {"conjugate", pycomplex1_conjugate, METH_NOARGS, "Complex conjugate."},
    {"lmat", pycomplex1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex1_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pycomplex1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex1_exp, METH_NOARGS, "Complex exponential."},
    {"log", pycomplex1_log, METH_NOARGS, "Principal complex logarithm."},
    {"sin", pycomplex1_sin, METH_NOARGS, "Complex sine."},
    {"cos", pycomplex1_cos, METH_NOARGS, "Complex cosine."},
    {"tan", pycomplex1_tan, METH_NOARGS, "Complex tangent."},
    {"sinh", pycomplex1_sinh, METH_NOARGS, "Complex hyperbolic sine."},
    {"cosh", pycomplex1_cosh, METH_NOARGS, "Complex hyperbolic cosine."},
    {"tanh", pycomplex1_tanh, METH_NOARGS, "Complex hyperbolic tangent."},
    {"sqrt", pycomplex1_sqrt, METH_NOARGS, "Principal complex square root."},
    {"__reduce__", pycomplex1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex1_length(PyObject *Py_UNUSED(self))
{
    return 2;
}

static PyObject *
pycomplex1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplex1 *)self)->v.data[i]);
}

static PyObject *
pycomplex1_subscript(PyObject *self, PyObject *key)
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
    return pycomplex1_item(self, i);
}

static int
pycomplex1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex1_slots[] = {
    {Py_tp_dealloc, pycomplex1_dealloc},
    {Py_tp_traverse, pycomplex1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex1_repr},
    {Py_tp_str, pycomplex1_repr},
    {Py_tp_hash, pycomplex1_hash},
    {Py_tp_doc, (void *)"Complex number with float64 components"},
    {Py_tp_richcompare, pycomplex1_richcompare},
    {Py_tp_methods, pycomplex1_methods},
    {Py_tp_getset, pycomplex1_getset},
    {Py_sq_length, pycomplex1_length},
    {Py_sq_item, pycomplex1_item},
    {Py_mp_subscript, pycomplex1_subscript},
    {Py_tp_new, pycomplex1_new},
    {Py_nb_add, pycomplex1_add},
    {Py_nb_subtract, pycomplex1_sub},
    {Py_nb_multiply, pycomplex1_mul},
    {Py_nb_true_divide, pycomplex1_div},
    {Py_nb_power, (void *)pycomplex1_pow},
    {Py_nb_positive, pycomplex1_pos},
    {Py_nb_negative, pycomplex1_neg},
    {Py_nb_bool, pycomplex1_bool},
    {0, NULL}
};

static PyType_Spec pycomplex1_spec = {
    .name = "algpy.algebras.multicomplex.complex1",
    .basicsize = sizeof(PyComplex1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex1_slots,
};

static PyObject *
npycomplex1_getitem(void *data, void *arr)
{
    complex1 v;
    memcpy(&v, data, sizeof(v));
    return NEW_COMPLEX1(v);
}

static int
npycomplex1_setitem(PyObject *item, void *data, void *arr)
{
    complex1 v;
    int r = complex1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex1 needs exactly 2 components");
            return -1;
        }
        v.data[0] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 0));
        v.data[1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 1));
        Py_DECREF(seq);
        if (PyErr_Occurred()) {
            return -1;
        }
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pycomplex1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 2; j++) {
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
pycomplex1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex1_copyswapn(dst, sizeof(complex1), src, sizeof(complex1), 1,
                         swap, arr);
}

static npy_bool
pycomplex1_nonzero(void *data, void *arr)
{
    complex1 v;
    memcpy(&v, data, sizeof(v));
    return (v.data[0] != 0 || v.data[1] != 0) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex1_arrfuncs;
typedef struct { char c; complex1 v; } pycomplex1_align_test;

static PyArray_DescrProto pycomplex1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex1), offsetof(pycomplex1_align_test, v), 0, 0, 0,
    &pycomplex1_arrfuncs,
};

#define COMPLEX1_CAST_FROM(name, type)                                       \
    static void                                                              \
    pycomplex1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex1 *dst = (complex1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex1_make((double)src[i], 0);                      \
        }                                                                    \
    }

#define COMPLEX1_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

COMPLEX1_REAL_TYPES(COMPLEX1_CAST_FROM)

static int
pycomplex1_register_casts(void)
{
#define REGISTER_COMPLEX1_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex1,                    \
                                     pycomplex1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX1_REAL_TYPES(REGISTER_COMPLEX1_CAST)
#undef REGISTER_COMPLEX1_CAST
    return 0;
}

typedef complex1 (*complex1_unary_op)(complex1);
typedef complex1 (*complex1_binary_op)(complex1, complex1);
typedef int (*complex1_predicate)(complex1, complex1);

static void
loop_complex1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex1_unary_op op = (complex1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex1_binary_op op = (complex1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex1 a, b, r;
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
loop_complex1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex1_predicate op = (complex1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex1_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex1_register_loop(numpy, "add", loop_complex1_binary, 3, 0, complex1_add) ||
        pycomplex1_register_loop(numpy, "subtract", loop_complex1_binary, 3, 0, complex1_sub) ||
        pycomplex1_register_loop(numpy, "multiply", loop_complex1_binary, 3, 0, complex1_mul) ||
        pycomplex1_register_loop(numpy, "divide", loop_complex1_binary, 3, 0, complex1_div) ||
        pycomplex1_register_loop(numpy, "power", loop_complex1_binary, 3, 0, complex1_pow) ||
        pycomplex1_register_loop(numpy, "positive", loop_complex1_unary, 2, 0, complex1_pos) ||
        pycomplex1_register_loop(numpy, "negative", loop_complex1_unary, 2, 0, complex1_neg) ||
        pycomplex1_register_loop(numpy, "conjugate", loop_complex1_unary, 2, 0, complex1_conjugate) ||
        pycomplex1_register_loop(numpy, "reciprocal", loop_complex1_unary, 2, 0, complex1_inv) ||
        pycomplex1_register_loop(numpy, "square", loop_complex1_unary, 2, 0, complex1_square) ||
        pycomplex1_register_loop(numpy, "exp", loop_complex1_unary, 2, 0, complex1_exp) ||
        pycomplex1_register_loop(numpy, "expm1", loop_complex1_unary, 2, 0, complex1_expm1) ||
        pycomplex1_register_loop(numpy, "log", loop_complex1_unary, 2, 0, complex1_log) ||
        pycomplex1_register_loop(numpy, "log1p", loop_complex1_unary, 2, 0, complex1_log1p) ||
        pycomplex1_register_loop(numpy, "sqrt", loop_complex1_unary, 2, 0, complex1_sqrt) ||
        pycomplex1_register_loop(numpy, "sin", loop_complex1_unary, 2, 0, complex1_sin) ||
        pycomplex1_register_loop(numpy, "cos", loop_complex1_unary, 2, 0, complex1_cos) ||
        pycomplex1_register_loop(numpy, "tan", loop_complex1_unary, 2, 0, complex1_tan) ||
        pycomplex1_register_loop(numpy, "arcsin", loop_complex1_unary, 2, 0, complex1_asin) ||
        pycomplex1_register_loop(numpy, "arccos", loop_complex1_unary, 2, 0, complex1_acos) ||
        pycomplex1_register_loop(numpy, "arctan", loop_complex1_unary, 2, 0, complex1_atan) ||
        pycomplex1_register_loop(numpy, "sinh", loop_complex1_unary, 2, 0, complex1_sinh) ||
        pycomplex1_register_loop(numpy, "cosh", loop_complex1_unary, 2, 0, complex1_cosh) ||
        pycomplex1_register_loop(numpy, "tanh", loop_complex1_unary, 2, 0, complex1_tanh) ||
        pycomplex1_register_loop(numpy, "arcsinh", loop_complex1_unary, 2, 0, complex1_asinh) ||
        pycomplex1_register_loop(numpy, "arccosh", loop_complex1_unary, 2, 0, complex1_acosh) ||
        pycomplex1_register_loop(numpy, "arctanh", loop_complex1_unary, 2, 0, complex1_atanh) ||
        pycomplex1_register_loop(numpy, "equal", loop_complex1_predicate, 3, 1, complex1_equal) ||
        pycomplex1_register_loop(numpy, "not_equal", loop_complex1_predicate, 3, 1, complex1_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex1", "Complex1 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex1(void)
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
    PyComplex1_Type = PyType_FromModuleAndSpec(
        m, &pycomplex1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplex1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex1_arrfuncs);
    pycomplex1_arrfuncs.getitem = npycomplex1_getitem;
    pycomplex1_arrfuncs.setitem = npycomplex1_setitem;
    pycomplex1_arrfuncs.copyswapn = pycomplex1_copyswapn;
    pycomplex1_arrfuncs.copyswap = pycomplex1_copyswap;
    pycomplex1_arrfuncs.nonzero = pycomplex1_nonzero;
    pycomplex1_descr_proto.typeobj = (PyTypeObject *)PyComplex1_Type;
    Py_SET_TYPE((PyObject *)&pycomplex1_descr_proto, &PyArrayDescr_Type);
    npy_complex1 = PyArray_RegisterDataType(&pycomplex1_descr_proto);
    if (npy_complex1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplex1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplex1_Type, 2) < 0) {
        goto fail;
    }
    if (pycomplex1_register_casts() < 0 ||
            pycomplex1_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex1", PyComplex1_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
