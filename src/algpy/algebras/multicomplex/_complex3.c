/*
 * complex3 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex3.h.
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

#include "complex3.h"

typedef struct {
    PyObject_HEAD
    complex3 v;
} PyComplex3;

static PyObject *PyComplex3_Type = NULL;
static int npy_complex3 = -1;

static PyObject *
PyComplex3_FromValue(PyTypeObject *type, complex3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplex3 *p = (PyComplex3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX3(v) \
    PyComplex3_FromValue((PyTypeObject *)PyComplex3_Type, (v))

static int
complex3_from_object(PyObject *o, complex3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplex3_Type)) {
        *out = ((PyComplex3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyComplex3_FromValue(type, complex3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]));
}

static PyObject *
pycomplex3_repr(PyObject *self)
{
    complex3 v = ((PyComplex3 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex3(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex3_hash(PyObject *self)
{
    complex3 v = ((PyComplex3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (complex3_is_real(v)) {
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
pycomplex3_richcompare(PyObject *a, PyObject *b, int op)
{
    complex3 x, y;
    int r;
    if ((r = complex3_from_object(a, &x)) <= 0 ||
            (r = complex3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(complex3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!complex3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define COMPLEX3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex3 x, y;                                                       \
        int r;                                                               \
        if ((r = complex3_from_object(a, &x)) <= 0 ||                        \
                (r = complex3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX3(expr(x, y));                                     \
    }

COMPLEX3_BINARY(add, complex3_add)
COMPLEX3_BINARY(sub, complex3_sub)
COMPLEX3_BINARY(mul, complex3_mul)
COMPLEX3_BINARY(div, complex3_div)

static PyObject *
pycomplex3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex3 does not support modular power");
        return NULL;
    }
    if ((r = complex3_from_object(a, &x)) <= 0 ||
            (r = complex3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX3(complex3_pow(x, y));
}

#define COMPLEX3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX3(expr(((PyComplex3 *)self)->v));                  \
    }

COMPLEX3_UNARY(pos, complex3_pos)
COMPLEX3_UNARY(neg, complex3_neg)

static PyObject *
pycomplex3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue((PyTypeObject *)cls, complex3_zero());
}

static PyObject *
pycomplex3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue((PyTypeObject *)cls, complex3_one());
}

static int
pycomplex3_bool(PyObject *self)
{
    complex3 v = ((PyComplex3 *)self)->v;
    return vec8_any_nonzero(v);
}

#define COMPLEX3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX3(op(((PyComplex3 *)self)->v));                   \
    }

COMPLEX3_METHOD(conjugate, complex3_conjugate)
COMPLEX3_METHOD(exp, complex3_exp)
COMPLEX3_METHOD(log, complex3_log)
COMPLEX3_METHOD(sqrt, complex3_sqrt)

static PyObject *
pycomplex3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex3_norm(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX3(complex3_inv(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex3_lmat(((PyComplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex3_rmat(((PyComplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {9};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex3_minpoly(((PyComplex3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[8];
    npy_intp dims[1] = {8};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex3_eivals(((PyComplex3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex3);
    complex3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex3 *)PyArray_DATA(arr);
    complex3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex3);
    complex3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex3 *)PyArray_DATA(arr);
    complex3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex3 v = ((PyComplex3 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pycomplex3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_square(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_expm1(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_log1p(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_sin(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_cos(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_tan(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_asin(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_acos(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_atan(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_sinh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_cosh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_tanh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_asinh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_acosh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_atanh(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_re(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex3_FromValue(Py_TYPE(self), complex3_im(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_real(PyObject *self, void *closure)
{
    return pycomplex3_re(self, NULL);
}

static PyObject *
pycomplex3_imag(PyObject *self, void *closure)
{
    return pycomplex3_im(self, NULL);
}

static PyGetSetDef pycomplex3_getset[] = {
    {"real", pycomplex3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex3_is_real(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex3_to_re(((PyComplex3 *)self)->v));
}

static PyObject *
pycomplex3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplex3_FromValue((PyTypeObject *)cls, complex3_from_re(a));
}

static PyObject *
pycomplex3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex3_FromValue(Py_TYPE(self),
        complex3_mix(((PyComplex3 *)self)->v, b, t));
}

static PyObject *
pycomplex3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex3, 8);
}

static PyObject *
pycomplex3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex3, 8);
}

static PyMethodDef pycomplex3_methods[] = {
    {"to_float", (PyCFunction)pycomplex3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycomplex3_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"expm1", pycomplex3_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log1p", pycomplex3_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sin", pycomplex3_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pycomplex3_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pycomplex3_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pycomplex3_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pycomplex3_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pycomplex3_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pycomplex3_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pycomplex3_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pycomplex3_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pycomplex3_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pycomplex3_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pycomplex3_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pycomplex3_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pycomplex3_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pycomplex3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycomplex3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycomplex3_norm, METH_NOARGS, "Norm as a Python float."},
    {"inv", pycomplex3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex3_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"lmat", pycomplex3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex3_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pycomplex3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex3_exp, METH_NOARGS, "exp."},
    {"log", pycomplex3_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex3_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pycomplex3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplex3 *)self)->v.data[i]);
}

static PyObject *
pycomplex3_subscript(PyObject *self, PyObject *key)
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
    return pycomplex3_item(self, i);
}

static int
pycomplex3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex3_slots[] = {
    {Py_tp_dealloc, pycomplex3_dealloc},
    {Py_tp_traverse, pycomplex3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex3_repr},
    {Py_tp_str, pycomplex3_repr},
    {Py_tp_hash, pycomplex3_hash},
    {Py_tp_doc, (void *)"complex3 number with float64 components"},
    {Py_tp_richcompare, pycomplex3_richcompare},
    {Py_tp_methods, pycomplex3_methods},
    {Py_tp_getset, pycomplex3_getset},
    {Py_sq_length, pycomplex3_length},
    {Py_sq_item, pycomplex3_item},
    {Py_mp_subscript, pycomplex3_subscript},
    {Py_tp_new, pycomplex3_new},
    {Py_nb_add, pycomplex3_add},
    {Py_nb_subtract, pycomplex3_sub},
    {Py_nb_multiply, pycomplex3_mul},
    {Py_nb_true_divide, pycomplex3_div},
    {Py_nb_power, (void *)pycomplex3_pow},
    {Py_nb_positive, pycomplex3_pos},
    {Py_nb_negative, pycomplex3_neg},
    {Py_nb_bool, pycomplex3_bool},
    {0, NULL}
};

static PyType_Spec pycomplex3_spec = {
    .name = "algpy.algebras.multicomplex.complex3",
    .basicsize = sizeof(PyComplex3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex3_slots,
};

static PyObject *
npycomplex3_getitem(void *data, void *arr)
{
    complex3 v;
    memcpy(&v, data, sizeof(v));
    return NEW_COMPLEX3(v);
}

static int
npycomplex3_setitem(PyObject *item, void *data, void *arr)
{
    complex3 v;
    int r = complex3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex3 needs exactly 8 components");
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
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pycomplex3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 8; j++) {
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
pycomplex3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex3_copyswapn(dst, sizeof(complex3), src, sizeof(complex3), 1,
                         swap, arr);
}

static npy_bool
pycomplex3_nonzero(void *data, void *arr)
{
    complex3 v;
    memcpy(&v, data, sizeof(v));
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex3_arrfuncs;
typedef struct { char c; complex3 v; } pycomplex3_align_test;

static PyArray_DescrProto pycomplex3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex3), offsetof(pycomplex3_align_test, v), 0, 0, 0,
    &pycomplex3_arrfuncs,
};

#define COMPLEX3_CAST_FROM(name, type)                                       \
    static void                                                              \
    pycomplex3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex3 *dst = (complex3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex3_from_re((double)src[i]);                      \
        }                                                                    \
    }

#define COMPLEX3_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

COMPLEX3_REAL_TYPES(COMPLEX3_CAST_FROM)

static int
pycomplex3_register_casts(void)
{
#define REGISTER_COMPLEX3_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex3,                    \
                                     pycomplex3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX3_REAL_TYPES(REGISTER_COMPLEX3_CAST)
#undef REGISTER_COMPLEX3_CAST
    return 0;
}

typedef complex3 (*complex3_unary_op)(complex3);
typedef complex3 (*complex3_binary_op)(complex3, complex3);
typedef int (*complex3_predicate)(complex3, complex3);

static void
loop_complex3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex3_unary_op op = (complex3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex3_binary_op op = (complex3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex3 a, b, r;
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
loop_complex3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex3_predicate op = (complex3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex3_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex3_register_loop(numpy, "add", loop_complex3_binary, 3, 0, complex3_add) ||
        pycomplex3_register_loop(numpy, "subtract", loop_complex3_binary, 3, 0, complex3_sub) ||
        pycomplex3_register_loop(numpy, "multiply", loop_complex3_binary, 3, 0, complex3_mul) ||
        pycomplex3_register_loop(numpy, "divide", loop_complex3_binary, 3, 0, complex3_div) ||
        pycomplex3_register_loop(numpy, "power", loop_complex3_binary, 3, 0, complex3_pow) ||
        pycomplex3_register_loop(numpy, "positive", loop_complex3_unary, 2, 0, complex3_pos) ||
        pycomplex3_register_loop(numpy, "negative", loop_complex3_unary, 2, 0, complex3_neg) ||
        pycomplex3_register_loop(numpy, "conjugate", loop_complex3_unary, 2, 0, complex3_conjugate) ||
        pycomplex3_register_loop(numpy, "reciprocal", loop_complex3_unary, 2, 0, complex3_inv) ||
        pycomplex3_register_loop(numpy, "square", loop_complex3_unary, 2, 0, complex3_square) ||
        pycomplex3_register_loop(numpy, "exp", loop_complex3_unary, 2, 0, complex3_exp) ||
        pycomplex3_register_loop(numpy, "expm1", loop_complex3_unary, 2, 0, complex3_expm1) ||
        pycomplex3_register_loop(numpy, "log", loop_complex3_unary, 2, 0, complex3_log) ||
        pycomplex3_register_loop(numpy, "log1p", loop_complex3_unary, 2, 0, complex3_log1p) ||
        pycomplex3_register_loop(numpy, "sqrt", loop_complex3_unary, 2, 0, complex3_sqrt) ||
        pycomplex3_register_loop(numpy, "sin", loop_complex3_unary, 2, 0, complex3_sin) ||
        pycomplex3_register_loop(numpy, "cos", loop_complex3_unary, 2, 0, complex3_cos) ||
        pycomplex3_register_loop(numpy, "tan", loop_complex3_unary, 2, 0, complex3_tan) ||
        pycomplex3_register_loop(numpy, "arcsin", loop_complex3_unary, 2, 0, complex3_asin) ||
        pycomplex3_register_loop(numpy, "arccos", loop_complex3_unary, 2, 0, complex3_acos) ||
        pycomplex3_register_loop(numpy, "arctan", loop_complex3_unary, 2, 0, complex3_atan) ||
        pycomplex3_register_loop(numpy, "sinh", loop_complex3_unary, 2, 0, complex3_sinh) ||
        pycomplex3_register_loop(numpy, "cosh", loop_complex3_unary, 2, 0, complex3_cosh) ||
        pycomplex3_register_loop(numpy, "tanh", loop_complex3_unary, 2, 0, complex3_tanh) ||
        pycomplex3_register_loop(numpy, "arcsinh", loop_complex3_unary, 2, 0, complex3_asinh) ||
        pycomplex3_register_loop(numpy, "arccosh", loop_complex3_unary, 2, 0, complex3_acosh) ||
        pycomplex3_register_loop(numpy, "arctanh", loop_complex3_unary, 2, 0, complex3_atanh) ||
        pycomplex3_register_loop(numpy, "equal", loop_complex3_predicate, 3, 1, complex3_equal) ||
        pycomplex3_register_loop(numpy, "not_equal", loop_complex3_predicate, 3, 1, complex3_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex3", "complex3 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex3(void)
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
    PyComplex3_Type = PyType_FromModuleAndSpec(
        m, &pycomplex3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplex3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex3_arrfuncs);
    pycomplex3_arrfuncs.getitem = npycomplex3_getitem;
    pycomplex3_arrfuncs.setitem = npycomplex3_setitem;
    pycomplex3_arrfuncs.copyswapn = pycomplex3_copyswapn;
    pycomplex3_arrfuncs.copyswap = pycomplex3_copyswap;
    pycomplex3_arrfuncs.nonzero = pycomplex3_nonzero;
    pycomplex3_descr_proto.typeobj = (PyTypeObject *)PyComplex3_Type;
    Py_SET_TYPE((PyObject *)&pycomplex3_descr_proto, &PyArrayDescr_Type);
    npy_complex3 = PyArray_RegisterDataType(&pycomplex3_descr_proto);
    if (npy_complex3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplex3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplex3_Type, 8) < 0) {
        goto fail;
    }
    if (pycomplex3_register_casts() < 0 ||
            pycomplex3_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex3", PyComplex3_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
