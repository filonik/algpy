/*
 * Complex jets of order 5, C[e]/(e^6), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex_jet5.h.
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
#include "../../_grades.h"

#include "complex_jet5.h"

typedef struct {
    PyObject_HEAD
    complex_jet5 v;
} PyComplexJet5;

static PyObject *PyComplexJet5_Type = NULL;
static int npy_complex_jet5 = -1;

static PyObject *
PyComplexJet5_FromValue(PyTypeObject *type, complex_jet5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplexJet5 *p = (PyComplexJet5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX_JET5(v) \
    PyComplexJet5_FromValue((PyTypeObject *)PyComplexJet5_Type, (v))

static int
complex_jet5_from_object(PyObject *o, complex_jet5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplexJet5_Type)) {
        *out = ((PyComplexJet5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex_jet5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex_jet5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "y", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "x5", "y5", NULL};
    double d[12] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddd", kwlist,
                                     &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11])) {
        return NULL;
    }
    return PyComplexJet5_FromValue(type, complex_jet5_make(d));
}

static PyObject *
pycomplex_jet5_repr(PyObject *self)
{
    complex_jet5 v = ((PyComplexJet5 *)self)->v;
    char *s[12] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 12; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex_jet5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11]);
done:
    for (i = 0; i < 12; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex_jet5_hash(PyObject *self)
{
    complex_jet5 v = ((PyComplexJet5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_COMPLEX_JET5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (complex_jet5_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pycomplex_jet5_richcompare(PyObject *a, PyObject *b, int op)
{
    complex_jet5 x, y;
    int r;
    if ((r = complex_jet5_from_object(a, &x)) <= 0 ||
            (r = complex_jet5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(complex_jet5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!complex_jet5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define COMPLEX_JET5_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex_jet5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex_jet5 x, y;                                                       \
        int r;                                                               \
        if ((r = complex_jet5_from_object(a, &x)) <= 0 ||                        \
                (r = complex_jet5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX_JET5(expr(x, y));                                     \
    }

COMPLEX_JET5_BINARY(add, complex_jet5_add)
COMPLEX_JET5_BINARY(sub, complex_jet5_sub)
COMPLEX_JET5_BINARY(mul, complex_jet5_mul)
COMPLEX_JET5_BINARY(div, complex_jet5_div)

static PyObject *
pycomplex_jet5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex_jet5 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex_jet5 does not support modular power");
        return NULL;
    }
    if ((r = complex_jet5_from_object(a, &x)) <= 0 ||
            (r = complex_jet5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX_JET5(complex_jet5_rpow(x, y));
}

#define COMPLEX_JET5_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex_jet5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX_JET5(expr(((PyComplexJet5 *)self)->v));                  \
    }

COMPLEX_JET5_UNARY(pos, complex_jet5_pos)
COMPLEX_JET5_UNARY(neg, complex_jet5_neg)

static PyObject *
pycomplex_jet5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue((PyTypeObject *)cls, complex_jet5_zero());
}

static PyObject *
pycomplex_jet5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue((PyTypeObject *)cls, complex_jet5_one());
}

static int
pycomplex_jet5_bool(PyObject *self)
{
    complex_jet5 v = ((PyComplexJet5 *)self)->v;
    return vec12_any_nonzero(v);
}

#define COMPLEX_JET5_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex_jet5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX_JET5(op(((PyComplexJet5 *)self)->v));                   \
    }

COMPLEX_JET5_METHOD(conjugate, complex_jet5_conjugate)
COMPLEX_JET5_METHOD(exp, complex_jet5_exp)
COMPLEX_JET5_METHOD(log, complex_jet5_log)
COMPLEX_JET5_METHOD(sqrt, complex_jet5_sqrt)

static PyObject *
pycomplex_jet5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex_jet5_norm(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX_JET5(complex_jet5_inv(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {12, 12};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet5_lmat(((PyComplexJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {12, 12};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet5_rmat(((PyComplexJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {13};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet5_minpoly(((PyComplexJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[12];
    npy_intp dims[1] = {12};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex_jet5_eivals(((PyComplexJet5 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 12; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex_jet5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {12};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex_jet5);
    complex_jet5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex_jet5 *)PyArray_DATA(arr);
    complex_jet5_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex_jet5);
    complex_jet5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex_jet5 *)PyArray_DATA(arr);
    complex_jet5_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet5_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {12, 12};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex_jet5);
    complex_jet5 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex_jet5 *)PyArray_DATA(arr);
    complex_jet5_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex_jet5 v = ((PyComplexJet5 *)self)->v;
    return Py_BuildValue("(O(dddddddddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11]);
}

static PyObject *
pycomplex_jet5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_square(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_expm1(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_log1p(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_sin(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_cos(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_tan(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_asin(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_acos(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_atan(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_sinh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_cosh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_tanh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_asinh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_acosh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_atanh(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_re(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_im(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_real(PyObject *self, void *closure)
{
    return pycomplex_jet5_re(self, NULL);
}

static PyObject *
pycomplex_jet5_imag(PyObject *self, void *closure)
{
    return pycomplex_jet5_im(self, NULL);
}

static PyObject *pycomplex_jet5_grades(PyObject *self, void *closure);

static PyGetSetDef pycomplex_jet5_getset[] = {
    {"grade", pycomplex_jet5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycomplex_jet5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex_jet5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex_jet5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex_jet5_is_real(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex_jet5_to_re(((PyComplexJet5 *)self)->v));
}

static PyObject *
pycomplex_jet5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplexJet5_FromValue((PyTypeObject *)cls, complex_jet5_from_re(a));
}

static PyObject *
pycomplex_jet5_lpow(PyObject *self, PyObject *other)
{
    complex_jet5 b;
    int r;
    r = complex_jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet5_FromValue(Py_TYPE(self),
        complex_jet5_lpow(((PyComplexJet5 *)self)->v, b));
}

static PyObject *
pycomplex_jet5_rpow(PyObject *self, PyObject *other)
{
    complex_jet5 b;
    int r;
    r = complex_jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet5_FromValue(Py_TYPE(self),
        complex_jet5_rpow(((PyComplexJet5 *)self)->v, b));
}

static PyObject *
pycomplex_jet5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex_jet5 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex_jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet5_FromValue(Py_TYPE(self),
        complex_jet5_mix(((PyComplexJet5 *)self)->v, b, t));
}

static PyObject *
pycomplex_jet5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex_jet5, 12);
}

static PyObject *
pycomplex_jet5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex_jet5, 12);
}

static PyObject *
pycomplex_jet5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_COMPLEX_JET5(complex_jet5_ipow(((PyComplexJet5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycomplex_jet5_grade(PyObject *self, Py_ssize_t k)
{
    complex_jet5 v = ((PyComplexJet5 *)self)->v;
    switch (k) {
        case 0:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade0(v));
        case 1:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade1(v));
        case 2:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade2(v));
        case 3:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade3(v));
        case 4:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade4(v));
        case 5:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_grade5(v));
        default:
            return PyComplexJet5_FromValue(Py_TYPE(self), complex_jet5_zero());
    }
}

static PyObject *
pycomplex_jet5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycomplex_jet5_grade, 6);
}

static PyMethodDef pycomplex_jet5_methods[] = {
    {"to_float", (PyCFunction)pycomplex_jet5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex_jet5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex_jet5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex_jet5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex_jet5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex_jet5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycomplex_jet5_square, METH_NOARGS, "square."},
    {"expm1", pycomplex_jet5_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycomplex_jet5_log1p, METH_NOARGS, "log1p."},
    {"sin", pycomplex_jet5_sin, METH_NOARGS, "sin."},
    {"cos", pycomplex_jet5_cos, METH_NOARGS, "cos."},
    {"tan", pycomplex_jet5_tan, METH_NOARGS, "tan."},
    {"asin", pycomplex_jet5_asin, METH_NOARGS, "asin."},
    {"acos", pycomplex_jet5_acos, METH_NOARGS, "acos."},
    {"atan", pycomplex_jet5_atan, METH_NOARGS, "atan."},
    {"sinh", pycomplex_jet5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycomplex_jet5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycomplex_jet5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycomplex_jet5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycomplex_jet5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycomplex_jet5_atanh, METH_NOARGS, "atanh."},
    {"re", pycomplex_jet5_re, METH_NOARGS, "re."},
    {"im", pycomplex_jet5_im, METH_NOARGS, "im."},
    {"arcsin", pycomplex_jet5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex_jet5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex_jet5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex_jet5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex_jet5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex_jet5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex_jet5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex_jet5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex_jet5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycomplex_jet5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pycomplex_jet5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycomplex_jet5_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pycomplex_jet5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycomplex_jet5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycomplex_jet5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycomplex_jet5_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex_jet5_conjugate, METH_NOARGS, "Complex conjugate i -> -i of every coefficient."},
    {"lmat", pycomplex_jet5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex_jet5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex_jet5_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pycomplex_jet5_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex_jet5_exp, METH_NOARGS, "exp."},
    {"log", pycomplex_jet5_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex_jet5_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex_jet5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex_jet5_length(PyObject *Py_UNUSED(self))
{
    return 12;
}

static PyObject *
pycomplex_jet5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 12) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplexJet5 *)self)->v.data[i]);
}

static PyObject *
pycomplex_jet5_subscript(PyObject *self, PyObject *key)
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
        i += 12;
    }
    return pycomplex_jet5_item(self, i);
}

static int
pycomplex_jet5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex_jet5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex_jet5_slots[] = {
    {Py_tp_dealloc, pycomplex_jet5_dealloc},
    {Py_tp_traverse, pycomplex_jet5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex_jet5_repr},
    {Py_tp_str, pycomplex_jet5_repr},
    {Py_tp_hash, pycomplex_jet5_hash},
    {Py_tp_doc, (void *)"Complex jet of order 5 in C[e]/(e^6) with 12 float64 components"},
    {Py_tp_richcompare, pycomplex_jet5_richcompare},
    {Py_tp_methods, pycomplex_jet5_methods},
    {Py_tp_getset, pycomplex_jet5_getset},
    {Py_sq_length, pycomplex_jet5_length},
    {Py_sq_item, pycomplex_jet5_item},
    {Py_mp_subscript, pycomplex_jet5_subscript},
    {Py_tp_new, pycomplex_jet5_new},
    {Py_nb_add, pycomplex_jet5_add},
    {Py_nb_subtract, pycomplex_jet5_sub},
    {Py_nb_multiply, pycomplex_jet5_mul},
    {Py_nb_true_divide, pycomplex_jet5_div},
    {Py_nb_power, (void *)pycomplex_jet5_pow},
    {Py_nb_positive, pycomplex_jet5_pos},
    {Py_nb_negative, pycomplex_jet5_neg},
    {Py_nb_bool, pycomplex_jet5_bool},
    {0, NULL}
};

static PyType_Spec pycomplex_jet5_spec = {
    .name = "algpy.algebras.jet.complex_jet5",
    .basicsize = sizeof(PyComplexJet5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex_jet5_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 12; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npycomplex_jet5_getitem(void *data, void *arr)
{
    complex_jet5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_COMPLEX_JET5(v);
}

static int
npycomplex_jet5_setitem(PyObject *item, void *data, void *arr)
{
    complex_jet5 v;
    int r = complex_jet5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex_jet5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 12) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex_jet5 needs exactly 12 components");
            return -1;
        }
        for (int i = 0; i < 12; i++) {
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
pycomplex_jet5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex_jet5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycomplex_jet5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex_jet5_copyswapn(dst, sizeof(complex_jet5), src, sizeof(complex_jet5), 1,
                         swap, arr);
}

static npy_bool
pycomplex_jet5_nonzero(void *data, void *arr)
{
    complex_jet5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec12_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex_jet5_arrfuncs;
typedef struct { char c; complex_jet5 v; } pycomplex_jet5_align_test;

static PyArray_DescrProto pycomplex_jet5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex_jet5), offsetof(pycomplex_jet5_align_test, v), 0, 0, 0,
    &pycomplex_jet5_arrfuncs,
};

#define COMPLEX_JET5_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycomplex_jet5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex_jet5 *dst = (complex_jet5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex_jet5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define COMPLEX_JET5_REAL_TYPES(X) \
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

COMPLEX_JET5_REAL_TYPES(COMPLEX_JET5_CAST_FROM)

static int
pycomplex_jet5_register_casts(void)
{
#define REGISTER_COMPLEX_JET5_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex_jet5,                    \
                                     pycomplex_jet5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex_jet5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX_JET5_REAL_TYPES(REGISTER_COMPLEX_JET5_CAST)
#undef REGISTER_COMPLEX_JET5_CAST
    return 0;
}

typedef complex_jet5 (*complex_jet5_unary_op)(complex_jet5);
typedef complex_jet5 (*complex_jet5_binary_op)(complex_jet5, complex_jet5);
typedef int (*complex_jet5_predicate)(complex_jet5, complex_jet5);

static void
loop_complex_jet5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex_jet5_unary_op op = (complex_jet5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex_jet5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex_jet5_binary_op op = (complex_jet5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet5 a, b, r;
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
loop_complex_jet5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex_jet5_predicate op = (complex_jet5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex_jet5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex_jet5;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex_jet5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex_jet5_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex_jet5_register_loop(numpy, "add", loop_complex_jet5_binary, 3, 0, complex_jet5_add) ||
        pycomplex_jet5_register_loop(numpy, "subtract", loop_complex_jet5_binary, 3, 0, complex_jet5_sub) ||
        pycomplex_jet5_register_loop(numpy, "multiply", loop_complex_jet5_binary, 3, 0, complex_jet5_mul) ||
        pycomplex_jet5_register_loop(numpy, "divide", loop_complex_jet5_binary, 3, 0, complex_jet5_div) ||
        pycomplex_jet5_register_loop(numpy, "power", loop_complex_jet5_binary, 3, 0, complex_jet5_rpow) ||
        pycomplex_jet5_register_loop(numpy, "positive", loop_complex_jet5_unary, 2, 0, complex_jet5_pos) ||
        pycomplex_jet5_register_loop(numpy, "negative", loop_complex_jet5_unary, 2, 0, complex_jet5_neg) ||
        pycomplex_jet5_register_loop(numpy, "conjugate", loop_complex_jet5_unary, 2, 0, complex_jet5_conjugate) ||
        pycomplex_jet5_register_loop(numpy, "reciprocal", loop_complex_jet5_unary, 2, 0, complex_jet5_inv) ||
        pycomplex_jet5_register_loop(numpy, "square", loop_complex_jet5_unary, 2, 0, complex_jet5_square) ||
        pycomplex_jet5_register_loop(numpy, "exp", loop_complex_jet5_unary, 2, 0, complex_jet5_exp) ||
        pycomplex_jet5_register_loop(numpy, "expm1", loop_complex_jet5_unary, 2, 0, complex_jet5_expm1) ||
        pycomplex_jet5_register_loop(numpy, "log", loop_complex_jet5_unary, 2, 0, complex_jet5_log) ||
        pycomplex_jet5_register_loop(numpy, "log1p", loop_complex_jet5_unary, 2, 0, complex_jet5_log1p) ||
        pycomplex_jet5_register_loop(numpy, "sqrt", loop_complex_jet5_unary, 2, 0, complex_jet5_sqrt) ||
        pycomplex_jet5_register_loop(numpy, "sin", loop_complex_jet5_unary, 2, 0, complex_jet5_sin) ||
        pycomplex_jet5_register_loop(numpy, "cos", loop_complex_jet5_unary, 2, 0, complex_jet5_cos) ||
        pycomplex_jet5_register_loop(numpy, "tan", loop_complex_jet5_unary, 2, 0, complex_jet5_tan) ||
        pycomplex_jet5_register_loop(numpy, "arcsin", loop_complex_jet5_unary, 2, 0, complex_jet5_asin) ||
        pycomplex_jet5_register_loop(numpy, "arccos", loop_complex_jet5_unary, 2, 0, complex_jet5_acos) ||
        pycomplex_jet5_register_loop(numpy, "arctan", loop_complex_jet5_unary, 2, 0, complex_jet5_atan) ||
        pycomplex_jet5_register_loop(numpy, "sinh", loop_complex_jet5_unary, 2, 0, complex_jet5_sinh) ||
        pycomplex_jet5_register_loop(numpy, "cosh", loop_complex_jet5_unary, 2, 0, complex_jet5_cosh) ||
        pycomplex_jet5_register_loop(numpy, "tanh", loop_complex_jet5_unary, 2, 0, complex_jet5_tanh) ||
        pycomplex_jet5_register_loop(numpy, "arcsinh", loop_complex_jet5_unary, 2, 0, complex_jet5_asinh) ||
        pycomplex_jet5_register_loop(numpy, "arccosh", loop_complex_jet5_unary, 2, 0, complex_jet5_acosh) ||
        pycomplex_jet5_register_loop(numpy, "arctanh", loop_complex_jet5_unary, 2, 0, complex_jet5_atanh) ||
        pycomplex_jet5_register_loop(numpy, "equal", loop_complex_jet5_predicate, 3, 1, complex_jet5_equal) ||
        pycomplex_jet5_register_loop(numpy, "not_equal", loop_complex_jet5_predicate, 3, 1, complex_jet5_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex_jet5", "Complex jet type of order 5 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex_jet5(void)
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
    PyComplexJet5_Type = PyType_FromModuleAndSpec(
        m, &pycomplex_jet5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplexJet5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex_jet5_arrfuncs);
    pycomplex_jet5_arrfuncs.getitem = npycomplex_jet5_getitem;
    pycomplex_jet5_arrfuncs.setitem = npycomplex_jet5_setitem;
    pycomplex_jet5_arrfuncs.copyswapn = pycomplex_jet5_copyswapn;
    pycomplex_jet5_arrfuncs.copyswap = pycomplex_jet5_copyswap;
    pycomplex_jet5_arrfuncs.nonzero = pycomplex_jet5_nonzero;
    pycomplex_jet5_descr_proto.typeobj = (PyTypeObject *)PyComplexJet5_Type;
    Py_SET_TYPE((PyObject *)&pycomplex_jet5_descr_proto, &PyArrayDescr_Type);
    npy_complex_jet5 = PyArray_RegisterDataType(&pycomplex_jet5_descr_proto);
    if (npy_complex_jet5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex_jet5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplexJet5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplexJet5_Type, 12) < 0) {
        goto fail;
    }
    if (pycomplex_jet5_register_casts() < 0 ||
            pycomplex_jet5_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex_jet5", PyComplexJet5_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
