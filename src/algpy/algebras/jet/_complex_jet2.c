/*
 * Complex jets of order 2, C[e]/(e^3), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex_jet2.h.
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

#include "complex_jet2.h"

typedef struct {
    PyObject_HEAD
    complex_jet2 v;
} PyComplexJet2;

static PyObject *PyComplexJet2_Type = NULL;
static int npy_complex_jet2 = -1;

static PyObject *
PyComplexJet2_FromValue(PyTypeObject *type, complex_jet2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplexJet2 *p = (PyComplexJet2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX_JET2(v) \
    PyComplexJet2_FromValue((PyTypeObject *)PyComplexJet2_Type, (v))

static int
complex_jet2_from_object(PyObject *o, complex_jet2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplexJet2_Type)) {
        *out = ((PyComplexJet2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex_jet2_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex_jet2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "y", "x1", "y1", "x2", "y2", NULL};
    double d[6] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddd", kwlist,
                                     &d[0], &d[1], &d[2], &d[3], &d[4], &d[5])) {
        return NULL;
    }
    return PyComplexJet2_FromValue(type, complex_jet2_make(d));
}

static PyObject *
pycomplex_jet2_repr(PyObject *self)
{
    complex_jet2 v = ((PyComplexJet2 *)self)->v;
    char *s[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 6; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex_jet2(%s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5]);
done:
    for (i = 0; i < 6; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex_jet2_hash(PyObject *self)
{
    complex_jet2 v = ((PyComplexJet2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_COMPLEX_JET2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (complex_jet2_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pycomplex_jet2_richcompare(PyObject *a, PyObject *b, int op)
{
    complex_jet2 x, y;
    int r;
    if ((r = complex_jet2_from_object(a, &x)) <= 0 ||
            (r = complex_jet2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(complex_jet2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!complex_jet2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define COMPLEX_JET2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex_jet2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex_jet2 x, y;                                                       \
        int r;                                                               \
        if ((r = complex_jet2_from_object(a, &x)) <= 0 ||                        \
                (r = complex_jet2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX_JET2(expr(x, y));                                     \
    }

COMPLEX_JET2_BINARY(add, complex_jet2_add)
COMPLEX_JET2_BINARY(sub, complex_jet2_sub)
COMPLEX_JET2_BINARY(mul, complex_jet2_mul)
COMPLEX_JET2_BINARY(div, complex_jet2_div)

static PyObject *
pycomplex_jet2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex_jet2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex_jet2 does not support modular power");
        return NULL;
    }
    if ((r = complex_jet2_from_object(a, &x)) <= 0 ||
            (r = complex_jet2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX_JET2(complex_jet2_rpow(x, y));
}

#define COMPLEX_JET2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex_jet2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX_JET2(expr(((PyComplexJet2 *)self)->v));                  \
    }

COMPLEX_JET2_UNARY(pos, complex_jet2_pos)
COMPLEX_JET2_UNARY(neg, complex_jet2_neg)

static PyObject *
pycomplex_jet2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue((PyTypeObject *)cls, complex_jet2_zero());
}

static PyObject *
pycomplex_jet2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue((PyTypeObject *)cls, complex_jet2_one());
}

static int
pycomplex_jet2_bool(PyObject *self)
{
    complex_jet2 v = ((PyComplexJet2 *)self)->v;
    return vec6_any_nonzero(v);
}

#define COMPLEX_JET2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex_jet2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX_JET2(op(((PyComplexJet2 *)self)->v));                   \
    }

COMPLEX_JET2_METHOD(conjugate, complex_jet2_conjugate)
COMPLEX_JET2_METHOD(exp, complex_jet2_exp)
COMPLEX_JET2_METHOD(log, complex_jet2_log)
COMPLEX_JET2_METHOD(sqrt, complex_jet2_sqrt)

static PyObject *
pycomplex_jet2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex_jet2_norm(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX_JET2(complex_jet2_inv(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet2_lmat(((PyComplexJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet2_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet2_rmat(((PyComplexJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {7};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex_jet2_minpoly(((PyComplexJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex_jet2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[6];
    npy_intp dims[1] = {6};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex_jet2_eivals(((PyComplexJet2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 6; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex_jet2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {6};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex_jet2);
    complex_jet2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex_jet2 *)PyArray_DATA(arr);
    complex_jet2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet2_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex_jet2);
    complex_jet2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex_jet2 *)PyArray_DATA(arr);
    complex_jet2_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex_jet2);
    complex_jet2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex_jet2 *)PyArray_DATA(arr);
    complex_jet2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex_jet2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex_jet2 v = ((PyComplexJet2 *)self)->v;
    return Py_BuildValue("(O(dddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5]);
}

static PyObject *
pycomplex_jet2_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_square(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_expm1(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_log1p(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_sin(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_cos(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_tan(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_asin(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_acos(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_atan(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_sinh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_cosh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_tanh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_asinh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_acosh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_atanh(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_re(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_im(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_real(PyObject *self, void *closure)
{
    return pycomplex_jet2_re(self, NULL);
}

static PyObject *
pycomplex_jet2_imag(PyObject *self, void *closure)
{
    return pycomplex_jet2_im(self, NULL);
}

static PyObject *pycomplex_jet2_grades(PyObject *self, void *closure);

static PyGetSetDef pycomplex_jet2_getset[] = {
    {"grade", pycomplex_jet2_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycomplex_jet2_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex_jet2_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex_jet2_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex_jet2_is_real(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex_jet2_to_re(((PyComplexJet2 *)self)->v));
}

static PyObject *
pycomplex_jet2_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplexJet2_FromValue((PyTypeObject *)cls, complex_jet2_from_re(a));
}

static PyObject *
pycomplex_jet2_lpow(PyObject *self, PyObject *other)
{
    complex_jet2 b;
    int r;
    r = complex_jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet2_FromValue(Py_TYPE(self),
        complex_jet2_lpow(((PyComplexJet2 *)self)->v, b));
}

static PyObject *
pycomplex_jet2_rpow(PyObject *self, PyObject *other)
{
    complex_jet2 b;
    int r;
    r = complex_jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet2_FromValue(Py_TYPE(self),
        complex_jet2_rpow(((PyComplexJet2 *)self)->v, b));
}

static PyObject *
pycomplex_jet2_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex_jet2 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex_jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplexJet2_FromValue(Py_TYPE(self),
        complex_jet2_mix(((PyComplexJet2 *)self)->v, b, t));
}

static PyObject *
pycomplex_jet2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex_jet2, 6);
}

static PyObject *
pycomplex_jet2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex_jet2, 6);
}

static PyObject *
pycomplex_jet2_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_COMPLEX_JET2(complex_jet2_ipow(((PyComplexJet2 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycomplex_jet2_grade(PyObject *self, Py_ssize_t k)
{
    complex_jet2 v = ((PyComplexJet2 *)self)->v;
    switch (k) {
        case 0:
            return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_grade0(v));
        case 1:
            return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_grade1(v));
        case 2:
            return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_grade2(v));
        default:
            return PyComplexJet2_FromValue(Py_TYPE(self), complex_jet2_zero());
    }
}

static PyObject *
pycomplex_jet2_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycomplex_jet2_grade, 3);
}

static PyMethodDef pycomplex_jet2_methods[] = {
    {"to_float", (PyCFunction)pycomplex_jet2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex_jet2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex_jet2_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex_jet2_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex_jet2_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex_jet2_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycomplex_jet2_square, METH_NOARGS, "square."},
    {"expm1", pycomplex_jet2_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycomplex_jet2_log1p, METH_NOARGS, "log1p."},
    {"sin", pycomplex_jet2_sin, METH_NOARGS, "sin."},
    {"cos", pycomplex_jet2_cos, METH_NOARGS, "cos."},
    {"tan", pycomplex_jet2_tan, METH_NOARGS, "tan."},
    {"asin", pycomplex_jet2_asin, METH_NOARGS, "asin."},
    {"acos", pycomplex_jet2_acos, METH_NOARGS, "acos."},
    {"atan", pycomplex_jet2_atan, METH_NOARGS, "atan."},
    {"sinh", pycomplex_jet2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycomplex_jet2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycomplex_jet2_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycomplex_jet2_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycomplex_jet2_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycomplex_jet2_atanh, METH_NOARGS, "atanh."},
    {"re", pycomplex_jet2_re, METH_NOARGS, "re."},
    {"im", pycomplex_jet2_im, METH_NOARGS, "im."},
    {"arcsin", pycomplex_jet2_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex_jet2_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex_jet2_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex_jet2_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex_jet2_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex_jet2_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex_jet2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex_jet2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex_jet2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycomplex_jet2_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pycomplex_jet2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycomplex_jet2_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pycomplex_jet2_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycomplex_jet2_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycomplex_jet2_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycomplex_jet2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex_jet2_conjugate, METH_NOARGS, "Complex conjugate i -> -i of every coefficient."},
    {"lmat", pycomplex_jet2_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex_jet2_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex_jet2_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pycomplex_jet2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex_jet2_exp, METH_NOARGS, "exp."},
    {"log", pycomplex_jet2_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex_jet2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex_jet2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex_jet2_length(PyObject *Py_UNUSED(self))
{
    return 6;
}

static PyObject *
pycomplex_jet2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 6) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplexJet2 *)self)->v.data[i]);
}

static PyObject *
pycomplex_jet2_subscript(PyObject *self, PyObject *key)
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
        i += 6;
    }
    return pycomplex_jet2_item(self, i);
}

static int
pycomplex_jet2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex_jet2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex_jet2_slots[] = {
    {Py_tp_dealloc, pycomplex_jet2_dealloc},
    {Py_tp_traverse, pycomplex_jet2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex_jet2_repr},
    {Py_tp_str, pycomplex_jet2_repr},
    {Py_tp_hash, pycomplex_jet2_hash},
    {Py_tp_doc, (void *)"Complex jet of order 2 in C[e]/(e^3) with 6 float64 components"},
    {Py_tp_richcompare, pycomplex_jet2_richcompare},
    {Py_tp_methods, pycomplex_jet2_methods},
    {Py_tp_getset, pycomplex_jet2_getset},
    {Py_sq_length, pycomplex_jet2_length},
    {Py_sq_item, pycomplex_jet2_item},
    {Py_mp_subscript, pycomplex_jet2_subscript},
    {Py_tp_new, pycomplex_jet2_new},
    {Py_nb_add, pycomplex_jet2_add},
    {Py_nb_subtract, pycomplex_jet2_sub},
    {Py_nb_multiply, pycomplex_jet2_mul},
    {Py_nb_true_divide, pycomplex_jet2_div},
    {Py_nb_power, (void *)pycomplex_jet2_pow},
    {Py_nb_positive, pycomplex_jet2_pos},
    {Py_nb_negative, pycomplex_jet2_neg},
    {Py_nb_bool, pycomplex_jet2_bool},
    {0, NULL}
};

static PyType_Spec pycomplex_jet2_spec = {
    .name = "algpy.algebras.jet.complex_jet2",
    .basicsize = sizeof(PyComplexJet2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex_jet2_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npycomplex_jet2_getitem(void *data, void *arr)
{
    complex_jet2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_COMPLEX_JET2(v);
}

static int
npycomplex_jet2_setitem(PyObject *item, void *data, void *arr)
{
    complex_jet2 v;
    int r = complex_jet2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex_jet2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 6) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex_jet2 needs exactly 6 components");
            return -1;
        }
        for (int i = 0; i < 6; i++) {
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
pycomplex_jet2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex_jet2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycomplex_jet2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex_jet2_copyswapn(dst, sizeof(complex_jet2), src, sizeof(complex_jet2), 1,
                         swap, arr);
}

static npy_bool
pycomplex_jet2_nonzero(void *data, void *arr)
{
    complex_jet2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec6_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex_jet2_arrfuncs;
typedef struct { char c; complex_jet2 v; } pycomplex_jet2_align_test;

static PyArray_DescrProto pycomplex_jet2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex_jet2), offsetof(pycomplex_jet2_align_test, v), 0, 0, 0,
    &pycomplex_jet2_arrfuncs,
};

#define COMPLEX_JET2_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycomplex_jet2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex_jet2 *dst = (complex_jet2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex_jet2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define COMPLEX_JET2_REAL_TYPES(X) \
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

COMPLEX_JET2_REAL_TYPES(COMPLEX_JET2_CAST_FROM)

static int
pycomplex_jet2_register_casts(void)
{
#define REGISTER_COMPLEX_JET2_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex_jet2,                    \
                                     pycomplex_jet2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex_jet2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX_JET2_REAL_TYPES(REGISTER_COMPLEX_JET2_CAST)
#undef REGISTER_COMPLEX_JET2_CAST
    return 0;
}

typedef complex_jet2 (*complex_jet2_unary_op)(complex_jet2);
typedef complex_jet2 (*complex_jet2_binary_op)(complex_jet2, complex_jet2);
typedef int (*complex_jet2_predicate)(complex_jet2, complex_jet2);

static void
loop_complex_jet2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex_jet2_unary_op op = (complex_jet2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex_jet2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex_jet2_binary_op op = (complex_jet2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet2 a, b, r;
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
loop_complex_jet2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex_jet2_predicate op = (complex_jet2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex_jet2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex_jet2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex_jet2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex_jet2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex_jet2_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex_jet2_register_loop(numpy, "add", loop_complex_jet2_binary, 3, 0, complex_jet2_add) ||
        pycomplex_jet2_register_loop(numpy, "subtract", loop_complex_jet2_binary, 3, 0, complex_jet2_sub) ||
        pycomplex_jet2_register_loop(numpy, "multiply", loop_complex_jet2_binary, 3, 0, complex_jet2_mul) ||
        pycomplex_jet2_register_loop(numpy, "divide", loop_complex_jet2_binary, 3, 0, complex_jet2_div) ||
        pycomplex_jet2_register_loop(numpy, "power", loop_complex_jet2_binary, 3, 0, complex_jet2_rpow) ||
        pycomplex_jet2_register_loop(numpy, "positive", loop_complex_jet2_unary, 2, 0, complex_jet2_pos) ||
        pycomplex_jet2_register_loop(numpy, "negative", loop_complex_jet2_unary, 2, 0, complex_jet2_neg) ||
        pycomplex_jet2_register_loop(numpy, "conjugate", loop_complex_jet2_unary, 2, 0, complex_jet2_conjugate) ||
        pycomplex_jet2_register_loop(numpy, "reciprocal", loop_complex_jet2_unary, 2, 0, complex_jet2_inv) ||
        pycomplex_jet2_register_loop(numpy, "square", loop_complex_jet2_unary, 2, 0, complex_jet2_square) ||
        pycomplex_jet2_register_loop(numpy, "exp", loop_complex_jet2_unary, 2, 0, complex_jet2_exp) ||
        pycomplex_jet2_register_loop(numpy, "expm1", loop_complex_jet2_unary, 2, 0, complex_jet2_expm1) ||
        pycomplex_jet2_register_loop(numpy, "log", loop_complex_jet2_unary, 2, 0, complex_jet2_log) ||
        pycomplex_jet2_register_loop(numpy, "log1p", loop_complex_jet2_unary, 2, 0, complex_jet2_log1p) ||
        pycomplex_jet2_register_loop(numpy, "sqrt", loop_complex_jet2_unary, 2, 0, complex_jet2_sqrt) ||
        pycomplex_jet2_register_loop(numpy, "sin", loop_complex_jet2_unary, 2, 0, complex_jet2_sin) ||
        pycomplex_jet2_register_loop(numpy, "cos", loop_complex_jet2_unary, 2, 0, complex_jet2_cos) ||
        pycomplex_jet2_register_loop(numpy, "tan", loop_complex_jet2_unary, 2, 0, complex_jet2_tan) ||
        pycomplex_jet2_register_loop(numpy, "arcsin", loop_complex_jet2_unary, 2, 0, complex_jet2_asin) ||
        pycomplex_jet2_register_loop(numpy, "arccos", loop_complex_jet2_unary, 2, 0, complex_jet2_acos) ||
        pycomplex_jet2_register_loop(numpy, "arctan", loop_complex_jet2_unary, 2, 0, complex_jet2_atan) ||
        pycomplex_jet2_register_loop(numpy, "sinh", loop_complex_jet2_unary, 2, 0, complex_jet2_sinh) ||
        pycomplex_jet2_register_loop(numpy, "cosh", loop_complex_jet2_unary, 2, 0, complex_jet2_cosh) ||
        pycomplex_jet2_register_loop(numpy, "tanh", loop_complex_jet2_unary, 2, 0, complex_jet2_tanh) ||
        pycomplex_jet2_register_loop(numpy, "arcsinh", loop_complex_jet2_unary, 2, 0, complex_jet2_asinh) ||
        pycomplex_jet2_register_loop(numpy, "arccosh", loop_complex_jet2_unary, 2, 0, complex_jet2_acosh) ||
        pycomplex_jet2_register_loop(numpy, "arctanh", loop_complex_jet2_unary, 2, 0, complex_jet2_atanh) ||
        pycomplex_jet2_register_loop(numpy, "equal", loop_complex_jet2_predicate, 3, 1, complex_jet2_equal) ||
        pycomplex_jet2_register_loop(numpy, "not_equal", loop_complex_jet2_predicate, 3, 1, complex_jet2_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex_jet2", "Complex jet type of order 2 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex_jet2(void)
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
    PyComplexJet2_Type = PyType_FromModuleAndSpec(
        m, &pycomplex_jet2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplexJet2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex_jet2_arrfuncs);
    pycomplex_jet2_arrfuncs.getitem = npycomplex_jet2_getitem;
    pycomplex_jet2_arrfuncs.setitem = npycomplex_jet2_setitem;
    pycomplex_jet2_arrfuncs.copyswapn = pycomplex_jet2_copyswapn;
    pycomplex_jet2_arrfuncs.copyswap = pycomplex_jet2_copyswap;
    pycomplex_jet2_arrfuncs.nonzero = pycomplex_jet2_nonzero;
    pycomplex_jet2_descr_proto.typeobj = (PyTypeObject *)PyComplexJet2_Type;
    Py_SET_TYPE((PyObject *)&pycomplex_jet2_descr_proto, &PyArrayDescr_Type);
    npy_complex_jet2 = PyArray_RegisterDataType(&pycomplex_jet2_descr_proto);
    if (npy_complex_jet2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex_jet2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplexJet2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplexJet2_Type, 6) < 0) {
        goto fail;
    }
    if (pycomplex_jet2_register_casts() < 0 ||
            pycomplex_jet2_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex_jet2", PyComplexJet2_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
