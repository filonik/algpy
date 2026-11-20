/*
 * complex5 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex5.h.
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

#include "complex5.h"

typedef struct {
    PyObject_HEAD
    complex5 v;
} PyComplex5;

static PyObject *PyComplex5_Type = NULL;
static int npy_complex5 = -1;

static PyObject *
PyComplex5_FromValue(PyTypeObject *type, complex5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplex5 *p = (PyComplex5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX5(v) \
    PyComplex5_FromValue((PyTypeObject *)PyComplex5_Type, (v))

static int
complex5_from_object(PyObject *o, complex5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplex5_Type)) {
        *out = ((PyComplex5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyComplex5_FromValue(type, complex5_make(d));
}

static PyObject *
pycomplex5_repr(PyObject *self)
{
    complex5 v = ((PyComplex5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex5_hash(PyObject *self)
{
    complex5 v = ((PyComplex5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_COMPLEX5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (complex5_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddddddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pycomplex5_richcompare(PyObject *a, PyObject *b, int op)
{
    complex5 x, y;
    int r;
    if ((r = complex5_from_object(a, &x)) <= 0 ||
            (r = complex5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(complex5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!complex5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define COMPLEX5_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex5 x, y;                                                       \
        int r;                                                               \
        if ((r = complex5_from_object(a, &x)) <= 0 ||                        \
                (r = complex5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX5(expr(x, y));                                     \
    }

COMPLEX5_BINARY(add, complex5_add)
COMPLEX5_BINARY(sub, complex5_sub)
COMPLEX5_BINARY(mul, complex5_mul)
COMPLEX5_BINARY(div, complex5_div)

static PyObject *
pycomplex5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex5 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex5 does not support modular power");
        return NULL;
    }
    if ((r = complex5_from_object(a, &x)) <= 0 ||
            (r = complex5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX5(complex5_rpow(x, y));
}

#define COMPLEX5_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX5(expr(((PyComplex5 *)self)->v));                  \
    }

COMPLEX5_UNARY(pos, complex5_pos)
COMPLEX5_UNARY(neg, complex5_neg)

static PyObject *
pycomplex5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue((PyTypeObject *)cls, complex5_zero());
}

static PyObject *
pycomplex5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue((PyTypeObject *)cls, complex5_one());
}

static int
pycomplex5_bool(PyObject *self)
{
    complex5 v = ((PyComplex5 *)self)->v;
    return vec32_any_nonzero(v);
}

#define COMPLEX5_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX5(op(((PyComplex5 *)self)->v));                   \
    }

COMPLEX5_METHOD(conjugate, complex5_conjugate)
COMPLEX5_METHOD(exp, complex5_exp)
COMPLEX5_METHOD(log, complex5_log)
COMPLEX5_METHOD(sqrt, complex5_sqrt)

static PyObject *
pycomplex5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex5_norm(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX5(complex5_inv(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex5_lmat(((PyComplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex5_rmat(((PyComplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {33};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex5_minpoly(((PyComplex5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[32];
    npy_intp dims[1] = {32};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex5_eivals(((PyComplex5 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 32; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex5);
    complex5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex5 *)PyArray_DATA(arr);
    complex5_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex5);
    complex5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex5 *)PyArray_DATA(arr);
    complex5_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex5_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex5);
    complex5 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex5 *)PyArray_DATA(arr);
    complex5_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex5 v = ((PyComplex5 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pycomplex5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_square(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_expm1(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_log1p(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_sin(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_cos(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_tan(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_asin(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_acos(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_atan(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_sinh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_cosh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_tanh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_asinh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_acosh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_atanh(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_re(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyComplex5_FromValue(Py_TYPE(self), complex5_im(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_real(PyObject *self, void *closure)
{
    return pycomplex5_re(self, NULL);
}

static PyObject *
pycomplex5_imag(PyObject *self, void *closure)
{
    return pycomplex5_im(self, NULL);
}

static PyObject *pycomplex5_grades(PyObject *self, void *closure);

static PyGetSetDef pycomplex5_getset[] = {
    {"grade", pycomplex5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycomplex5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycomplex5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycomplex5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(complex5_is_real(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex5_to_re(((PyComplex5 *)self)->v));
}

static PyObject *
pycomplex5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyComplex5_FromValue((PyTypeObject *)cls, complex5_from_re(a));
}

static PyObject *
pycomplex5_lpow(PyObject *self, PyObject *other)
{
    complex5 b;
    int r;
    r = complex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex5_FromValue(Py_TYPE(self),
        complex5_lpow(((PyComplex5 *)self)->v, b));
}

static PyObject *
pycomplex5_rpow(PyObject *self, PyObject *other)
{
    complex5 b;
    int r;
    r = complex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex5_FromValue(Py_TYPE(self),
        complex5_rpow(((PyComplex5 *)self)->v, b));
}

static PyObject *
pycomplex5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    complex5 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = complex5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyComplex5_FromValue(Py_TYPE(self),
        complex5_mix(((PyComplex5 *)self)->v, b, t));
}

static PyObject *
pycomplex5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_complex5, 32);
}

static PyObject *
pycomplex5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_complex5, 32);
}

static PyObject *
pycomplex5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_COMPLEX5(complex5_ipow(((PyComplex5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycomplex5_grade(PyObject *self, Py_ssize_t k)
{
    complex5 v = ((PyComplex5 *)self)->v;
    switch (k) {
        case 0:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade0(v));
        case 1:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade1(v));
        case 2:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade2(v));
        case 3:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade3(v));
        case 4:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade4(v));
        case 5:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_grade5(v));
        default:
            return PyComplex5_FromValue(Py_TYPE(self), complex5_zero());
    }
}

static PyObject *
pycomplex5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycomplex5_grade, 6);
}

static PyMethodDef pycomplex5_methods[] = {
    {"to_float", (PyCFunction)pycomplex5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycomplex5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycomplex5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycomplex5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycomplex5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycomplex5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycomplex5_square, METH_NOARGS, "square."},
    {"expm1", pycomplex5_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycomplex5_log1p, METH_NOARGS, "log1p."},
    {"sin", pycomplex5_sin, METH_NOARGS, "sin."},
    {"cos", pycomplex5_cos, METH_NOARGS, "cos."},
    {"tan", pycomplex5_tan, METH_NOARGS, "tan."},
    {"asin", pycomplex5_asin, METH_NOARGS, "asin."},
    {"acos", pycomplex5_acos, METH_NOARGS, "acos."},
    {"atan", pycomplex5_atan, METH_NOARGS, "atan."},
    {"sinh", pycomplex5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycomplex5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycomplex5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycomplex5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycomplex5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycomplex5_atanh, METH_NOARGS, "atanh."},
    {"re", pycomplex5_re, METH_NOARGS, "re."},
    {"im", pycomplex5_im, METH_NOARGS, "im."},
    {"arcsin", pycomplex5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycomplex5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycomplex5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycomplex5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycomplex5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycomplex5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycomplex5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycomplex5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pycomplex5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycomplex5_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pycomplex5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycomplex5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycomplex5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycomplex5_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex5_conjugate, METH_NOARGS, "Conjugation of the outer complex generator."},
    {"lmat", pycomplex5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycomplex5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycomplex5_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pycomplex5_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex5_exp, METH_NOARGS, "exp."},
    {"log", pycomplex5_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex5_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex5_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pycomplex5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplex5 *)self)->v.data[i]);
}

static PyObject *
pycomplex5_subscript(PyObject *self, PyObject *key)
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
        i += 32;
    }
    return pycomplex5_item(self, i);
}

static int
pycomplex5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex5_slots[] = {
    {Py_tp_dealloc, pycomplex5_dealloc},
    {Py_tp_traverse, pycomplex5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex5_repr},
    {Py_tp_str, pycomplex5_repr},
    {Py_tp_hash, pycomplex5_hash},
    {Py_tp_doc, (void *)"Fivefold complex number with 32 float64 coefficients"},
    {Py_tp_richcompare, pycomplex5_richcompare},
    {Py_tp_methods, pycomplex5_methods},
    {Py_tp_getset, pycomplex5_getset},
    {Py_sq_length, pycomplex5_length},
    {Py_sq_item, pycomplex5_item},
    {Py_mp_subscript, pycomplex5_subscript},
    {Py_tp_new, pycomplex5_new},
    {Py_nb_add, pycomplex5_add},
    {Py_nb_subtract, pycomplex5_sub},
    {Py_nb_multiply, pycomplex5_mul},
    {Py_nb_true_divide, pycomplex5_div},
    {Py_nb_power, (void *)pycomplex5_pow},
    {Py_nb_positive, pycomplex5_pos},
    {Py_nb_negative, pycomplex5_neg},
    {Py_nb_bool, pycomplex5_bool},
    {0, NULL}
};

static PyType_Spec pycomplex5_spec = {
    .name = "algpy.algebras.multicomplex.complex5",
    .basicsize = sizeof(PyComplex5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex5_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npycomplex5_getitem(void *data, void *arr)
{
    complex5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_COMPLEX5(v);
}

static int
npycomplex5_setitem(PyObject *item, void *data, void *arr)
{
    complex5 v;
    int r = complex5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex5 needs exactly 32 components");
            return -1;
        }
        for (int i = 0; i < 32; i++) {
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
pycomplex5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycomplex5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex5_copyswapn(dst, sizeof(complex5), src, sizeof(complex5), 1,
                         swap, arr);
}

static npy_bool
pycomplex5_nonzero(void *data, void *arr)
{
    complex5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex5_arrfuncs;
typedef struct { char c; complex5 v; } pycomplex5_align_test;

static PyArray_DescrProto pycomplex5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex5), offsetof(pycomplex5_align_test, v), 0, 0, 0,
    &pycomplex5_arrfuncs,
};

#define COMPLEX5_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycomplex5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex5 *dst = (complex5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define COMPLEX5_REAL_TYPES(X) \
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

COMPLEX5_REAL_TYPES(COMPLEX5_CAST_FROM)

static int
pycomplex5_register_casts(void)
{
#define REGISTER_COMPLEX5_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex5,                    \
                                     pycomplex5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX5_REAL_TYPES(REGISTER_COMPLEX5_CAST)
#undef REGISTER_COMPLEX5_CAST
    return 0;
}

typedef complex5 (*complex5_unary_op)(complex5);
typedef complex5 (*complex5_binary_op)(complex5, complex5);
typedef int (*complex5_predicate)(complex5, complex5);

static void
loop_complex5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex5_unary_op op = (complex5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex5_binary_op op = (complex5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex5 a, b, r;
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
loop_complex5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex5_predicate op = (complex5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex5;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex5_register_ufuncs(PyObject *numpy)
{
    return
        pycomplex5_register_loop(numpy, "add", loop_complex5_binary, 3, 0, complex5_add) ||
        pycomplex5_register_loop(numpy, "subtract", loop_complex5_binary, 3, 0, complex5_sub) ||
        pycomplex5_register_loop(numpy, "multiply", loop_complex5_binary, 3, 0, complex5_mul) ||
        pycomplex5_register_loop(numpy, "divide", loop_complex5_binary, 3, 0, complex5_div) ||
        pycomplex5_register_loop(numpy, "power", loop_complex5_binary, 3, 0, complex5_rpow) ||
        pycomplex5_register_loop(numpy, "positive", loop_complex5_unary, 2, 0, complex5_pos) ||
        pycomplex5_register_loop(numpy, "negative", loop_complex5_unary, 2, 0, complex5_neg) ||
        pycomplex5_register_loop(numpy, "conjugate", loop_complex5_unary, 2, 0, complex5_conjugate) ||
        pycomplex5_register_loop(numpy, "reciprocal", loop_complex5_unary, 2, 0, complex5_inv) ||
        pycomplex5_register_loop(numpy, "square", loop_complex5_unary, 2, 0, complex5_square) ||
        pycomplex5_register_loop(numpy, "exp", loop_complex5_unary, 2, 0, complex5_exp) ||
        pycomplex5_register_loop(numpy, "expm1", loop_complex5_unary, 2, 0, complex5_expm1) ||
        pycomplex5_register_loop(numpy, "log", loop_complex5_unary, 2, 0, complex5_log) ||
        pycomplex5_register_loop(numpy, "log1p", loop_complex5_unary, 2, 0, complex5_log1p) ||
        pycomplex5_register_loop(numpy, "sqrt", loop_complex5_unary, 2, 0, complex5_sqrt) ||
        pycomplex5_register_loop(numpy, "sin", loop_complex5_unary, 2, 0, complex5_sin) ||
        pycomplex5_register_loop(numpy, "cos", loop_complex5_unary, 2, 0, complex5_cos) ||
        pycomplex5_register_loop(numpy, "tan", loop_complex5_unary, 2, 0, complex5_tan) ||
        pycomplex5_register_loop(numpy, "arcsin", loop_complex5_unary, 2, 0, complex5_asin) ||
        pycomplex5_register_loop(numpy, "arccos", loop_complex5_unary, 2, 0, complex5_acos) ||
        pycomplex5_register_loop(numpy, "arctan", loop_complex5_unary, 2, 0, complex5_atan) ||
        pycomplex5_register_loop(numpy, "sinh", loop_complex5_unary, 2, 0, complex5_sinh) ||
        pycomplex5_register_loop(numpy, "cosh", loop_complex5_unary, 2, 0, complex5_cosh) ||
        pycomplex5_register_loop(numpy, "tanh", loop_complex5_unary, 2, 0, complex5_tanh) ||
        pycomplex5_register_loop(numpy, "arcsinh", loop_complex5_unary, 2, 0, complex5_asinh) ||
        pycomplex5_register_loop(numpy, "arccosh", loop_complex5_unary, 2, 0, complex5_acosh) ||
        pycomplex5_register_loop(numpy, "arctanh", loop_complex5_unary, 2, 0, complex5_atanh) ||
        pycomplex5_register_loop(numpy, "equal", loop_complex5_predicate, 3, 1, complex5_equal) ||
        pycomplex5_register_loop(numpy, "not_equal", loop_complex5_predicate, 3, 1, complex5_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex5", "complex5 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex5(void)
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
    PyComplex5_Type = PyType_FromModuleAndSpec(
        m, &pycomplex5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplex5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex5_arrfuncs);
    pycomplex5_arrfuncs.getitem = npycomplex5_getitem;
    pycomplex5_arrfuncs.setitem = npycomplex5_setitem;
    pycomplex5_arrfuncs.copyswapn = pycomplex5_copyswapn;
    pycomplex5_arrfuncs.copyswap = pycomplex5_copyswap;
    pycomplex5_arrfuncs.nonzero = pycomplex5_nonzero;
    pycomplex5_descr_proto.typeobj = (PyTypeObject *)PyComplex5_Type;
    Py_SET_TYPE((PyObject *)&pycomplex5_descr_proto, &PyArrayDescr_Type);
    npy_complex5 = PyArray_RegisterDataType(&pycomplex5_descr_proto);
    if (npy_complex5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplex5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyComplex5_Type, 32) < 0) {
        goto fail;
    }
    if (pycomplex5_register_casts() < 0 ||
            pycomplex5_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex5", PyComplex5_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
