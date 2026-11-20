/*
 * dual5 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in dual5.h.
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

#include "dual5.h"

typedef struct {
    PyObject_HEAD
    dual5 v;
} PyDual5;

static PyObject *PyDual5_Type = NULL;
static int npy_dual5 = -1;

static PyObject *
PyDual5_FromValue(PyTypeObject *type, dual5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyDual5 *p = (PyDual5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_DUAL5(v) \
    PyDual5_FromValue((PyTypeObject *)PyDual5_Type, (v))

static int
dual5_from_object(PyObject *o, dual5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyDual5_Type)) {
        *out = ((PyDual5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = dual5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pydual5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyDual5_FromValue(type, dual5_make(d));
}

static PyObject *
pydual5_repr(PyObject *self)
{
    dual5 v = ((PyDual5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("dual5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pydual5_hash(PyObject *self)
{
    dual5 v = ((PyDual5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_DUAL5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (dual5_is_real(v)) {
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
pydual5_richcompare(PyObject *a, PyObject *b, int op)
{
    dual5 x, y;
    int r;
    if ((r = dual5_from_object(a, &x)) <= 0 ||
            (r = dual5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(dual5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!dual5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define DUAL5_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pydual5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        dual5 x, y;                                                       \
        int r;                                                               \
        if ((r = dual5_from_object(a, &x)) <= 0 ||                        \
                (r = dual5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_DUAL5(expr(x, y));                                     \
    }

DUAL5_BINARY(add, dual5_add)
DUAL5_BINARY(sub, dual5_sub)
DUAL5_BINARY(mul, dual5_mul)
DUAL5_BINARY(div, dual5_div)

static PyObject *
pydual5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    dual5 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "dual5 does not support modular power");
        return NULL;
    }
    if ((r = dual5_from_object(a, &x)) <= 0 ||
            (r = dual5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_DUAL5(dual5_rpow(x, y));
}

#define DUAL5_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pydual5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_DUAL5(expr(((PyDual5 *)self)->v));                  \
    }

DUAL5_UNARY(pos, dual5_pos)
DUAL5_UNARY(neg, dual5_neg)

static PyObject *
pydual5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue((PyTypeObject *)cls, dual5_zero());
}

static PyObject *
pydual5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue((PyTypeObject *)cls, dual5_one());
}

static int
pydual5_bool(PyObject *self)
{
    dual5 v = ((PyDual5 *)self)->v;
    return vec32_any_nonzero(v);
}

#define DUAL5_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pydual5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_DUAL5(op(((PyDual5 *)self)->v));                   \
    }

DUAL5_METHOD(conjugate, dual5_conjugate)
DUAL5_METHOD(exp, dual5_exp)
DUAL5_METHOD(log, dual5_log)
DUAL5_METHOD(sqrt, dual5_sqrt)

static PyObject *
pydual5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual5_norm(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_DUAL5(dual5_inv(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual5_lmat(((PyDual5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual5_rmat(((PyDual5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {33};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    dual5_minpoly(((PyDual5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pydual5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[32];
    npy_intp dims[1] = {32};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    dual5_eivals(((PyDual5 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 32; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pydual5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual5);
    dual5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual5 *)PyArray_DATA(arr);
    dual5_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_dual5);
    dual5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (dual5 *)PyArray_DATA(arr);
    dual5_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pydual5_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_dual5);
    dual5 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (dual5 *)PyArray_DATA(arr);
    dual5_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pydual5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual5 v = ((PyDual5 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pydual5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_square(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_expm1(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_log1p(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_sin(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_cos(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_tan(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_asin(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_acos(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_atan(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_sinh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_cosh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_tanh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_asinh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_acosh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_atanh(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_re(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyDual5_FromValue(Py_TYPE(self), dual5_im(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_real(PyObject *self, void *closure)
{
    return pydual5_re(self, NULL);
}

static PyObject *
pydual5_imag(PyObject *self, void *closure)
{
    return pydual5_im(self, NULL);
}

static PyObject *pydual5_grades(PyObject *self, void *closure);

static PyGetSetDef pydual5_getset[] = {
    {"grade", pydual5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pydual5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pydual5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pydual5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(dual5_is_real(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(dual5_to_re(((PyDual5 *)self)->v));
}

static PyObject *
pydual5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyDual5_FromValue((PyTypeObject *)cls, dual5_from_re(a));
}

static PyObject *
pydual5_lpow(PyObject *self, PyObject *other)
{
    dual5 b;
    int r;
    r = dual5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual5_FromValue(Py_TYPE(self),
        dual5_lpow(((PyDual5 *)self)->v, b));
}

static PyObject *
pydual5_rpow(PyObject *self, PyObject *other)
{
    dual5 b;
    int r;
    r = dual5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual5_FromValue(Py_TYPE(self),
        dual5_rpow(((PyDual5 *)self)->v, b));
}

static PyObject *
pydual5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    dual5 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = dual5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyDual5_FromValue(Py_TYPE(self),
        dual5_mix(((PyDual5 *)self)->v, b, t));
}

static PyObject *
pydual5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_dual5, 32);
}

static PyObject *
pydual5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_dual5, 32);
}

static PyObject *
pydual5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_DUAL5(dual5_ipow(((PyDual5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pydual5_grade(PyObject *self, Py_ssize_t k)
{
    dual5 v = ((PyDual5 *)self)->v;
    switch (k) {
        case 0:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade0(v));
        case 1:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade1(v));
        case 2:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade2(v));
        case 3:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade3(v));
        case 4:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade4(v));
        case 5:
            return PyDual5_FromValue(Py_TYPE(self), dual5_grade5(v));
        default:
            return PyDual5_FromValue(Py_TYPE(self), dual5_zero());
    }
}

static PyObject *
pydual5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pydual5_grade, 6);
}

static PyMethodDef pydual5_methods[] = {
    {"to_float", (PyCFunction)pydual5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pydual5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pydual5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pydual5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pydual5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pydual5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pydual5_square, METH_NOARGS, "square."},
    {"expm1", pydual5_expm1, METH_NOARGS, "expm1."},
    {"log1p", pydual5_log1p, METH_NOARGS, "log1p."},
    {"sin", pydual5_sin, METH_NOARGS, "sin."},
    {"cos", pydual5_cos, METH_NOARGS, "cos."},
    {"tan", pydual5_tan, METH_NOARGS, "tan."},
    {"asin", pydual5_asin, METH_NOARGS, "asin."},
    {"acos", pydual5_acos, METH_NOARGS, "acos."},
    {"atan", pydual5_atan, METH_NOARGS, "atan."},
    {"sinh", pydual5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pydual5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pydual5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pydual5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pydual5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pydual5_atanh, METH_NOARGS, "atanh."},
    {"re", pydual5_re, METH_NOARGS, "re."},
    {"im", pydual5_im, METH_NOARGS, "im."},
    {"arcsin", pydual5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pydual5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pydual5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pydual5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pydual5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pydual5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pydual5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pydual5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pydual5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pydual5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pydual5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pydual5_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pydual5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pydual5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pydual5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pydual5_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pydual5_conjugate, METH_NOARGS, "Conjugation of the outer dual generator."},
    {"lmat", pydual5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pydual5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pydual5_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pydual5_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pydual5_exp, METH_NOARGS, "exp."},
    {"log", pydual5_log, METH_NOARGS, "log."},
    {"sqrt", pydual5_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pydual5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pydual5_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pydual5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyDual5 *)self)->v.data[i]);
}

static PyObject *
pydual5_subscript(PyObject *self, PyObject *key)
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
    return pydual5_item(self, i);
}

static int
pydual5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pydual5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pydual5_slots[] = {
    {Py_tp_dealloc, pydual5_dealloc},
    {Py_tp_traverse, pydual5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pydual5_repr},
    {Py_tp_str, pydual5_repr},
    {Py_tp_hash, pydual5_hash},
    {Py_tp_doc, (void *)"Fivefold dual number with 32 float64 coefficients"},
    {Py_tp_richcompare, pydual5_richcompare},
    {Py_tp_methods, pydual5_methods},
    {Py_tp_getset, pydual5_getset},
    {Py_sq_length, pydual5_length},
    {Py_sq_item, pydual5_item},
    {Py_mp_subscript, pydual5_subscript},
    {Py_tp_new, pydual5_new},
    {Py_nb_add, pydual5_add},
    {Py_nb_subtract, pydual5_sub},
    {Py_nb_multiply, pydual5_mul},
    {Py_nb_true_divide, pydual5_div},
    {Py_nb_power, (void *)pydual5_pow},
    {Py_nb_positive, pydual5_pos},
    {Py_nb_negative, pydual5_neg},
    {Py_nb_bool, pydual5_bool},
    {0, NULL}
};

static PyType_Spec pydual5_spec = {
    .name = "algpy.algebras.multicomplex.dual5",
    .basicsize = sizeof(PyDual5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pydual5_slots,
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
npydual5_getitem(void *data, void *arr)
{
    dual5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_DUAL5(v);
}

static int
npydual5_setitem(PyObject *item, void *data, void *arr)
{
    dual5 v;
    int r = dual5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to dual5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "dual5 needs exactly 32 components");
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
pydual5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(dual5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pydual5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pydual5_copyswapn(dst, sizeof(dual5), src, sizeof(dual5), 1,
                         swap, arr);
}

static npy_bool
pydual5_nonzero(void *data, void *arr)
{
    dual5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pydual5_arrfuncs;
typedef struct { char c; dual5 v; } pydual5_align_test;

static PyArray_DescrProto pydual5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(dual5), offsetof(pydual5_align_test, v), 0, 0, 0,
    &pydual5_arrfuncs,
};

#define DUAL5_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pydual5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        dual5 *dst = (dual5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = dual5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define DUAL5_REAL_TYPES(X) \
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

DUAL5_REAL_TYPES(DUAL5_CAST_FROM)

static int
pydual5_register_casts(void)
{
#define REGISTER_DUAL5_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_dual5,                    \
                                     pydual5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_dual5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    DUAL5_REAL_TYPES(REGISTER_DUAL5_CAST)
#undef REGISTER_DUAL5_CAST
    return 0;
}

typedef dual5 (*dual5_unary_op)(dual5);
typedef dual5 (*dual5_binary_op)(dual5, dual5);
typedef int (*dual5_predicate)(dual5, dual5);

static void
loop_dual5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    dual5_unary_op op = (dual5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_dual5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    dual5_binary_op op = (dual5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual5 a, b, r;
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
loop_dual5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    dual5_predicate op = (dual5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        dual5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pydual5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_dual5;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_dual5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pydual5_register_ufuncs(PyObject *numpy)
{
    return
        pydual5_register_loop(numpy, "add", loop_dual5_binary, 3, 0, dual5_add) ||
        pydual5_register_loop(numpy, "subtract", loop_dual5_binary, 3, 0, dual5_sub) ||
        pydual5_register_loop(numpy, "multiply", loop_dual5_binary, 3, 0, dual5_mul) ||
        pydual5_register_loop(numpy, "divide", loop_dual5_binary, 3, 0, dual5_div) ||
        pydual5_register_loop(numpy, "power", loop_dual5_binary, 3, 0, dual5_rpow) ||
        pydual5_register_loop(numpy, "positive", loop_dual5_unary, 2, 0, dual5_pos) ||
        pydual5_register_loop(numpy, "negative", loop_dual5_unary, 2, 0, dual5_neg) ||
        pydual5_register_loop(numpy, "conjugate", loop_dual5_unary, 2, 0, dual5_conjugate) ||
        pydual5_register_loop(numpy, "reciprocal", loop_dual5_unary, 2, 0, dual5_inv) ||
        pydual5_register_loop(numpy, "square", loop_dual5_unary, 2, 0, dual5_square) ||
        pydual5_register_loop(numpy, "exp", loop_dual5_unary, 2, 0, dual5_exp) ||
        pydual5_register_loop(numpy, "expm1", loop_dual5_unary, 2, 0, dual5_expm1) ||
        pydual5_register_loop(numpy, "log", loop_dual5_unary, 2, 0, dual5_log) ||
        pydual5_register_loop(numpy, "log1p", loop_dual5_unary, 2, 0, dual5_log1p) ||
        pydual5_register_loop(numpy, "sqrt", loop_dual5_unary, 2, 0, dual5_sqrt) ||
        pydual5_register_loop(numpy, "sin", loop_dual5_unary, 2, 0, dual5_sin) ||
        pydual5_register_loop(numpy, "cos", loop_dual5_unary, 2, 0, dual5_cos) ||
        pydual5_register_loop(numpy, "tan", loop_dual5_unary, 2, 0, dual5_tan) ||
        pydual5_register_loop(numpy, "arcsin", loop_dual5_unary, 2, 0, dual5_asin) ||
        pydual5_register_loop(numpy, "arccos", loop_dual5_unary, 2, 0, dual5_acos) ||
        pydual5_register_loop(numpy, "arctan", loop_dual5_unary, 2, 0, dual5_atan) ||
        pydual5_register_loop(numpy, "sinh", loop_dual5_unary, 2, 0, dual5_sinh) ||
        pydual5_register_loop(numpy, "cosh", loop_dual5_unary, 2, 0, dual5_cosh) ||
        pydual5_register_loop(numpy, "tanh", loop_dual5_unary, 2, 0, dual5_tanh) ||
        pydual5_register_loop(numpy, "arcsinh", loop_dual5_unary, 2, 0, dual5_asinh) ||
        pydual5_register_loop(numpy, "arccosh", loop_dual5_unary, 2, 0, dual5_acosh) ||
        pydual5_register_loop(numpy, "arctanh", loop_dual5_unary, 2, 0, dual5_atanh) ||
        pydual5_register_loop(numpy, "equal", loop_dual5_predicate, 3, 1, dual5_equal) ||
        pydual5_register_loop(numpy, "not_equal", loop_dual5_predicate, 3, 1, dual5_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_dual5", "dual5 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__dual5(void)
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
    PyDual5_Type = PyType_FromModuleAndSpec(
        m, &pydual5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyDual5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pydual5_arrfuncs);
    pydual5_arrfuncs.getitem = npydual5_getitem;
    pydual5_arrfuncs.setitem = npydual5_setitem;
    pydual5_arrfuncs.copyswapn = pydual5_copyswapn;
    pydual5_arrfuncs.copyswap = pydual5_copyswap;
    pydual5_arrfuncs.nonzero = pydual5_nonzero;
    pydual5_descr_proto.typeobj = (PyTypeObject *)PyDual5_Type;
    Py_SET_TYPE((PyObject *)&pydual5_descr_proto, &PyArrayDescr_Type);
    npy_dual5 = PyArray_RegisterDataType(&pydual5_descr_proto);
    if (npy_dual5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_dual5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyDual5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyDual5_Type, 32) < 0) {
        goto fail;
    }
    if (pydual5_register_casts() < 0 ||
            pydual5_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "dual5", PyDual5_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
