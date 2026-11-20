/*
 * Jets of order 5, R[e]/(e^6), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in jet5.h.
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

#include "jet5.h"

typedef struct {
    PyObject_HEAD
    jet5 v;
} PyJet5;

static PyObject *PyJet5_Type = NULL;
static int npy_jet5 = -1;

static PyObject *
PyJet5_FromValue(PyTypeObject *type, jet5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyJet5 *p = (PyJet5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_JET5(v) \
    PyJet5_FromValue((PyTypeObject *)PyJet5_Type, (v))

static int
jet5_from_object(PyObject *o, jet5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyJet5_Type)) {
        *out = ((PyJet5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = jet5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyjet5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x1", "x2", "x3", "x4", "x5", NULL};
    double d[6] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddd", kwlist,
                                     &d[0], &d[1], &d[2], &d[3], &d[4], &d[5])) {
        return NULL;
    }
    return PyJet5_FromValue(type, jet5_make(d));
}

static PyObject *
pyjet5_repr(PyObject *self)
{
    jet5 v = ((PyJet5 *)self)->v;
    char *s[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 6; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("jet5(%s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5]);
done:
    for (i = 0; i < 6; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyjet5_hash(PyObject *self)
{
    jet5 v = ((PyJet5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_JET5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (jet5_is_real(v)) {
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
pyjet5_richcompare(PyObject *a, PyObject *b, int op)
{
    jet5 x, y;
    int r;
    if ((r = jet5_from_object(a, &x)) <= 0 ||
            (r = jet5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(jet5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!jet5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define JET5_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyjet5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        jet5 x, y;                                                       \
        int r;                                                               \
        if ((r = jet5_from_object(a, &x)) <= 0 ||                        \
                (r = jet5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_JET5(expr(x, y));                                     \
    }

JET5_BINARY(add, jet5_add)
JET5_BINARY(sub, jet5_sub)
JET5_BINARY(mul, jet5_mul)
JET5_BINARY(div, jet5_div)

static PyObject *
pyjet5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    jet5 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "jet5 does not support modular power");
        return NULL;
    }
    if ((r = jet5_from_object(a, &x)) <= 0 ||
            (r = jet5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_JET5(jet5_rpow(x, y));
}

#define JET5_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyjet5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_JET5(expr(((PyJet5 *)self)->v));                  \
    }

JET5_UNARY(pos, jet5_pos)
JET5_UNARY(neg, jet5_neg)

static PyObject *
pyjet5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue((PyTypeObject *)cls, jet5_zero());
}

static PyObject *
pyjet5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue((PyTypeObject *)cls, jet5_one());
}

static int
pyjet5_bool(PyObject *self)
{
    jet5 v = ((PyJet5 *)self)->v;
    return vec6_any_nonzero(v);
}

#define JET5_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyjet5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_JET5(op(((PyJet5 *)self)->v));                   \
    }

JET5_METHOD(conjugate, jet5_conjugate)
JET5_METHOD(exp, jet5_exp)
JET5_METHOD(log, jet5_log)
JET5_METHOD(sqrt, jet5_sqrt)

static PyObject *
pyjet5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet5_norm(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_JET5(jet5_inv(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet5_lmat(((PyJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet5_rmat(((PyJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {7};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet5_minpoly(((PyJet5 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[6];
    npy_intp dims[1] = {6};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    jet5_eivals(((PyJet5 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 6; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyjet5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {6};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet5);
    jet5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet5 *)PyArray_DATA(arr);
    jet5_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet5);
    jet5 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet5 *)PyArray_DATA(arr);
    jet5_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet5_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {6, 6};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_jet5);
    jet5 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (jet5 *)PyArray_DATA(arr);
    jet5_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyjet5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    jet5 v = ((PyJet5 *)self)->v;
    return Py_BuildValue("(O(dddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5]);
}

static PyObject *
pyjet5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_square(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_expm1(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_log1p(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_sin(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_cos(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_tan(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_asin(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_acos(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_atan(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_sinh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_cosh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_tanh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_asinh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_acosh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_atanh(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_re(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet5_FromValue(Py_TYPE(self), jet5_im(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_real(PyObject *self, void *closure)
{
    return pyjet5_re(self, NULL);
}

static PyObject *
pyjet5_imag(PyObject *self, void *closure)
{
    return pyjet5_im(self, NULL);
}

static PyObject *pyjet5_grades(PyObject *self, void *closure);

static PyGetSetDef pyjet5_getset[] = {
    {"grade", pyjet5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pyjet5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyjet5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyjet5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(jet5_is_real(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet5_to_re(((PyJet5 *)self)->v));
}

static PyObject *
pyjet5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyJet5_FromValue((PyTypeObject *)cls, jet5_from_re(a));
}

static PyObject *
pyjet5_lpow(PyObject *self, PyObject *other)
{
    jet5 b;
    int r;
    r = jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet5_FromValue(Py_TYPE(self),
        jet5_lpow(((PyJet5 *)self)->v, b));
}

static PyObject *
pyjet5_rpow(PyObject *self, PyObject *other)
{
    jet5 b;
    int r;
    r = jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet5_FromValue(Py_TYPE(self),
        jet5_rpow(((PyJet5 *)self)->v, b));
}

static PyObject *
pyjet5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    jet5 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = jet5_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet5_FromValue(Py_TYPE(self),
        jet5_mix(((PyJet5 *)self)->v, b, t));
}

static PyObject *
pyjet5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_jet5, 6);
}

static PyObject *
pyjet5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_jet5, 6);
}

static PyObject *
pyjet5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_JET5(jet5_ipow(((PyJet5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyjet5_grade(PyObject *self, Py_ssize_t k)
{
    jet5 v = ((PyJet5 *)self)->v;
    switch (k) {
        case 0:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade0(v));
        case 1:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade1(v));
        case 2:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade2(v));
        case 3:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade3(v));
        case 4:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade4(v));
        case 5:
            return PyJet5_FromValue(Py_TYPE(self), jet5_grade5(v));
        default:
            return PyJet5_FromValue(Py_TYPE(self), jet5_zero());
    }
}

static PyObject *
pyjet5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyjet5_grade, 6);
}

static PyMethodDef pyjet5_methods[] = {
    {"to_float", (PyCFunction)pyjet5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyjet5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyjet5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyjet5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyjet5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyjet5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyjet5_square, METH_NOARGS, "square."},
    {"expm1", pyjet5_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyjet5_log1p, METH_NOARGS, "log1p."},
    {"sin", pyjet5_sin, METH_NOARGS, "sin."},
    {"cos", pyjet5_cos, METH_NOARGS, "cos."},
    {"tan", pyjet5_tan, METH_NOARGS, "tan."},
    {"asin", pyjet5_asin, METH_NOARGS, "asin."},
    {"acos", pyjet5_acos, METH_NOARGS, "acos."},
    {"atan", pyjet5_atan, METH_NOARGS, "atan."},
    {"sinh", pyjet5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyjet5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyjet5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyjet5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyjet5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyjet5_atanh, METH_NOARGS, "atanh."},
    {"re", pyjet5_re, METH_NOARGS, "re."},
    {"im", pyjet5_im, METH_NOARGS, "im."},
    {"arcsin", pyjet5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyjet5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyjet5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyjet5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyjet5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyjet5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyjet5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyjet5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyjet5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pyjet5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pyjet5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyjet5_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyjet5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyjet5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyjet5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyjet5_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyjet5_conjugate, METH_NOARGS, "Conjugate e -> -e."},
    {"lmat", pyjet5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyjet5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyjet5_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyjet5_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyjet5_exp, METH_NOARGS, "exp."},
    {"log", pyjet5_log, METH_NOARGS, "log."},
    {"sqrt", pyjet5_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyjet5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyjet5_length(PyObject *Py_UNUSED(self))
{
    return 6;
}

static PyObject *
pyjet5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 6) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyJet5 *)self)->v.data[i]);
}

static PyObject *
pyjet5_subscript(PyObject *self, PyObject *key)
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
    return pyjet5_item(self, i);
}

static int
pyjet5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyjet5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyjet5_slots[] = {
    {Py_tp_dealloc, pyjet5_dealloc},
    {Py_tp_traverse, pyjet5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyjet5_repr},
    {Py_tp_str, pyjet5_repr},
    {Py_tp_hash, pyjet5_hash},
    {Py_tp_doc, (void *)"Jet of order 5 in R[e]/(e^6) with 6 float64 components"},
    {Py_tp_richcompare, pyjet5_richcompare},
    {Py_tp_methods, pyjet5_methods},
    {Py_tp_getset, pyjet5_getset},
    {Py_sq_length, pyjet5_length},
    {Py_sq_item, pyjet5_item},
    {Py_mp_subscript, pyjet5_subscript},
    {Py_tp_new, pyjet5_new},
    {Py_nb_add, pyjet5_add},
    {Py_nb_subtract, pyjet5_sub},
    {Py_nb_multiply, pyjet5_mul},
    {Py_nb_true_divide, pyjet5_div},
    {Py_nb_power, (void *)pyjet5_pow},
    {Py_nb_positive, pyjet5_pos},
    {Py_nb_negative, pyjet5_neg},
    {Py_nb_bool, pyjet5_bool},
    {0, NULL}
};

static PyType_Spec pyjet5_spec = {
    .name = "algpy.algebras.jet.jet5",
    .basicsize = sizeof(PyJet5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyjet5_slots,
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
npyjet5_getitem(void *data, void *arr)
{
    jet5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_JET5(v);
}

static int
npyjet5_setitem(PyObject *item, void *data, void *arr)
{
    jet5 v;
    int r = jet5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to jet5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 6) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "jet5 needs exactly 6 components");
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
pyjet5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(jet5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyjet5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyjet5_copyswapn(dst, sizeof(jet5), src, sizeof(jet5), 1,
                         swap, arr);
}

static npy_bool
pyjet5_nonzero(void *data, void *arr)
{
    jet5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec6_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyjet5_arrfuncs;
typedef struct { char c; jet5 v; } pyjet5_align_test;

static PyArray_DescrProto pyjet5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(jet5), offsetof(pyjet5_align_test, v), 0, 0, 0,
    &pyjet5_arrfuncs,
};

#define JET5_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyjet5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        jet5 *dst = (jet5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = jet5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define JET5_REAL_TYPES(X) \
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

JET5_REAL_TYPES(JET5_CAST_FROM)

static int
pyjet5_register_casts(void)
{
#define REGISTER_JET5_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_jet5,                    \
                                     pyjet5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_jet5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    JET5_REAL_TYPES(REGISTER_JET5_CAST)
#undef REGISTER_JET5_CAST
    return 0;
}

typedef jet5 (*jet5_unary_op)(jet5);
typedef jet5 (*jet5_binary_op)(jet5, jet5);
typedef int (*jet5_predicate)(jet5, jet5);

static void
loop_jet5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    jet5_unary_op op = (jet5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_jet5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    jet5_binary_op op = (jet5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet5 a, b, r;
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
loop_jet5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    jet5_predicate op = (jet5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyjet5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_jet5;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_jet5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyjet5_register_ufuncs(PyObject *numpy)
{
    return
        pyjet5_register_loop(numpy, "add", loop_jet5_binary, 3, 0, jet5_add) ||
        pyjet5_register_loop(numpy, "subtract", loop_jet5_binary, 3, 0, jet5_sub) ||
        pyjet5_register_loop(numpy, "multiply", loop_jet5_binary, 3, 0, jet5_mul) ||
        pyjet5_register_loop(numpy, "divide", loop_jet5_binary, 3, 0, jet5_div) ||
        pyjet5_register_loop(numpy, "power", loop_jet5_binary, 3, 0, jet5_rpow) ||
        pyjet5_register_loop(numpy, "positive", loop_jet5_unary, 2, 0, jet5_pos) ||
        pyjet5_register_loop(numpy, "negative", loop_jet5_unary, 2, 0, jet5_neg) ||
        pyjet5_register_loop(numpy, "conjugate", loop_jet5_unary, 2, 0, jet5_conjugate) ||
        pyjet5_register_loop(numpy, "reciprocal", loop_jet5_unary, 2, 0, jet5_inv) ||
        pyjet5_register_loop(numpy, "square", loop_jet5_unary, 2, 0, jet5_square) ||
        pyjet5_register_loop(numpy, "exp", loop_jet5_unary, 2, 0, jet5_exp) ||
        pyjet5_register_loop(numpy, "expm1", loop_jet5_unary, 2, 0, jet5_expm1) ||
        pyjet5_register_loop(numpy, "log", loop_jet5_unary, 2, 0, jet5_log) ||
        pyjet5_register_loop(numpy, "log1p", loop_jet5_unary, 2, 0, jet5_log1p) ||
        pyjet5_register_loop(numpy, "sqrt", loop_jet5_unary, 2, 0, jet5_sqrt) ||
        pyjet5_register_loop(numpy, "sin", loop_jet5_unary, 2, 0, jet5_sin) ||
        pyjet5_register_loop(numpy, "cos", loop_jet5_unary, 2, 0, jet5_cos) ||
        pyjet5_register_loop(numpy, "tan", loop_jet5_unary, 2, 0, jet5_tan) ||
        pyjet5_register_loop(numpy, "arcsin", loop_jet5_unary, 2, 0, jet5_asin) ||
        pyjet5_register_loop(numpy, "arccos", loop_jet5_unary, 2, 0, jet5_acos) ||
        pyjet5_register_loop(numpy, "arctan", loop_jet5_unary, 2, 0, jet5_atan) ||
        pyjet5_register_loop(numpy, "sinh", loop_jet5_unary, 2, 0, jet5_sinh) ||
        pyjet5_register_loop(numpy, "cosh", loop_jet5_unary, 2, 0, jet5_cosh) ||
        pyjet5_register_loop(numpy, "tanh", loop_jet5_unary, 2, 0, jet5_tanh) ||
        pyjet5_register_loop(numpy, "arcsinh", loop_jet5_unary, 2, 0, jet5_asinh) ||
        pyjet5_register_loop(numpy, "arccosh", loop_jet5_unary, 2, 0, jet5_acosh) ||
        pyjet5_register_loop(numpy, "arctanh", loop_jet5_unary, 2, 0, jet5_atanh) ||
        pyjet5_register_loop(numpy, "equal", loop_jet5_predicate, 3, 1, jet5_equal) ||
        pyjet5_register_loop(numpy, "not_equal", loop_jet5_predicate, 3, 1, jet5_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_jet5", "Jet type of order 5 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__jet5(void)
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
    PyJet5_Type = PyType_FromModuleAndSpec(
        m, &pyjet5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyJet5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyjet5_arrfuncs);
    pyjet5_arrfuncs.getitem = npyjet5_getitem;
    pyjet5_arrfuncs.setitem = npyjet5_setitem;
    pyjet5_arrfuncs.copyswapn = pyjet5_copyswapn;
    pyjet5_arrfuncs.copyswap = pyjet5_copyswap;
    pyjet5_arrfuncs.nonzero = pyjet5_nonzero;
    pyjet5_descr_proto.typeobj = (PyTypeObject *)PyJet5_Type;
    Py_SET_TYPE((PyObject *)&pyjet5_descr_proto, &PyArrayDescr_Type);
    npy_jet5 = PyArray_RegisterDataType(&pyjet5_descr_proto);
    if (npy_jet5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_jet5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyJet5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyJet5_Type, 6) < 0) {
        goto fail;
    }
    if (pyjet5_register_casts() < 0 ||
            pyjet5_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "jet5", PyJet5_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
