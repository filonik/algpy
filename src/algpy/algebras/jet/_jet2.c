/*
 * Jets of order 2, R[e]/(e^3), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in jet2.h.
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

#include "jet2.h"

typedef struct {
    PyObject_HEAD
    jet2 v;
} PyJet2;

static PyObject *PyJet2_Type = NULL;
static int npy_jet2 = -1;

static PyObject *
PyJet2_FromValue(PyTypeObject *type, jet2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyJet2 *p = (PyJet2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_JET2(v) \
    PyJet2_FromValue((PyTypeObject *)PyJet2_Type, (v))

static int
jet2_from_object(PyObject *o, jet2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyJet2_Type)) {
        *out = ((PyJet2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = jet2_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyjet2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x1", "x2", NULL};
    double d[3] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ddd", kwlist,
                                     &d[0], &d[1], &d[2])) {
        return NULL;
    }
    return PyJet2_FromValue(type, jet2_make(d));
}

static PyObject *
pyjet2_repr(PyObject *self)
{
    jet2 v = ((PyJet2 *)self)->v;
    char *s[3] = {NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 3; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("jet2(%s, %s, %s)",
                                  s[0], s[1], s[2]);
done:
    for (i = 0; i < 3; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyjet2_hash(PyObject *self)
{
    jet2 v = ((PyJet2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_JET2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (jet2_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(ddd)", v.data[0], v.data[1], v.data[2]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyjet2_richcompare(PyObject *a, PyObject *b, int op)
{
    jet2 x, y;
    int r;
    if ((r = jet2_from_object(a, &x)) <= 0 ||
            (r = jet2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(jet2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!jet2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define JET2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyjet2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        jet2 x, y;                                                       \
        int r;                                                               \
        if ((r = jet2_from_object(a, &x)) <= 0 ||                        \
                (r = jet2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_JET2(expr(x, y));                                     \
    }

JET2_BINARY(add, jet2_add)
JET2_BINARY(sub, jet2_sub)
JET2_BINARY(mul, jet2_mul)
JET2_BINARY(div, jet2_div)

static PyObject *
pyjet2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    jet2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "jet2 does not support modular power");
        return NULL;
    }
    if ((r = jet2_from_object(a, &x)) <= 0 ||
            (r = jet2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_JET2(jet2_rpow(x, y));
}

#define JET2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyjet2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_JET2(expr(((PyJet2 *)self)->v));                  \
    }

JET2_UNARY(pos, jet2_pos)
JET2_UNARY(neg, jet2_neg)

static PyObject *
pyjet2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue((PyTypeObject *)cls, jet2_zero());
}

static PyObject *
pyjet2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue((PyTypeObject *)cls, jet2_one());
}

static int
pyjet2_bool(PyObject *self)
{
    jet2 v = ((PyJet2 *)self)->v;
    return vec3_any_nonzero(v);
}

#define JET2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyjet2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_JET2(op(((PyJet2 *)self)->v));                   \
    }

JET2_METHOD(conjugate, jet2_conjugate)
JET2_METHOD(exp, jet2_exp)
JET2_METHOD(log, jet2_log)
JET2_METHOD(sqrt, jet2_sqrt)

static PyObject *
pyjet2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet2_norm(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_JET2(jet2_inv(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {3, 3};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet2_lmat(((PyJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet2_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {3, 3};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet2_rmat(((PyJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet2_minpoly(((PyJet2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[3];
    npy_intp dims[1] = {3};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    jet2_eivals(((PyJet2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 3; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyjet2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet2);
    jet2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet2 *)PyArray_DATA(arr);
    jet2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet2_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet2);
    jet2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet2 *)PyArray_DATA(arr);
    jet2_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {3, 3};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_jet2);
    jet2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (jet2 *)PyArray_DATA(arr);
    jet2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyjet2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    jet2 v = ((PyJet2 *)self)->v;
    return Py_BuildValue("(O(ddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2]);
}

static PyObject *
pyjet2_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_square(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_expm1(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_log1p(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_sin(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_cos(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_tan(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_asin(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_acos(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_atan(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_sinh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_cosh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_tanh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_asinh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_acosh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_atanh(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_re(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet2_FromValue(Py_TYPE(self), jet2_im(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_real(PyObject *self, void *closure)
{
    return pyjet2_re(self, NULL);
}

static PyObject *
pyjet2_imag(PyObject *self, void *closure)
{
    return pyjet2_im(self, NULL);
}

static PyObject *pyjet2_grades(PyObject *self, void *closure);

static PyGetSetDef pyjet2_getset[] = {
    {"grade", pyjet2_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pyjet2_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyjet2_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyjet2_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(jet2_is_real(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet2_to_re(((PyJet2 *)self)->v));
}

static PyObject *
pyjet2_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyJet2_FromValue((PyTypeObject *)cls, jet2_from_re(a));
}

static PyObject *
pyjet2_lpow(PyObject *self, PyObject *other)
{
    jet2 b;
    int r;
    r = jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet2_FromValue(Py_TYPE(self),
        jet2_lpow(((PyJet2 *)self)->v, b));
}

static PyObject *
pyjet2_rpow(PyObject *self, PyObject *other)
{
    jet2 b;
    int r;
    r = jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet2_FromValue(Py_TYPE(self),
        jet2_rpow(((PyJet2 *)self)->v, b));
}

static PyObject *
pyjet2_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    jet2 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = jet2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet2_FromValue(Py_TYPE(self),
        jet2_mix(((PyJet2 *)self)->v, b, t));
}

static PyObject *
pyjet2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_jet2, 3);
}

static PyObject *
pyjet2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_jet2, 3);
}

static PyObject *
pyjet2_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_JET2(jet2_ipow(((PyJet2 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyjet2_grade(PyObject *self, Py_ssize_t k)
{
    jet2 v = ((PyJet2 *)self)->v;
    switch (k) {
        case 0:
            return PyJet2_FromValue(Py_TYPE(self), jet2_grade0(v));
        case 1:
            return PyJet2_FromValue(Py_TYPE(self), jet2_grade1(v));
        case 2:
            return PyJet2_FromValue(Py_TYPE(self), jet2_grade2(v));
        default:
            return PyJet2_FromValue(Py_TYPE(self), jet2_zero());
    }
}

static PyObject *
pyjet2_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyjet2_grade, 3);
}

static PyMethodDef pyjet2_methods[] = {
    {"to_float", (PyCFunction)pyjet2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyjet2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyjet2_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyjet2_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyjet2_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyjet2_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyjet2_square, METH_NOARGS, "square."},
    {"expm1", pyjet2_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyjet2_log1p, METH_NOARGS, "log1p."},
    {"sin", pyjet2_sin, METH_NOARGS, "sin."},
    {"cos", pyjet2_cos, METH_NOARGS, "cos."},
    {"tan", pyjet2_tan, METH_NOARGS, "tan."},
    {"asin", pyjet2_asin, METH_NOARGS, "asin."},
    {"acos", pyjet2_acos, METH_NOARGS, "acos."},
    {"atan", pyjet2_atan, METH_NOARGS, "atan."},
    {"sinh", pyjet2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyjet2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyjet2_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyjet2_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyjet2_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyjet2_atanh, METH_NOARGS, "atanh."},
    {"re", pyjet2_re, METH_NOARGS, "re."},
    {"im", pyjet2_im, METH_NOARGS, "im."},
    {"arcsin", pyjet2_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyjet2_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyjet2_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyjet2_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyjet2_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyjet2_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyjet2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyjet2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyjet2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pyjet2_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pyjet2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyjet2_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyjet2_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyjet2_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyjet2_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyjet2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyjet2_conjugate, METH_NOARGS, "Conjugate e -> -e."},
    {"lmat", pyjet2_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyjet2_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyjet2_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyjet2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyjet2_exp, METH_NOARGS, "exp."},
    {"log", pyjet2_log, METH_NOARGS, "log."},
    {"sqrt", pyjet2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyjet2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyjet2_length(PyObject *Py_UNUSED(self))
{
    return 3;
}

static PyObject *
pyjet2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 3) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyJet2 *)self)->v.data[i]);
}

static PyObject *
pyjet2_subscript(PyObject *self, PyObject *key)
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
        i += 3;
    }
    return pyjet2_item(self, i);
}

static int
pyjet2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyjet2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyjet2_slots[] = {
    {Py_tp_dealloc, pyjet2_dealloc},
    {Py_tp_traverse, pyjet2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyjet2_repr},
    {Py_tp_str, pyjet2_repr},
    {Py_tp_hash, pyjet2_hash},
    {Py_tp_doc, (void *)"Jet of order 2 in R[e]/(e^3) with 3 float64 components"},
    {Py_tp_richcompare, pyjet2_richcompare},
    {Py_tp_methods, pyjet2_methods},
    {Py_tp_getset, pyjet2_getset},
    {Py_sq_length, pyjet2_length},
    {Py_sq_item, pyjet2_item},
    {Py_mp_subscript, pyjet2_subscript},
    {Py_tp_new, pyjet2_new},
    {Py_nb_add, pyjet2_add},
    {Py_nb_subtract, pyjet2_sub},
    {Py_nb_multiply, pyjet2_mul},
    {Py_nb_true_divide, pyjet2_div},
    {Py_nb_power, (void *)pyjet2_pow},
    {Py_nb_positive, pyjet2_pos},
    {Py_nb_negative, pyjet2_neg},
    {Py_nb_bool, pyjet2_bool},
    {0, NULL}
};

static PyType_Spec pyjet2_spec = {
    .name = "algpy.algebras.jet.jet2",
    .basicsize = sizeof(PyJet2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyjet2_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyjet2_getitem(void *data, void *arr)
{
    jet2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_JET2(v);
}

static int
npyjet2_setitem(PyObject *item, void *data, void *arr)
{
    jet2 v;
    int r = jet2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to jet2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 3) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "jet2 needs exactly 3 components");
            return -1;
        }
        for (int i = 0; i < 3; i++) {
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
pyjet2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(jet2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyjet2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyjet2_copyswapn(dst, sizeof(jet2), src, sizeof(jet2), 1,
                         swap, arr);
}

static npy_bool
pyjet2_nonzero(void *data, void *arr)
{
    jet2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec3_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyjet2_arrfuncs;
typedef struct { char c; jet2 v; } pyjet2_align_test;

static PyArray_DescrProto pyjet2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(jet2), offsetof(pyjet2_align_test, v), 0, 0, 0,
    &pyjet2_arrfuncs,
};

#define JET2_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyjet2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        jet2 *dst = (jet2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = jet2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define JET2_REAL_TYPES(X) \
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

JET2_REAL_TYPES(JET2_CAST_FROM)

static int
pyjet2_register_casts(void)
{
#define REGISTER_JET2_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_jet2,                    \
                                     pyjet2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_jet2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    JET2_REAL_TYPES(REGISTER_JET2_CAST)
#undef REGISTER_JET2_CAST
    return 0;
}

typedef jet2 (*jet2_unary_op)(jet2);
typedef jet2 (*jet2_binary_op)(jet2, jet2);
typedef int (*jet2_predicate)(jet2, jet2);

static void
loop_jet2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    jet2_unary_op op = (jet2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_jet2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    jet2_binary_op op = (jet2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet2 a, b, r;
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
loop_jet2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    jet2_predicate op = (jet2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyjet2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_jet2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_jet2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyjet2_register_ufuncs(PyObject *numpy)
{
    return
        pyjet2_register_loop(numpy, "add", loop_jet2_binary, 3, 0, jet2_add) ||
        pyjet2_register_loop(numpy, "subtract", loop_jet2_binary, 3, 0, jet2_sub) ||
        pyjet2_register_loop(numpy, "multiply", loop_jet2_binary, 3, 0, jet2_mul) ||
        pyjet2_register_loop(numpy, "divide", loop_jet2_binary, 3, 0, jet2_div) ||
        pyjet2_register_loop(numpy, "power", loop_jet2_binary, 3, 0, jet2_rpow) ||
        pyjet2_register_loop(numpy, "positive", loop_jet2_unary, 2, 0, jet2_pos) ||
        pyjet2_register_loop(numpy, "negative", loop_jet2_unary, 2, 0, jet2_neg) ||
        pyjet2_register_loop(numpy, "conjugate", loop_jet2_unary, 2, 0, jet2_conjugate) ||
        pyjet2_register_loop(numpy, "reciprocal", loop_jet2_unary, 2, 0, jet2_inv) ||
        pyjet2_register_loop(numpy, "square", loop_jet2_unary, 2, 0, jet2_square) ||
        pyjet2_register_loop(numpy, "exp", loop_jet2_unary, 2, 0, jet2_exp) ||
        pyjet2_register_loop(numpy, "expm1", loop_jet2_unary, 2, 0, jet2_expm1) ||
        pyjet2_register_loop(numpy, "log", loop_jet2_unary, 2, 0, jet2_log) ||
        pyjet2_register_loop(numpy, "log1p", loop_jet2_unary, 2, 0, jet2_log1p) ||
        pyjet2_register_loop(numpy, "sqrt", loop_jet2_unary, 2, 0, jet2_sqrt) ||
        pyjet2_register_loop(numpy, "sin", loop_jet2_unary, 2, 0, jet2_sin) ||
        pyjet2_register_loop(numpy, "cos", loop_jet2_unary, 2, 0, jet2_cos) ||
        pyjet2_register_loop(numpy, "tan", loop_jet2_unary, 2, 0, jet2_tan) ||
        pyjet2_register_loop(numpy, "arcsin", loop_jet2_unary, 2, 0, jet2_asin) ||
        pyjet2_register_loop(numpy, "arccos", loop_jet2_unary, 2, 0, jet2_acos) ||
        pyjet2_register_loop(numpy, "arctan", loop_jet2_unary, 2, 0, jet2_atan) ||
        pyjet2_register_loop(numpy, "sinh", loop_jet2_unary, 2, 0, jet2_sinh) ||
        pyjet2_register_loop(numpy, "cosh", loop_jet2_unary, 2, 0, jet2_cosh) ||
        pyjet2_register_loop(numpy, "tanh", loop_jet2_unary, 2, 0, jet2_tanh) ||
        pyjet2_register_loop(numpy, "arcsinh", loop_jet2_unary, 2, 0, jet2_asinh) ||
        pyjet2_register_loop(numpy, "arccosh", loop_jet2_unary, 2, 0, jet2_acosh) ||
        pyjet2_register_loop(numpy, "arctanh", loop_jet2_unary, 2, 0, jet2_atanh) ||
        pyjet2_register_loop(numpy, "equal", loop_jet2_predicate, 3, 1, jet2_equal) ||
        pyjet2_register_loop(numpy, "not_equal", loop_jet2_predicate, 3, 1, jet2_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_jet2", "Jet type of order 2 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__jet2(void)
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
    PyJet2_Type = PyType_FromModuleAndSpec(
        m, &pyjet2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyJet2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyjet2_arrfuncs);
    pyjet2_arrfuncs.getitem = npyjet2_getitem;
    pyjet2_arrfuncs.setitem = npyjet2_setitem;
    pyjet2_arrfuncs.copyswapn = pyjet2_copyswapn;
    pyjet2_arrfuncs.copyswap = pyjet2_copyswap;
    pyjet2_arrfuncs.nonzero = pyjet2_nonzero;
    pyjet2_descr_proto.typeobj = (PyTypeObject *)PyJet2_Type;
    Py_SET_TYPE((PyObject *)&pyjet2_descr_proto, &PyArrayDescr_Type);
    npy_jet2 = PyArray_RegisterDataType(&pyjet2_descr_proto);
    if (npy_jet2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_jet2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyJet2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyJet2_Type, 3) < 0) {
        goto fail;
    }
    if (pyjet2_register_casts() < 0 ||
            pyjet2_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "jet2", PyJet2_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
