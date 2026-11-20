/*
 * Jets of order 1, R[e]/(e^2), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in jet1.h.
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

#include "jet1.h"

typedef struct {
    PyObject_HEAD
    jet1 v;
} PyJet1;

static PyObject *PyJet1_Type = NULL;
static int npy_jet1 = -1;

static PyObject *
PyJet1_FromValue(PyTypeObject *type, jet1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyJet1 *p = (PyJet1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_JET1(v) \
    PyJet1_FromValue((PyTypeObject *)PyJet1_Type, (v))

static int
jet1_from_object(PyObject *o, jet1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyJet1_Type)) {
        *out = ((PyJet1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = jet1_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyjet1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x1", NULL};
    double d[2] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist,
                                     &d[0], &d[1])) {
        return NULL;
    }
    return PyJet1_FromValue(type, jet1_make(d));
}

static PyObject *
pyjet1_repr(PyObject *self)
{
    jet1 v = ((PyJet1 *)self)->v;
    char *s[2] = {NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 2; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("jet1(%s, %s)",
                                  s[0], s[1]);
done:
    for (i = 0; i < 2; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyjet1_hash(PyObject *self)
{
    jet1 v = ((PyJet1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_JET1_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (jet1_is_real(v)) {
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
pyjet1_richcompare(PyObject *a, PyObject *b, int op)
{
    jet1 x, y;
    int r;
    if ((r = jet1_from_object(a, &x)) <= 0 ||
            (r = jet1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(jet1_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!jet1_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define JET1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyjet1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        jet1 x, y;                                                       \
        int r;                                                               \
        if ((r = jet1_from_object(a, &x)) <= 0 ||                        \
                (r = jet1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_JET1(expr(x, y));                                     \
    }

JET1_BINARY(add, jet1_add)
JET1_BINARY(sub, jet1_sub)
JET1_BINARY(mul, jet1_mul)
JET1_BINARY(div, jet1_div)

static PyObject *
pyjet1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    jet1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "jet1 does not support modular power");
        return NULL;
    }
    if ((r = jet1_from_object(a, &x)) <= 0 ||
            (r = jet1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_JET1(jet1_rpow(x, y));
}

#define JET1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyjet1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_JET1(expr(((PyJet1 *)self)->v));                  \
    }

JET1_UNARY(pos, jet1_pos)
JET1_UNARY(neg, jet1_neg)

static PyObject *
pyjet1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue((PyTypeObject *)cls, jet1_zero());
}

static PyObject *
pyjet1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue((PyTypeObject *)cls, jet1_one());
}

static int
pyjet1_bool(PyObject *self)
{
    jet1 v = ((PyJet1 *)self)->v;
    return vec2_any_nonzero(v);
}

#define JET1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyjet1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_JET1(op(((PyJet1 *)self)->v));                   \
    }

JET1_METHOD(conjugate, jet1_conjugate)
JET1_METHOD(exp, jet1_exp)
JET1_METHOD(log, jet1_log)
JET1_METHOD(sqrt, jet1_sqrt)

static PyObject *
pyjet1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet1_norm(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_JET1(jet1_inv(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet1_lmat(((PyJet1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet1_rmat(((PyJet1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet1_minpoly(((PyJet1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[2];
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    jet1_eivals(((PyJet1 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyjet1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet1);
    jet1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet1 *)PyArray_DATA(arr);
    jet1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet1_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet1);
    jet1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet1 *)PyArray_DATA(arr);
    jet1_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_jet1);
    jet1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (jet1 *)PyArray_DATA(arr);
    jet1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyjet1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    jet1 v = ((PyJet1 *)self)->v;
    return Py_BuildValue("(O(dd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1]);
}

static PyObject *
pyjet1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_square(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_expm1(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_log1p(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_sin(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_cos(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_tan(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_asin(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_acos(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_atan(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_sinh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_cosh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_tanh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_asinh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_acosh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_atanh(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_re(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet1_FromValue(Py_TYPE(self), jet1_im(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_real(PyObject *self, void *closure)
{
    return pyjet1_re(self, NULL);
}

static PyObject *
pyjet1_imag(PyObject *self, void *closure)
{
    return pyjet1_im(self, NULL);
}

static PyObject *pyjet1_grades(PyObject *self, void *closure);

static PyGetSetDef pyjet1_getset[] = {
    {"grade", pyjet1_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pyjet1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyjet1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyjet1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(jet1_is_real(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet1_to_re(((PyJet1 *)self)->v));
}

static PyObject *
pyjet1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyJet1_FromValue((PyTypeObject *)cls, jet1_from_re(a));
}

static PyObject *
pyjet1_lpow(PyObject *self, PyObject *other)
{
    jet1 b;
    int r;
    r = jet1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet1_FromValue(Py_TYPE(self),
        jet1_lpow(((PyJet1 *)self)->v, b));
}

static PyObject *
pyjet1_rpow(PyObject *self, PyObject *other)
{
    jet1 b;
    int r;
    r = jet1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet1_FromValue(Py_TYPE(self),
        jet1_rpow(((PyJet1 *)self)->v, b));
}

static PyObject *
pyjet1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    jet1 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = jet1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet1_FromValue(Py_TYPE(self),
        jet1_mix(((PyJet1 *)self)->v, b, t));
}

static PyObject *
pyjet1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_jet1, 2);
}

static PyObject *
pyjet1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_jet1, 2);
}

static PyObject *
pyjet1_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_JET1(jet1_ipow(((PyJet1 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyjet1_grade(PyObject *self, Py_ssize_t k)
{
    jet1 v = ((PyJet1 *)self)->v;
    switch (k) {
        case 0:
            return PyJet1_FromValue(Py_TYPE(self), jet1_grade0(v));
        case 1:
            return PyJet1_FromValue(Py_TYPE(self), jet1_grade1(v));
        default:
            return PyJet1_FromValue(Py_TYPE(self), jet1_zero());
    }
}

static PyObject *
pyjet1_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyjet1_grade, 2);
}

static PyMethodDef pyjet1_methods[] = {
    {"to_float", (PyCFunction)pyjet1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyjet1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyjet1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyjet1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyjet1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyjet1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyjet1_square, METH_NOARGS, "square."},
    {"expm1", pyjet1_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyjet1_log1p, METH_NOARGS, "log1p."},
    {"sin", pyjet1_sin, METH_NOARGS, "sin."},
    {"cos", pyjet1_cos, METH_NOARGS, "cos."},
    {"tan", pyjet1_tan, METH_NOARGS, "tan."},
    {"asin", pyjet1_asin, METH_NOARGS, "asin."},
    {"acos", pyjet1_acos, METH_NOARGS, "acos."},
    {"atan", pyjet1_atan, METH_NOARGS, "atan."},
    {"sinh", pyjet1_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyjet1_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyjet1_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyjet1_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyjet1_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyjet1_atanh, METH_NOARGS, "atanh."},
    {"re", pyjet1_re, METH_NOARGS, "re."},
    {"im", pyjet1_im, METH_NOARGS, "im."},
    {"arcsin", pyjet1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyjet1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyjet1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyjet1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyjet1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyjet1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyjet1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyjet1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyjet1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pyjet1_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pyjet1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyjet1_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyjet1_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyjet1_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyjet1_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyjet1_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyjet1_conjugate, METH_NOARGS, "Conjugate e -> -e."},
    {"lmat", pyjet1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyjet1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyjet1_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyjet1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyjet1_exp, METH_NOARGS, "exp."},
    {"log", pyjet1_log, METH_NOARGS, "log."},
    {"sqrt", pyjet1_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyjet1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyjet1_length(PyObject *Py_UNUSED(self))
{
    return 2;
}

static PyObject *
pyjet1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyJet1 *)self)->v.data[i]);
}

static PyObject *
pyjet1_subscript(PyObject *self, PyObject *key)
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
    return pyjet1_item(self, i);
}

static int
pyjet1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyjet1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyjet1_slots[] = {
    {Py_tp_dealloc, pyjet1_dealloc},
    {Py_tp_traverse, pyjet1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyjet1_repr},
    {Py_tp_str, pyjet1_repr},
    {Py_tp_hash, pyjet1_hash},
    {Py_tp_doc, (void *)"Jet of order 1 in R[e]/(e^2) with 2 float64 components"},
    {Py_tp_richcompare, pyjet1_richcompare},
    {Py_tp_methods, pyjet1_methods},
    {Py_tp_getset, pyjet1_getset},
    {Py_sq_length, pyjet1_length},
    {Py_sq_item, pyjet1_item},
    {Py_mp_subscript, pyjet1_subscript},
    {Py_tp_new, pyjet1_new},
    {Py_nb_add, pyjet1_add},
    {Py_nb_subtract, pyjet1_sub},
    {Py_nb_multiply, pyjet1_mul},
    {Py_nb_true_divide, pyjet1_div},
    {Py_nb_power, (void *)pyjet1_pow},
    {Py_nb_positive, pyjet1_pos},
    {Py_nb_negative, pyjet1_neg},
    {Py_nb_bool, pyjet1_bool},
    {0, NULL}
};

static PyType_Spec pyjet1_spec = {
    .name = "algpy.algebras.jet.jet1",
    .basicsize = sizeof(PyJet1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyjet1_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyjet1_getitem(void *data, void *arr)
{
    jet1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_JET1(v);
}

static int
npyjet1_setitem(PyObject *item, void *data, void *arr)
{
    jet1 v;
    int r = jet1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to jet1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "jet1 needs exactly 2 components");
            return -1;
        }
        for (int i = 0; i < 2; i++) {
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
pyjet1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(jet1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyjet1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyjet1_copyswapn(dst, sizeof(jet1), src, sizeof(jet1), 1,
                         swap, arr);
}

static npy_bool
pyjet1_nonzero(void *data, void *arr)
{
    jet1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec2_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyjet1_arrfuncs;
typedef struct { char c; jet1 v; } pyjet1_align_test;

static PyArray_DescrProto pyjet1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(jet1), offsetof(pyjet1_align_test, v), 0, 0, 0,
    &pyjet1_arrfuncs,
};

#define JET1_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyjet1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        jet1 *dst = (jet1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = jet1_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define JET1_REAL_TYPES(X) \
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

JET1_REAL_TYPES(JET1_CAST_FROM)

static int
pyjet1_register_casts(void)
{
#define REGISTER_JET1_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_jet1,                    \
                                     pyjet1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_jet1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    JET1_REAL_TYPES(REGISTER_JET1_CAST)
#undef REGISTER_JET1_CAST
    return 0;
}

typedef jet1 (*jet1_unary_op)(jet1);
typedef jet1 (*jet1_binary_op)(jet1, jet1);
typedef int (*jet1_predicate)(jet1, jet1);

static void
loop_jet1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    jet1_unary_op op = (jet1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_jet1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    jet1_binary_op op = (jet1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet1 a, b, r;
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
loop_jet1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    jet1_predicate op = (jet1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyjet1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_jet1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_jet1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyjet1_register_ufuncs(PyObject *numpy)
{
    return
        pyjet1_register_loop(numpy, "add", loop_jet1_binary, 3, 0, jet1_add) ||
        pyjet1_register_loop(numpy, "subtract", loop_jet1_binary, 3, 0, jet1_sub) ||
        pyjet1_register_loop(numpy, "multiply", loop_jet1_binary, 3, 0, jet1_mul) ||
        pyjet1_register_loop(numpy, "divide", loop_jet1_binary, 3, 0, jet1_div) ||
        pyjet1_register_loop(numpy, "power", loop_jet1_binary, 3, 0, jet1_rpow) ||
        pyjet1_register_loop(numpy, "positive", loop_jet1_unary, 2, 0, jet1_pos) ||
        pyjet1_register_loop(numpy, "negative", loop_jet1_unary, 2, 0, jet1_neg) ||
        pyjet1_register_loop(numpy, "conjugate", loop_jet1_unary, 2, 0, jet1_conjugate) ||
        pyjet1_register_loop(numpy, "reciprocal", loop_jet1_unary, 2, 0, jet1_inv) ||
        pyjet1_register_loop(numpy, "square", loop_jet1_unary, 2, 0, jet1_square) ||
        pyjet1_register_loop(numpy, "exp", loop_jet1_unary, 2, 0, jet1_exp) ||
        pyjet1_register_loop(numpy, "expm1", loop_jet1_unary, 2, 0, jet1_expm1) ||
        pyjet1_register_loop(numpy, "log", loop_jet1_unary, 2, 0, jet1_log) ||
        pyjet1_register_loop(numpy, "log1p", loop_jet1_unary, 2, 0, jet1_log1p) ||
        pyjet1_register_loop(numpy, "sqrt", loop_jet1_unary, 2, 0, jet1_sqrt) ||
        pyjet1_register_loop(numpy, "sin", loop_jet1_unary, 2, 0, jet1_sin) ||
        pyjet1_register_loop(numpy, "cos", loop_jet1_unary, 2, 0, jet1_cos) ||
        pyjet1_register_loop(numpy, "tan", loop_jet1_unary, 2, 0, jet1_tan) ||
        pyjet1_register_loop(numpy, "arcsin", loop_jet1_unary, 2, 0, jet1_asin) ||
        pyjet1_register_loop(numpy, "arccos", loop_jet1_unary, 2, 0, jet1_acos) ||
        pyjet1_register_loop(numpy, "arctan", loop_jet1_unary, 2, 0, jet1_atan) ||
        pyjet1_register_loop(numpy, "sinh", loop_jet1_unary, 2, 0, jet1_sinh) ||
        pyjet1_register_loop(numpy, "cosh", loop_jet1_unary, 2, 0, jet1_cosh) ||
        pyjet1_register_loop(numpy, "tanh", loop_jet1_unary, 2, 0, jet1_tanh) ||
        pyjet1_register_loop(numpy, "arcsinh", loop_jet1_unary, 2, 0, jet1_asinh) ||
        pyjet1_register_loop(numpy, "arccosh", loop_jet1_unary, 2, 0, jet1_acosh) ||
        pyjet1_register_loop(numpy, "arctanh", loop_jet1_unary, 2, 0, jet1_atanh) ||
        pyjet1_register_loop(numpy, "equal", loop_jet1_predicate, 3, 1, jet1_equal) ||
        pyjet1_register_loop(numpy, "not_equal", loop_jet1_predicate, 3, 1, jet1_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_jet1", "Jet type of order 1 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__jet1(void)
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
    PyJet1_Type = PyType_FromModuleAndSpec(
        m, &pyjet1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyJet1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyjet1_arrfuncs);
    pyjet1_arrfuncs.getitem = npyjet1_getitem;
    pyjet1_arrfuncs.setitem = npyjet1_setitem;
    pyjet1_arrfuncs.copyswapn = pyjet1_copyswapn;
    pyjet1_arrfuncs.copyswap = pyjet1_copyswap;
    pyjet1_arrfuncs.nonzero = pyjet1_nonzero;
    pyjet1_descr_proto.typeobj = (PyTypeObject *)PyJet1_Type;
    Py_SET_TYPE((PyObject *)&pyjet1_descr_proto, &PyArrayDescr_Type);
    npy_jet1 = PyArray_RegisterDataType(&pyjet1_descr_proto);
    if (npy_jet1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_jet1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyJet1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyJet1_Type, 2) < 0) {
        goto fail;
    }
    if (pyjet1_register_casts() < 0 ||
            pyjet1_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "jet1", PyJet1_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
