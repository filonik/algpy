/*
 * Jets of order 3, R[e]/(e^4), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in jet3.h.
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

#include "jet3.h"

typedef struct {
    PyObject_HEAD
    jet3 v;
} PyJet3;

static PyObject *PyJet3_Type = NULL;
static int npy_jet3 = -1;

static PyObject *
PyJet3_FromValue(PyTypeObject *type, jet3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyJet3 *p = (PyJet3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_JET3(v) \
    PyJet3_FromValue((PyTypeObject *)PyJet3_Type, (v))

static int
jet3_from_object(PyObject *o, jet3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyJet3_Type)) {
        *out = ((PyJet3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = jet3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyjet3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x1", "x2", "x3", NULL};
    double d[4] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                     &d[0], &d[1], &d[2], &d[3])) {
        return NULL;
    }
    return PyJet3_FromValue(type, jet3_make(d));
}

static PyObject *
pyjet3_repr(PyObject *self)
{
    jet3 v = ((PyJet3 *)self)->v;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("jet3(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyjet3_hash(PyObject *self)
{
    jet3 v = ((PyJet3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_JET3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (jet3_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddd)", v.data[0], v.data[1], v.data[2], v.data[3]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyjet3_richcompare(PyObject *a, PyObject *b, int op)
{
    jet3 x, y;
    int r;
    if ((r = jet3_from_object(a, &x)) <= 0 ||
            (r = jet3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(jet3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!jet3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define JET3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyjet3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        jet3 x, y;                                                       \
        int r;                                                               \
        if ((r = jet3_from_object(a, &x)) <= 0 ||                        \
                (r = jet3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_JET3(expr(x, y));                                     \
    }

JET3_BINARY(add, jet3_add)
JET3_BINARY(sub, jet3_sub)
JET3_BINARY(mul, jet3_mul)
JET3_BINARY(div, jet3_div)

static PyObject *
pyjet3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    jet3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "jet3 does not support modular power");
        return NULL;
    }
    if ((r = jet3_from_object(a, &x)) <= 0 ||
            (r = jet3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_JET3(jet3_rpow(x, y));
}

#define JET3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyjet3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_JET3(expr(((PyJet3 *)self)->v));                  \
    }

JET3_UNARY(pos, jet3_pos)
JET3_UNARY(neg, jet3_neg)

static PyObject *
pyjet3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue((PyTypeObject *)cls, jet3_zero());
}

static PyObject *
pyjet3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue((PyTypeObject *)cls, jet3_one());
}

static int
pyjet3_bool(PyObject *self)
{
    jet3 v = ((PyJet3 *)self)->v;
    return vec4_any_nonzero(v);
}

#define JET3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyjet3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_JET3(op(((PyJet3 *)self)->v));                   \
    }

JET3_METHOD(conjugate, jet3_conjugate)
JET3_METHOD(exp, jet3_exp)
JET3_METHOD(log, jet3_log)
JET3_METHOD(sqrt, jet3_sqrt)

static PyObject *
pyjet3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet3_norm(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_JET3(jet3_inv(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet3_lmat(((PyJet3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet3_rmat(((PyJet3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet3_minpoly(((PyJet3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    jet3_eivals(((PyJet3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyjet3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet3);
    jet3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet3 *)PyArray_DATA(arr);
    jet3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet3_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet3);
    jet3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet3 *)PyArray_DATA(arr);
    jet3_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_jet3);
    jet3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (jet3 *)PyArray_DATA(arr);
    jet3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyjet3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    jet3 v = ((PyJet3 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyObject *
pyjet3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_square(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_expm1(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_log1p(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_sin(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_cos(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_tan(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_asin(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_acos(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_atan(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_sinh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_cosh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_tanh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_asinh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_acosh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_atanh(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_re(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet3_FromValue(Py_TYPE(self), jet3_im(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_real(PyObject *self, void *closure)
{
    return pyjet3_re(self, NULL);
}

static PyObject *
pyjet3_imag(PyObject *self, void *closure)
{
    return pyjet3_im(self, NULL);
}

static PyObject *pyjet3_grades(PyObject *self, void *closure);

static PyGetSetDef pyjet3_getset[] = {
    {"grade", pyjet3_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pyjet3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyjet3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyjet3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(jet3_is_real(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet3_to_re(((PyJet3 *)self)->v));
}

static PyObject *
pyjet3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyJet3_FromValue((PyTypeObject *)cls, jet3_from_re(a));
}

static PyObject *
pyjet3_lpow(PyObject *self, PyObject *other)
{
    jet3 b;
    int r;
    r = jet3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet3_FromValue(Py_TYPE(self),
        jet3_lpow(((PyJet3 *)self)->v, b));
}

static PyObject *
pyjet3_rpow(PyObject *self, PyObject *other)
{
    jet3 b;
    int r;
    r = jet3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet3_FromValue(Py_TYPE(self),
        jet3_rpow(((PyJet3 *)self)->v, b));
}

static PyObject *
pyjet3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    jet3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = jet3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet3_FromValue(Py_TYPE(self),
        jet3_mix(((PyJet3 *)self)->v, b, t));
}

static PyObject *
pyjet3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_jet3, 4);
}

static PyObject *
pyjet3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_jet3, 4);
}

static PyObject *
pyjet3_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_JET3(jet3_ipow(((PyJet3 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyjet3_grade(PyObject *self, Py_ssize_t k)
{
    jet3 v = ((PyJet3 *)self)->v;
    switch (k) {
        case 0:
            return PyJet3_FromValue(Py_TYPE(self), jet3_grade0(v));
        case 1:
            return PyJet3_FromValue(Py_TYPE(self), jet3_grade1(v));
        case 2:
            return PyJet3_FromValue(Py_TYPE(self), jet3_grade2(v));
        case 3:
            return PyJet3_FromValue(Py_TYPE(self), jet3_grade3(v));
        default:
            return PyJet3_FromValue(Py_TYPE(self), jet3_zero());
    }
}

static PyObject *
pyjet3_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyjet3_grade, 4);
}

static PyMethodDef pyjet3_methods[] = {
    {"to_float", (PyCFunction)pyjet3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyjet3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyjet3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyjet3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyjet3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyjet3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyjet3_square, METH_NOARGS, "square."},
    {"expm1", pyjet3_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyjet3_log1p, METH_NOARGS, "log1p."},
    {"sin", pyjet3_sin, METH_NOARGS, "sin."},
    {"cos", pyjet3_cos, METH_NOARGS, "cos."},
    {"tan", pyjet3_tan, METH_NOARGS, "tan."},
    {"asin", pyjet3_asin, METH_NOARGS, "asin."},
    {"acos", pyjet3_acos, METH_NOARGS, "acos."},
    {"atan", pyjet3_atan, METH_NOARGS, "atan."},
    {"sinh", pyjet3_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyjet3_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyjet3_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyjet3_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyjet3_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyjet3_atanh, METH_NOARGS, "atanh."},
    {"re", pyjet3_re, METH_NOARGS, "re."},
    {"im", pyjet3_im, METH_NOARGS, "im."},
    {"arcsin", pyjet3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyjet3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyjet3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyjet3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyjet3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyjet3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyjet3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyjet3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyjet3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pyjet3_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pyjet3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyjet3_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyjet3_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyjet3_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyjet3_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyjet3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyjet3_conjugate, METH_NOARGS, "Conjugate e -> -e."},
    {"lmat", pyjet3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyjet3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyjet3_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyjet3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyjet3_exp, METH_NOARGS, "exp."},
    {"log", pyjet3_log, METH_NOARGS, "log."},
    {"sqrt", pyjet3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyjet3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyjet3_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pyjet3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyJet3 *)self)->v.data[i]);
}

static PyObject *
pyjet3_subscript(PyObject *self, PyObject *key)
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
        i += 4;
    }
    return pyjet3_item(self, i);
}

static int
pyjet3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyjet3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyjet3_slots[] = {
    {Py_tp_dealloc, pyjet3_dealloc},
    {Py_tp_traverse, pyjet3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyjet3_repr},
    {Py_tp_str, pyjet3_repr},
    {Py_tp_hash, pyjet3_hash},
    {Py_tp_doc, (void *)"Jet of order 3 in R[e]/(e^4) with 4 float64 components"},
    {Py_tp_richcompare, pyjet3_richcompare},
    {Py_tp_methods, pyjet3_methods},
    {Py_tp_getset, pyjet3_getset},
    {Py_sq_length, pyjet3_length},
    {Py_sq_item, pyjet3_item},
    {Py_mp_subscript, pyjet3_subscript},
    {Py_tp_new, pyjet3_new},
    {Py_nb_add, pyjet3_add},
    {Py_nb_subtract, pyjet3_sub},
    {Py_nb_multiply, pyjet3_mul},
    {Py_nb_true_divide, pyjet3_div},
    {Py_nb_power, (void *)pyjet3_pow},
    {Py_nb_positive, pyjet3_pos},
    {Py_nb_negative, pyjet3_neg},
    {Py_nb_bool, pyjet3_bool},
    {0, NULL}
};

static PyType_Spec pyjet3_spec = {
    .name = "algpy.algebras.jet.jet3",
    .basicsize = sizeof(PyJet3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyjet3_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyjet3_getitem(void *data, void *arr)
{
    jet3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_JET3(v);
}

static int
npyjet3_setitem(PyObject *item, void *data, void *arr)
{
    jet3 v;
    int r = jet3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to jet3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "jet3 needs exactly 4 components");
            return -1;
        }
        for (int i = 0; i < 4; i++) {
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
pyjet3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(jet3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyjet3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyjet3_copyswapn(dst, sizeof(jet3), src, sizeof(jet3), 1,
                         swap, arr);
}

static npy_bool
pyjet3_nonzero(void *data, void *arr)
{
    jet3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyjet3_arrfuncs;
typedef struct { char c; jet3 v; } pyjet3_align_test;

static PyArray_DescrProto pyjet3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(jet3), offsetof(pyjet3_align_test, v), 0, 0, 0,
    &pyjet3_arrfuncs,
};

#define JET3_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyjet3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        jet3 *dst = (jet3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = jet3_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define JET3_REAL_TYPES(X) \
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

JET3_REAL_TYPES(JET3_CAST_FROM)

static int
pyjet3_register_casts(void)
{
#define REGISTER_JET3_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_jet3,                    \
                                     pyjet3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_jet3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    JET3_REAL_TYPES(REGISTER_JET3_CAST)
#undef REGISTER_JET3_CAST
    return 0;
}

typedef jet3 (*jet3_unary_op)(jet3);
typedef jet3 (*jet3_binary_op)(jet3, jet3);
typedef int (*jet3_predicate)(jet3, jet3);

static void
loop_jet3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    jet3_unary_op op = (jet3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_jet3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    jet3_binary_op op = (jet3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet3 a, b, r;
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
loop_jet3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    jet3_predicate op = (jet3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyjet3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_jet3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_jet3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyjet3_register_ufuncs(PyObject *numpy)
{
    return
        pyjet3_register_loop(numpy, "add", loop_jet3_binary, 3, 0, jet3_add) ||
        pyjet3_register_loop(numpy, "subtract", loop_jet3_binary, 3, 0, jet3_sub) ||
        pyjet3_register_loop(numpy, "multiply", loop_jet3_binary, 3, 0, jet3_mul) ||
        pyjet3_register_loop(numpy, "divide", loop_jet3_binary, 3, 0, jet3_div) ||
        pyjet3_register_loop(numpy, "power", loop_jet3_binary, 3, 0, jet3_rpow) ||
        pyjet3_register_loop(numpy, "positive", loop_jet3_unary, 2, 0, jet3_pos) ||
        pyjet3_register_loop(numpy, "negative", loop_jet3_unary, 2, 0, jet3_neg) ||
        pyjet3_register_loop(numpy, "conjugate", loop_jet3_unary, 2, 0, jet3_conjugate) ||
        pyjet3_register_loop(numpy, "reciprocal", loop_jet3_unary, 2, 0, jet3_inv) ||
        pyjet3_register_loop(numpy, "square", loop_jet3_unary, 2, 0, jet3_square) ||
        pyjet3_register_loop(numpy, "exp", loop_jet3_unary, 2, 0, jet3_exp) ||
        pyjet3_register_loop(numpy, "expm1", loop_jet3_unary, 2, 0, jet3_expm1) ||
        pyjet3_register_loop(numpy, "log", loop_jet3_unary, 2, 0, jet3_log) ||
        pyjet3_register_loop(numpy, "log1p", loop_jet3_unary, 2, 0, jet3_log1p) ||
        pyjet3_register_loop(numpy, "sqrt", loop_jet3_unary, 2, 0, jet3_sqrt) ||
        pyjet3_register_loop(numpy, "sin", loop_jet3_unary, 2, 0, jet3_sin) ||
        pyjet3_register_loop(numpy, "cos", loop_jet3_unary, 2, 0, jet3_cos) ||
        pyjet3_register_loop(numpy, "tan", loop_jet3_unary, 2, 0, jet3_tan) ||
        pyjet3_register_loop(numpy, "arcsin", loop_jet3_unary, 2, 0, jet3_asin) ||
        pyjet3_register_loop(numpy, "arccos", loop_jet3_unary, 2, 0, jet3_acos) ||
        pyjet3_register_loop(numpy, "arctan", loop_jet3_unary, 2, 0, jet3_atan) ||
        pyjet3_register_loop(numpy, "sinh", loop_jet3_unary, 2, 0, jet3_sinh) ||
        pyjet3_register_loop(numpy, "cosh", loop_jet3_unary, 2, 0, jet3_cosh) ||
        pyjet3_register_loop(numpy, "tanh", loop_jet3_unary, 2, 0, jet3_tanh) ||
        pyjet3_register_loop(numpy, "arcsinh", loop_jet3_unary, 2, 0, jet3_asinh) ||
        pyjet3_register_loop(numpy, "arccosh", loop_jet3_unary, 2, 0, jet3_acosh) ||
        pyjet3_register_loop(numpy, "arctanh", loop_jet3_unary, 2, 0, jet3_atanh) ||
        pyjet3_register_loop(numpy, "equal", loop_jet3_predicate, 3, 1, jet3_equal) ||
        pyjet3_register_loop(numpy, "not_equal", loop_jet3_predicate, 3, 1, jet3_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_jet3", "Jet type of order 3 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__jet3(void)
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
    PyJet3_Type = PyType_FromModuleAndSpec(
        m, &pyjet3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyJet3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyjet3_arrfuncs);
    pyjet3_arrfuncs.getitem = npyjet3_getitem;
    pyjet3_arrfuncs.setitem = npyjet3_setitem;
    pyjet3_arrfuncs.copyswapn = pyjet3_copyswapn;
    pyjet3_arrfuncs.copyswap = pyjet3_copyswap;
    pyjet3_arrfuncs.nonzero = pyjet3_nonzero;
    pyjet3_descr_proto.typeobj = (PyTypeObject *)PyJet3_Type;
    Py_SET_TYPE((PyObject *)&pyjet3_descr_proto, &PyArrayDescr_Type);
    npy_jet3 = PyArray_RegisterDataType(&pyjet3_descr_proto);
    if (npy_jet3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_jet3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyJet3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyJet3_Type, 4) < 0) {
        goto fail;
    }
    if (pyjet3_register_casts() < 0 ||
            pyjet3_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "jet3", PyJet3_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
