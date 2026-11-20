/*
 * Jets of order 4, R[e]/(e^5), exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in jet4.h.
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

#include "jet4.h"

typedef struct {
    PyObject_HEAD
    jet4 v;
} PyJet4;

static PyObject *PyJet4_Type = NULL;
static int npy_jet4 = -1;

static PyObject *
PyJet4_FromValue(PyTypeObject *type, jet4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyJet4 *p = (PyJet4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_JET4(v) \
    PyJet4_FromValue((PyTypeObject *)PyJet4_Type, (v))

static int
jet4_from_object(PyObject *o, jet4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyJet4_Type)) {
        *out = ((PyJet4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = jet4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pyjet4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x1", "x2", "x3", "x4", NULL};
    double d[5] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ddddd", kwlist,
                                     &d[0], &d[1], &d[2], &d[3], &d[4])) {
        return NULL;
    }
    return PyJet4_FromValue(type, jet4_make(d));
}

static PyObject *
pyjet4_repr(PyObject *self)
{
    jet4 v = ((PyJet4 *)self)->v;
    char *s[5] = {NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 5; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("jet4(%s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4]);
done:
    for (i = 0; i < 5; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyjet4_hash(PyObject *self)
{
    jet4 v = ((PyJet4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_JET4_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (jet4_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(ddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyjet4_richcompare(PyObject *a, PyObject *b, int op)
{
    jet4 x, y;
    int r;
    if ((r = jet4_from_object(a, &x)) <= 0 ||
            (r = jet4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(jet4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!jet4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define JET4_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pyjet4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        jet4 x, y;                                                       \
        int r;                                                               \
        if ((r = jet4_from_object(a, &x)) <= 0 ||                        \
                (r = jet4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_JET4(expr(x, y));                                     \
    }

JET4_BINARY(add, jet4_add)
JET4_BINARY(sub, jet4_sub)
JET4_BINARY(mul, jet4_mul)
JET4_BINARY(div, jet4_div)

static PyObject *
pyjet4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    jet4 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "jet4 does not support modular power");
        return NULL;
    }
    if ((r = jet4_from_object(a, &x)) <= 0 ||
            (r = jet4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_JET4(jet4_rpow(x, y));
}

#define JET4_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pyjet4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_JET4(expr(((PyJet4 *)self)->v));                  \
    }

JET4_UNARY(pos, jet4_pos)
JET4_UNARY(neg, jet4_neg)

static PyObject *
pyjet4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue((PyTypeObject *)cls, jet4_zero());
}

static PyObject *
pyjet4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue((PyTypeObject *)cls, jet4_one());
}

static int
pyjet4_bool(PyObject *self)
{
    jet4 v = ((PyJet4 *)self)->v;
    return vec5_any_nonzero(v);
}

#define JET4_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pyjet4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_JET4(op(((PyJet4 *)self)->v));                   \
    }

JET4_METHOD(conjugate, jet4_conjugate)
JET4_METHOD(exp, jet4_exp)
JET4_METHOD(log, jet4_log)
JET4_METHOD(sqrt, jet4_sqrt)

static PyObject *
pyjet4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet4_norm(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_JET4(jet4_inv(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {5, 5};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet4_lmat(((PyJet4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {5, 5};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet4_rmat(((PyJet4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {6};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    jet4_minpoly(((PyJet4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pyjet4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[5];
    npy_intp dims[1] = {5};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    jet4_eivals(((PyJet4 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 5; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pyjet4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet4);
    jet4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet4 *)PyArray_DATA(arr);
    jet4_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet4_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_jet4);
    jet4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (jet4 *)PyArray_DATA(arr);
    jet4_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pyjet4_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {5, 5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_jet4);
    jet4 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (jet4 *)PyArray_DATA(arr);
    jet4_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pyjet4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    jet4 v = ((PyJet4 *)self)->v;
    return Py_BuildValue("(O(ddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4]);
}

static PyObject *
pyjet4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_square(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_expm1(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_log1p(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_sin(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_cos(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_tan(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_asin(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_acos(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_atan(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_sinh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_cosh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_tanh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_asinh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_acosh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_atanh(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_re(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyJet4_FromValue(Py_TYPE(self), jet4_im(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_real(PyObject *self, void *closure)
{
    return pyjet4_re(self, NULL);
}

static PyObject *
pyjet4_imag(PyObject *self, void *closure)
{
    return pyjet4_im(self, NULL);
}

static PyObject *pyjet4_grades(PyObject *self, void *closure);

static PyGetSetDef pyjet4_getset[] = {
    {"grade", pyjet4_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pyjet4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pyjet4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyjet4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(jet4_is_real(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(jet4_to_re(((PyJet4 *)self)->v));
}

static PyObject *
pyjet4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyJet4_FromValue((PyTypeObject *)cls, jet4_from_re(a));
}

static PyObject *
pyjet4_lpow(PyObject *self, PyObject *other)
{
    jet4 b;
    int r;
    r = jet4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet4_FromValue(Py_TYPE(self),
        jet4_lpow(((PyJet4 *)self)->v, b));
}

static PyObject *
pyjet4_rpow(PyObject *self, PyObject *other)
{
    jet4 b;
    int r;
    r = jet4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet4_FromValue(Py_TYPE(self),
        jet4_rpow(((PyJet4 *)self)->v, b));
}

static PyObject *
pyjet4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    jet4 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = jet4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyJet4_FromValue(Py_TYPE(self),
        jet4_mix(((PyJet4 *)self)->v, b, t));
}

static PyObject *
pyjet4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_jet4, 5);
}

static PyObject *
pyjet4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_jet4, 5);
}

static PyObject *
pyjet4_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_JET4(jet4_ipow(((PyJet4 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyjet4_grade(PyObject *self, Py_ssize_t k)
{
    jet4 v = ((PyJet4 *)self)->v;
    switch (k) {
        case 0:
            return PyJet4_FromValue(Py_TYPE(self), jet4_grade0(v));
        case 1:
            return PyJet4_FromValue(Py_TYPE(self), jet4_grade1(v));
        case 2:
            return PyJet4_FromValue(Py_TYPE(self), jet4_grade2(v));
        case 3:
            return PyJet4_FromValue(Py_TYPE(self), jet4_grade3(v));
        case 4:
            return PyJet4_FromValue(Py_TYPE(self), jet4_grade4(v));
        default:
            return PyJet4_FromValue(Py_TYPE(self), jet4_zero());
    }
}

static PyObject *
pyjet4_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyjet4_grade, 5);
}

static PyMethodDef pyjet4_methods[] = {
    {"to_float", (PyCFunction)pyjet4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyjet4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pyjet4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pyjet4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pyjet4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pyjet4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pyjet4_square, METH_NOARGS, "square."},
    {"expm1", pyjet4_expm1, METH_NOARGS, "expm1."},
    {"log1p", pyjet4_log1p, METH_NOARGS, "log1p."},
    {"sin", pyjet4_sin, METH_NOARGS, "sin."},
    {"cos", pyjet4_cos, METH_NOARGS, "cos."},
    {"tan", pyjet4_tan, METH_NOARGS, "tan."},
    {"asin", pyjet4_asin, METH_NOARGS, "asin."},
    {"acos", pyjet4_acos, METH_NOARGS, "acos."},
    {"atan", pyjet4_atan, METH_NOARGS, "atan."},
    {"sinh", pyjet4_sinh, METH_NOARGS, "sinh."},
    {"cosh", pyjet4_cosh, METH_NOARGS, "cosh."},
    {"tanh", pyjet4_tanh, METH_NOARGS, "tanh."},
    {"asinh", pyjet4_asinh, METH_NOARGS, "asinh."},
    {"acosh", pyjet4_acosh, METH_NOARGS, "acosh."},
    {"atanh", pyjet4_atanh, METH_NOARGS, "atanh."},
    {"re", pyjet4_re, METH_NOARGS, "re."},
    {"im", pyjet4_im, METH_NOARGS, "im."},
    {"arcsin", pyjet4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pyjet4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pyjet4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pyjet4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pyjet4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pyjet4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pyjet4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyjet4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyjet4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pyjet4_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pyjet4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pyjet4_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pyjet4_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pyjet4_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pyjet4_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pyjet4_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pyjet4_conjugate, METH_NOARGS, "Conjugate e -> -e."},
    {"lmat", pyjet4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pyjet4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pyjet4_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pyjet4_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pyjet4_exp, METH_NOARGS, "exp."},
    {"log", pyjet4_log, METH_NOARGS, "log."},
    {"sqrt", pyjet4_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pyjet4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyjet4_length(PyObject *Py_UNUSED(self))
{
    return 5;
}

static PyObject *
pyjet4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 5) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyJet4 *)self)->v.data[i]);
}

static PyObject *
pyjet4_subscript(PyObject *self, PyObject *key)
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
        i += 5;
    }
    return pyjet4_item(self, i);
}

static int
pyjet4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyjet4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyjet4_slots[] = {
    {Py_tp_dealloc, pyjet4_dealloc},
    {Py_tp_traverse, pyjet4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyjet4_repr},
    {Py_tp_str, pyjet4_repr},
    {Py_tp_hash, pyjet4_hash},
    {Py_tp_doc, (void *)"Jet of order 4 in R[e]/(e^5) with 5 float64 components"},
    {Py_tp_richcompare, pyjet4_richcompare},
    {Py_tp_methods, pyjet4_methods},
    {Py_tp_getset, pyjet4_getset},
    {Py_sq_length, pyjet4_length},
    {Py_sq_item, pyjet4_item},
    {Py_mp_subscript, pyjet4_subscript},
    {Py_tp_new, pyjet4_new},
    {Py_nb_add, pyjet4_add},
    {Py_nb_subtract, pyjet4_sub},
    {Py_nb_multiply, pyjet4_mul},
    {Py_nb_true_divide, pyjet4_div},
    {Py_nb_power, (void *)pyjet4_pow},
    {Py_nb_positive, pyjet4_pos},
    {Py_nb_negative, pyjet4_neg},
    {Py_nb_bool, pyjet4_bool},
    {0, NULL}
};

static PyType_Spec pyjet4_spec = {
    .name = "algpy.algebras.jet.jet4",
    .basicsize = sizeof(PyJet4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyjet4_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyjet4_getitem(void *data, void *arr)
{
    jet4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_JET4(v);
}

static int
npyjet4_setitem(PyObject *item, void *data, void *arr)
{
    jet4 v;
    int r = jet4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to jet4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 5) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "jet4 needs exactly 5 components");
            return -1;
        }
        for (int i = 0; i < 5; i++) {
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
pyjet4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(jet4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyjet4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyjet4_copyswapn(dst, sizeof(jet4), src, sizeof(jet4), 1,
                         swap, arr);
}

static npy_bool
pyjet4_nonzero(void *data, void *arr)
{
    jet4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec5_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyjet4_arrfuncs;
typedef struct { char c; jet4 v; } pyjet4_align_test;

static PyArray_DescrProto pyjet4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(jet4), offsetof(pyjet4_align_test, v), 0, 0, 0,
    &pyjet4_arrfuncs,
};

#define JET4_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pyjet4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        jet4 *dst = (jet4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = jet4_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define JET4_REAL_TYPES(X) \
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

JET4_REAL_TYPES(JET4_CAST_FROM)

static int
pyjet4_register_casts(void)
{
#define REGISTER_JET4_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_jet4,                    \
                                     pyjet4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_jet4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    JET4_REAL_TYPES(REGISTER_JET4_CAST)
#undef REGISTER_JET4_CAST
    return 0;
}

typedef jet4 (*jet4_unary_op)(jet4);
typedef jet4 (*jet4_binary_op)(jet4, jet4);
typedef int (*jet4_predicate)(jet4, jet4);

static void
loop_jet4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    jet4_unary_op op = (jet4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_jet4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    jet4_binary_op op = (jet4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet4 a, b, r;
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
loop_jet4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    jet4_predicate op = (jet4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        jet4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pyjet4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_jet4;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_jet4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyjet4_register_ufuncs(PyObject *numpy)
{
    return
        pyjet4_register_loop(numpy, "add", loop_jet4_binary, 3, 0, jet4_add) ||
        pyjet4_register_loop(numpy, "subtract", loop_jet4_binary, 3, 0, jet4_sub) ||
        pyjet4_register_loop(numpy, "multiply", loop_jet4_binary, 3, 0, jet4_mul) ||
        pyjet4_register_loop(numpy, "divide", loop_jet4_binary, 3, 0, jet4_div) ||
        pyjet4_register_loop(numpy, "power", loop_jet4_binary, 3, 0, jet4_rpow) ||
        pyjet4_register_loop(numpy, "positive", loop_jet4_unary, 2, 0, jet4_pos) ||
        pyjet4_register_loop(numpy, "negative", loop_jet4_unary, 2, 0, jet4_neg) ||
        pyjet4_register_loop(numpy, "conjugate", loop_jet4_unary, 2, 0, jet4_conjugate) ||
        pyjet4_register_loop(numpy, "reciprocal", loop_jet4_unary, 2, 0, jet4_inv) ||
        pyjet4_register_loop(numpy, "square", loop_jet4_unary, 2, 0, jet4_square) ||
        pyjet4_register_loop(numpy, "exp", loop_jet4_unary, 2, 0, jet4_exp) ||
        pyjet4_register_loop(numpy, "expm1", loop_jet4_unary, 2, 0, jet4_expm1) ||
        pyjet4_register_loop(numpy, "log", loop_jet4_unary, 2, 0, jet4_log) ||
        pyjet4_register_loop(numpy, "log1p", loop_jet4_unary, 2, 0, jet4_log1p) ||
        pyjet4_register_loop(numpy, "sqrt", loop_jet4_unary, 2, 0, jet4_sqrt) ||
        pyjet4_register_loop(numpy, "sin", loop_jet4_unary, 2, 0, jet4_sin) ||
        pyjet4_register_loop(numpy, "cos", loop_jet4_unary, 2, 0, jet4_cos) ||
        pyjet4_register_loop(numpy, "tan", loop_jet4_unary, 2, 0, jet4_tan) ||
        pyjet4_register_loop(numpy, "arcsin", loop_jet4_unary, 2, 0, jet4_asin) ||
        pyjet4_register_loop(numpy, "arccos", loop_jet4_unary, 2, 0, jet4_acos) ||
        pyjet4_register_loop(numpy, "arctan", loop_jet4_unary, 2, 0, jet4_atan) ||
        pyjet4_register_loop(numpy, "sinh", loop_jet4_unary, 2, 0, jet4_sinh) ||
        pyjet4_register_loop(numpy, "cosh", loop_jet4_unary, 2, 0, jet4_cosh) ||
        pyjet4_register_loop(numpy, "tanh", loop_jet4_unary, 2, 0, jet4_tanh) ||
        pyjet4_register_loop(numpy, "arcsinh", loop_jet4_unary, 2, 0, jet4_asinh) ||
        pyjet4_register_loop(numpy, "arccosh", loop_jet4_unary, 2, 0, jet4_acosh) ||
        pyjet4_register_loop(numpy, "arctanh", loop_jet4_unary, 2, 0, jet4_atanh) ||
        pyjet4_register_loop(numpy, "equal", loop_jet4_predicate, 3, 1, jet4_equal) ||
        pyjet4_register_loop(numpy, "not_equal", loop_jet4_predicate, 3, 1, jet4_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_jet4", "Jet type of order 4 for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__jet4(void)
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
    PyJet4_Type = PyType_FromModuleAndSpec(
        m, &pyjet4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyJet4_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyjet4_arrfuncs);
    pyjet4_arrfuncs.getitem = npyjet4_getitem;
    pyjet4_arrfuncs.setitem = npyjet4_setitem;
    pyjet4_arrfuncs.copyswapn = pyjet4_copyswapn;
    pyjet4_arrfuncs.copyswap = pyjet4_copyswap;
    pyjet4_arrfuncs.nonzero = pyjet4_nonzero;
    pyjet4_descr_proto.typeobj = (PyTypeObject *)PyJet4_Type;
    Py_SET_TYPE((PyObject *)&pyjet4_descr_proto, &PyArrayDescr_Type);
    npy_jet4 = PyArray_RegisterDataType(&pyjet4_descr_proto);
    if (npy_jet4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_jet4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyJet4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyJet4_Type, 5) < 0) {
        goto fail;
    }
    if (pyjet4_register_casts() < 0 ||
            pyjet4_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "jet4", PyJet4_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
