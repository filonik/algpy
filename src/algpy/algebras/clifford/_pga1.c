/*
 * Projective geometric algebra in 1D exposed as a NumPy user-defined dtype.
 * The arithmetic itself lives in pga1.h.
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

#include "pga1.h"

typedef struct {
    PyObject_HEAD
    pga1 v;
} PyPga1;

static PyObject *PyPga1_Type = NULL;
static int npy_pga1 = -1;

static PyObject *
PyPga1_FromValue(PyTypeObject *type, pga1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyPga1 *p = (PyPga1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_PGA1(v) \
    PyPga1_FromValue((PyTypeObject *)PyPga1_Type, (v))

static int
pga1_from_object(PyObject *o, pga1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyPga1_Type)) {
        *out = ((PyPga1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = pga1_make((const double[]){d, 0, 0, 0});
        return 1;
    }
    return 0;
}

static PyObject *
pypga1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[4] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
            &d[0], &d[1], &d[2], &d[3])) {
        return NULL;
    }
    return PyPga1_FromValue(type, pga1_make(d));
}

static PyObject *
pypga1_repr(PyObject *self)
{
    pga1 v = ((PyPga1 *)self)->v;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("pga1(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pypga1_hash(PyObject *self)
{
    pga1 v = ((PyPga1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_PGA1_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (v.data[1] == 0 && v.data[2] == 0 && v.data[3] == 0) {
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
pypga1_richcompare(PyObject *a, PyObject *b, int op)
{
    pga1 x, y;
    int r;
    if ((r = pga1_from_object(a, &x)) <= 0 ||
            (r = pga1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(pga1_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!pga1_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define PGA1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pypga1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        pga1 x, y;                                                       \
        int r;                                                               \
        if ((r = pga1_from_object(a, &x)) <= 0 ||                        \
                (r = pga1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_PGA1(expr(x, y));                                     \
    }

PGA1_BINARY(add, pga1_add)
PGA1_BINARY(sub, pga1_sub)
PGA1_BINARY(mul, pga1_mul)
PGA1_BINARY(div, pga1_div)
PGA1_BINARY(join, pga1_join)
PGA1_BINARY(meet, pga1_meet)

static PyObject *
pypga1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    pga1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "pga1 does not support modular power");
        return NULL;
    }
    if ((r = pga1_from_object(a, &x)) <= 0 ||
            (r = pga1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_PGA1(pga1_rpow(x, y));
}

#define PGA1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pypga1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_PGA1(expr(((PyPga1 *)self)->v));                  \
    }

PGA1_UNARY(pos, pga1_pos)
PGA1_UNARY(neg, pga1_neg)

/* ldual and rdual differ in even dimensions, so ~ is deliberately undefined. */
static PyObject *
pypga1_invert(PyObject *Py_UNUSED(self))
{
    PyErr_SetString(PyExc_TypeError,
                    "~ is undefined for pga1: ldual and rdual differ in "
                    "even dimensions; use ldual() or rdual()");
    return NULL;
}

static PyObject *
pypga1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue((PyTypeObject *)cls, pga1_zero());
}

static PyObject *
pypga1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue((PyTypeObject *)cls, pga1_one());
}

static PyObject *
pypga1_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue((PyTypeObject *)cls, pga1_bottom());
}

static PyObject *
pypga1_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue((PyTypeObject *)cls, pga1_top());
}

static PyObject *
pypga1_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_ldual(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_rdual(((PyPga1 *)self)->v));
}

static int
pypga1_bool(PyObject *self)
{
    pga1 v = ((PyPga1 *)self)->v;
    return vec4_any_nonzero(v);
}

#define PGA1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pypga1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_PGA1(op(((PyPga1 *)self)->v));                   \
    }

PGA1_METHOD(conjugate, pga1_conjugate)
PGA1_METHOD(exp, pga1_exp)
PGA1_METHOD(log, pga1_log)
PGA1_METHOD(sqrt, pga1_sqrt)
PGA1_METHOD(sin, pga1_sin)
PGA1_METHOD(cos, pga1_cos)
PGA1_METHOD(tan, pga1_tan)
PGA1_METHOD(sinh, pga1_sinh)
PGA1_METHOD(cosh, pga1_cosh)
PGA1_METHOD(tanh, pga1_tanh)

static PyObject *
pypga1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = pga1_norm(((PyPga1 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pypga1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_PGA1(pga1_inv(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga1_lmat(((PyPga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga1_rmat(((PyPga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga1_minpoly(((PyPga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[2];
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    pga1_eivals(((PyPga1 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pypga1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga1);
    pga1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga1 *)PyArray_DATA(arr);
    pga1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga1_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga1);
    pga1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga1 *)PyArray_DATA(arr);
    pga1_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_pga1);
    pga1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (pga1 *)PyArray_DATA(arr);
    pga1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pypga1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    pga1 v = ((PyPga1 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyObject *
pypga1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_square(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_expm1(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_log1p(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_asin(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_acos(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_atan(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_asinh(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_acosh(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_atanh(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_re(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga1_FromValue(Py_TYPE(self), pga1_im(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_real(PyObject *self, void *closure)
{
    return pypga1_re(self, NULL);
}

static PyObject *
pypga1_imag(PyObject *self, void *closure)
{
    return pypga1_im(self, NULL);
}

static PyObject *pypga1_grades(PyObject *self, void *closure);

static PyGetSetDef pypga1_getset[] = {
    {"grade", pypga1_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pypga1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pypga1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pypga1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(pga1_is_real(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga1_to_re(((PyPga1 *)self)->v));
}

static PyObject *
pypga1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyPga1_FromValue((PyTypeObject *)cls, pga1_from_re(a));
}

static PyObject *
pypga1_lpow(PyObject *self, PyObject *other)
{
    pga1 b;
    int r;
    r = pga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga1_FromValue(Py_TYPE(self),
        pga1_lpow(((PyPga1 *)self)->v, b));
}

static PyObject *
pypga1_rpow(PyObject *self, PyObject *other)
{
    pga1 b;
    int r;
    r = pga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga1_FromValue(Py_TYPE(self),
        pga1_rpow(((PyPga1 *)self)->v, b));
}

static PyObject *
pypga1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    pga1 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = pga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga1_FromValue(Py_TYPE(self),
        pga1_mix(((PyPga1 *)self)->v, b, t));
}

static PyObject *
pypga1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_pga1, 4);
}

static PyObject *
pypga1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_pga1, 4);
}

static PyObject *
pypga1_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_PGA1(pga1_ipow(((PyPga1 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pypga1_grade(PyObject *self, Py_ssize_t k)
{
    pga1 v = ((PyPga1 *)self)->v;
    switch (k) {
        case 0:
            return PyPga1_FromValue(Py_TYPE(self), pga1_grade0(v));
        case 1:
            return PyPga1_FromValue(Py_TYPE(self), pga1_grade1(v));
        case 2:
            return PyPga1_FromValue(Py_TYPE(self), pga1_grade2(v));
        default:
            return PyPga1_FromValue(Py_TYPE(self), pga1_zero());
    }
}

static PyObject *
pypga1_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pypga1_grade, 3);
}

static PyMethodDef pypga1_methods[] = {
    {"to_float", (PyCFunction)pypga1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pypga1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pypga1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pypga1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pypga1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pypga1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pypga1_square, METH_NOARGS, "square."},
    {"expm1", pypga1_expm1, METH_NOARGS, "expm1."},
    {"log1p", pypga1_log1p, METH_NOARGS, "log1p."},
    {"asin", pypga1_asin, METH_NOARGS, "asin."},
    {"acos", pypga1_acos, METH_NOARGS, "acos."},
    {"atan", pypga1_atan, METH_NOARGS, "atan."},
    {"asinh", pypga1_asinh, METH_NOARGS, "asinh."},
    {"acosh", pypga1_acosh, METH_NOARGS, "acosh."},
    {"atanh", pypga1_atanh, METH_NOARGS, "atanh."},
    {"re", pypga1_re, METH_NOARGS, "re."},
    {"im", pypga1_im, METH_NOARGS, "im."},
    {"arcsin", pypga1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pypga1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pypga1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pypga1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pypga1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pypga1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pypga1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pypga1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pypga1_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pypga1_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pypga1_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pypga1_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"basis", pypga1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pypga1_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pypga1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pypga1_norm, METH_NOARGS, "Norm as a Python complex."},
    {"lpow", pypga1_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pypga1_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pypga1_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pypga1_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pypga1_conjugate, METH_NOARGS, "Conjugate of the e1 part."},
    {"lmat", pypga1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pypga1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pypga1_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pypga1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pypga1_exp, METH_NOARGS, "exp."},
    {"log", pypga1_log, METH_NOARGS, "log."},
    {"sqrt", pypga1_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pypga1_sin, METH_NOARGS, "sin."},
    {"cos", pypga1_cos, METH_NOARGS, "cos."},
    {"tan", pypga1_tan, METH_NOARGS, "tan."},
    {"sinh", pypga1_sinh, METH_NOARGS, "sinh."},
    {"cosh", pypga1_cosh, METH_NOARGS, "cosh."},
    {"tanh", pypga1_tanh, METH_NOARGS, "tanh."},
    {"__reduce__", pypga1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pypga1_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pypga1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyPga1 *)self)->v.data[i]);
}

static PyObject *
pypga1_subscript(PyObject *self, PyObject *key)
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
    return pypga1_item(self, i);
}

static int
pypga1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pypga1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pypga1_slots[] = {
    {Py_tp_dealloc, pypga1_dealloc},
    {Py_tp_traverse, pypga1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pypga1_repr},
    {Py_tp_str, pypga1_repr},
    {Py_tp_hash, pypga1_hash},
    {Py_tp_doc, (void *)"Projective geometric algebra in 1D with 4 float64 coefficients"},
    {Py_tp_richcompare, pypga1_richcompare},
    {Py_tp_methods, pypga1_methods},
    {Py_tp_getset, pypga1_getset},
    {Py_sq_length, pypga1_length},
    {Py_sq_item, pypga1_item},
    {Py_mp_subscript, pypga1_subscript},
    {Py_tp_new, pypga1_new},
    {Py_nb_add, pypga1_add},
    {Py_nb_subtract, pypga1_sub},
    {Py_nb_multiply, pypga1_mul},
    {Py_nb_true_divide, pypga1_div},
    {Py_nb_power, (void *)pypga1_pow},
    {Py_nb_positive, pypga1_pos},
    {Py_nb_negative, pypga1_neg},
    {Py_nb_or, pypga1_join},
    {Py_nb_and, pypga1_meet},
    {Py_nb_invert, pypga1_invert},
    {Py_nb_bool, pypga1_bool},
    {0, NULL}
};

static PyType_Spec pypga1_spec = {
    .name = "algpy.algebras.clifford.pga1",
    .basicsize = sizeof(PyPga1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pypga1_slots,
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
npypga1_getitem(void *data, void *arr)
{
    pga1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_PGA1(v);
}

static int
npypga1_setitem(PyObject *item, void *data, void *arr)
{
    pga1 v;
    int r = pga1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to pga1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "pga1 needs exactly 4 components");
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
pypga1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(pga1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pypga1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pypga1_copyswapn(dst, sizeof(pga1), src, sizeof(pga1), 1,
                         swap, arr);
}

static npy_bool
pypga1_nonzero(void *data, void *arr)
{
    pga1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pypga1_arrfuncs;
typedef struct { char c; pga1 v; } pypga1_align_test;

static PyArray_DescrProto pypga1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(pga1), offsetof(pypga1_align_test, v), 0, 0, 0,
    &pypga1_arrfuncs,
};

#define PGA1_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pypga1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        pga1 *dst = (pga1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = pga1_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define PGA1_REAL_TYPES(X) \
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

PGA1_REAL_TYPES(PGA1_CAST_FROM)

static int
pypga1_register_casts(void)
{
#define REGISTER_PGA1_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_pga1,                    \
                                     pypga1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_pga1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    PGA1_REAL_TYPES(REGISTER_PGA1_CAST)
#undef REGISTER_PGA1_CAST
    return 0;
}

typedef pga1 (*pga1_unary_op)(pga1);
typedef pga1 (*pga1_binary_op)(pga1, pga1);
typedef int (*pga1_predicate)(pga1, pga1);

static void
loop_pga1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    pga1_unary_op op = (pga1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_pga1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    pga1_binary_op op = (pga1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga1 a, b, r;
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
loop_pga1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    pga1_predicate op = (pga1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pypga1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_pga1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_pga1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pypga1_register_ufuncs(PyObject *numpy)
{
    return
        pypga1_register_loop(numpy, "add", loop_pga1_binary, 3, 0, pga1_add) ||
        pypga1_register_loop(numpy, "subtract", loop_pga1_binary, 3, 0, pga1_sub) ||
        pypga1_register_loop(numpy, "multiply", loop_pga1_binary, 3, 0, pga1_mul) ||
        pypga1_register_loop(numpy, "divide", loop_pga1_binary, 3, 0, pga1_div) ||
        pypga1_register_loop(numpy, "power", loop_pga1_binary, 3, 0, pga1_rpow) ||
        pypga1_register_loop(numpy, "positive", loop_pga1_unary, 2, 0, pga1_pos) ||
        pypga1_register_loop(numpy, "negative", loop_pga1_unary, 2, 0, pga1_neg) ||
        pypga1_register_loop(numpy, "conjugate", loop_pga1_unary, 2, 0, pga1_conjugate) ||
        pypga1_register_loop(numpy, "reciprocal", loop_pga1_unary, 2, 0, pga1_inv) ||
        pypga1_register_loop(numpy, "square", loop_pga1_unary, 2, 0, pga1_square) ||
        pypga1_register_loop(numpy, "exp", loop_pga1_unary, 2, 0, pga1_exp) ||
        pypga1_register_loop(numpy, "expm1", loop_pga1_unary, 2, 0, pga1_expm1) ||
        pypga1_register_loop(numpy, "log", loop_pga1_unary, 2, 0, pga1_log) ||
        pypga1_register_loop(numpy, "log1p", loop_pga1_unary, 2, 0, pga1_log1p) ||
        pypga1_register_loop(numpy, "sqrt", loop_pga1_unary, 2, 0, pga1_sqrt) ||
        pypga1_register_loop(numpy, "sin", loop_pga1_unary, 2, 0, pga1_sin) ||
        pypga1_register_loop(numpy, "cos", loop_pga1_unary, 2, 0, pga1_cos) ||
        pypga1_register_loop(numpy, "tan", loop_pga1_unary, 2, 0, pga1_tan) ||
        pypga1_register_loop(numpy, "arcsin", loop_pga1_unary, 2, 0, pga1_asin) ||
        pypga1_register_loop(numpy, "arccos", loop_pga1_unary, 2, 0, pga1_acos) ||
        pypga1_register_loop(numpy, "arctan", loop_pga1_unary, 2, 0, pga1_atan) ||
        pypga1_register_loop(numpy, "sinh", loop_pga1_unary, 2, 0, pga1_sinh) ||
        pypga1_register_loop(numpy, "cosh", loop_pga1_unary, 2, 0, pga1_cosh) ||
        pypga1_register_loop(numpy, "tanh", loop_pga1_unary, 2, 0, pga1_tanh) ||
        pypga1_register_loop(numpy, "arcsinh", loop_pga1_unary, 2, 0, pga1_asinh) ||
        pypga1_register_loop(numpy, "arccosh", loop_pga1_unary, 2, 0, pga1_acosh) ||
        pypga1_register_loop(numpy, "arctanh", loop_pga1_unary, 2, 0, pga1_atanh) ||
        pypga1_register_loop(numpy, "equal", loop_pga1_predicate, 3, 1, pga1_equal) ||
        pypga1_register_loop(numpy, "bitwise_or", loop_pga1_binary, 3, 0, pga1_join) ||
        pypga1_register_loop(numpy, "bitwise_and", loop_pga1_binary, 3, 0, pga1_meet) ||
        pypga1_register_loop(numpy, "not_equal", loop_pga1_predicate, 3, 1, pga1_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pypga1_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pypga1_register_loop(ufuncs, "ldual", loop_pga1_unary, 2, 0, pga1_ldual) ||
        pypga1_register_loop(ufuncs, "rdual", loop_pga1_unary, 2, 0, pga1_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_pga1", "pga1 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__pga1(void)
{
    PyObject *m, *numpy, *ufuncs = NULL;
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
    ufuncs = PyImport_ImportModule("algpy._ufuncs");
    if (ufuncs == NULL) {
        goto fail;
    }
    PyPga1_Type = PyType_FromModuleAndSpec(
        m, &pypga1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyPga1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pypga1_arrfuncs);
    pypga1_arrfuncs.getitem = npypga1_getitem;
    pypga1_arrfuncs.setitem = npypga1_setitem;
    pypga1_arrfuncs.copyswapn = pypga1_copyswapn;
    pypga1_arrfuncs.copyswap = pypga1_copyswap;
    pypga1_arrfuncs.nonzero = pypga1_nonzero;
    pypga1_descr_proto.typeobj = (PyTypeObject *)PyPga1_Type;
    Py_SET_TYPE((PyObject *)&pypga1_descr_proto, &PyArrayDescr_Type);
    npy_pga1 = PyArray_RegisterDataType(&pypga1_descr_proto);
    if (npy_pga1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_pga1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyPga1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyPga1_Type, 4) < 0) {
        goto fail;
    }
    if (pypga1_register_casts() < 0 ||
            pypga1_register_ufuncs(numpy) < 0 ||
            pypga1_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "pga1", PyPga1_Type) < 0) {
        goto fail;
    }
    Py_DECREF(ufuncs);
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(ufuncs);
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
