/*
 * pga3 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in pga3.h.
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

#include "pga3.h"

typedef struct {
    PyObject_HEAD
    pga3 v;
} PyPga3;

static PyObject *PyPga3_Type = NULL;
static int npy_pga3 = -1;

static PyObject *
PyPga3_FromValue(PyTypeObject *type, pga3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyPga3 *p = (PyPga3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_PGA3(v) \
    PyPga3_FromValue((PyTypeObject *)PyPga3_Type, (v))

static int
pga3_from_object(PyObject *o, pga3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyPga3_Type)) {
        *out = ((PyPga3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = pga3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pypga3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PyPga3_FromValue(type, pga3_make(d));
}

static PyObject *
pypga3_repr(PyObject *self)
{
    pga3 v = ((PyPga3 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("pga3(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pypga3_hash(PyObject *self)
{
    pga3 v = ((PyPga3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_PGA3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (pga3_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pypga3_richcompare(PyObject *a, PyObject *b, int op)
{
    pga3 x, y;
    int r;
    if ((r = pga3_from_object(a, &x)) <= 0 ||
            (r = pga3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(pga3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!pga3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define PGA3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pypga3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        pga3 x, y;                                                       \
        int r;                                                               \
        if ((r = pga3_from_object(a, &x)) <= 0 ||                        \
                (r = pga3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_PGA3(expr(x, y));                                     \
    }

PGA3_BINARY(add, pga3_add)
PGA3_BINARY(sub, pga3_sub)
PGA3_BINARY(mul, pga3_mul)
PGA3_BINARY(div, pga3_div)
PGA3_BINARY(join, pga3_join)
PGA3_BINARY(meet, pga3_meet)

static PyObject *
pypga3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    pga3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "pga3 does not support modular power");
        return NULL;
    }
    if ((r = pga3_from_object(a, &x)) <= 0 ||
            (r = pga3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_PGA3(pga3_rpow(x, y));
}

#define PGA3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pypga3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_PGA3(expr(((PyPga3 *)self)->v));                  \
    }

PGA3_UNARY(pos, pga3_pos)
PGA3_UNARY(neg, pga3_neg)

/* ldual and rdual differ in even dimensions, so ~ is deliberately undefined. */
static PyObject *
pypga3_invert(PyObject *Py_UNUSED(self))
{
    PyErr_SetString(PyExc_TypeError,
                    "~ is undefined for pga3: ldual and rdual differ in "
                    "even dimensions; use ldual() or rdual()");
    return NULL;
}

static PyObject *
pypga3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue((PyTypeObject *)cls, pga3_zero());
}

static PyObject *
pypga3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue((PyTypeObject *)cls, pga3_one());
}

static PyObject *
pypga3_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue((PyTypeObject *)cls, pga3_bottom());
}

static PyObject *
pypga3_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue((PyTypeObject *)cls, pga3_top());
}

static PyObject *
pypga3_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_ldual(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_rdual(((PyPga3 *)self)->v));
}

static int
pypga3_bool(PyObject *self)
{
    pga3 v = ((PyPga3 *)self)->v;
    return vec16_any_nonzero(v);
}

#define PGA3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pypga3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_PGA3(op(((PyPga3 *)self)->v));                   \
    }

PGA3_METHOD(conjugate, pga3_conjugate)
PGA3_METHOD(exp, pga3_exp)
PGA3_METHOD(log, pga3_log)
PGA3_METHOD(sqrt, pga3_sqrt)

static PyObject *
pypga3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga3_norm(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_PGA3(pga3_inv(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga3_lmat(((PyPga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga3_rmat(((PyPga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga3_minpoly(((PyPga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    pga3_eivals(((PyPga3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pypga3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga3);
    pga3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga3 *)PyArray_DATA(arr);
    pga3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga3_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga3);
    pga3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga3 *)PyArray_DATA(arr);
    pga3_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_pga3);
    pga3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (pga3 *)PyArray_DATA(arr);
    pga3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pypga3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    pga3 v = ((PyPga3 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pypga3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_square(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_expm1(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_log1p(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_sin(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_cos(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_tan(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_asin(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_acos(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_atan(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_sinh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_cosh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_tanh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_asinh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_acosh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_atanh(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_re(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga3_FromValue(Py_TYPE(self), pga3_im(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_real(PyObject *self, void *closure)
{
    return pypga3_re(self, NULL);
}

static PyObject *
pypga3_imag(PyObject *self, void *closure)
{
    return pypga3_im(self, NULL);
}

static PyObject *pypga3_grades(PyObject *self, void *closure);

static PyGetSetDef pypga3_getset[] = {
    {"grade", pypga3_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pypga3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pypga3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pypga3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(pga3_is_real(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga3_to_re(((PyPga3 *)self)->v));
}

static PyObject *
pypga3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyPga3_FromValue((PyTypeObject *)cls, pga3_from_re(a));
}

static PyObject *
pypga3_lpow(PyObject *self, PyObject *other)
{
    pga3 b;
    int r;
    r = pga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga3_FromValue(Py_TYPE(self),
        pga3_lpow(((PyPga3 *)self)->v, b));
}

static PyObject *
pypga3_rpow(PyObject *self, PyObject *other)
{
    pga3 b;
    int r;
    r = pga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga3_FromValue(Py_TYPE(self),
        pga3_rpow(((PyPga3 *)self)->v, b));
}

static PyObject *
pypga3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    pga3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = pga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga3_FromValue(Py_TYPE(self),
        pga3_mix(((PyPga3 *)self)->v, b, t));
}

static PyObject *
pypga3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_pga3, 16);
}

static PyObject *
pypga3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_pga3, 16);
}

static PyObject *
pypga3_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_PGA3(pga3_ipow(((PyPga3 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pypga3_grade(PyObject *self, Py_ssize_t k)
{
    pga3 v = ((PyPga3 *)self)->v;
    switch (k) {
        case 0:
            return PyPga3_FromValue(Py_TYPE(self), pga3_grade0(v));
        case 1:
            return PyPga3_FromValue(Py_TYPE(self), pga3_grade1(v));
        case 2:
            return PyPga3_FromValue(Py_TYPE(self), pga3_grade2(v));
        case 3:
            return PyPga3_FromValue(Py_TYPE(self), pga3_grade3(v));
        case 4:
            return PyPga3_FromValue(Py_TYPE(self), pga3_grade4(v));
        default:
            return PyPga3_FromValue(Py_TYPE(self), pga3_zero());
    }
}

static PyObject *
pypga3_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pypga3_grade, 5);
}

static PyMethodDef pypga3_methods[] = {
    {"to_float", (PyCFunction)pypga3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pypga3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pypga3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pypga3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pypga3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pypga3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pypga3_square, METH_NOARGS, "square."},
    {"expm1", pypga3_expm1, METH_NOARGS, "expm1."},
    {"log1p", pypga3_log1p, METH_NOARGS, "log1p."},
    {"sin", pypga3_sin, METH_NOARGS, "sin."},
    {"cos", pypga3_cos, METH_NOARGS, "cos."},
    {"tan", pypga3_tan, METH_NOARGS, "tan."},
    {"asin", pypga3_asin, METH_NOARGS, "asin."},
    {"acos", pypga3_acos, METH_NOARGS, "acos."},
    {"atan", pypga3_atan, METH_NOARGS, "atan."},
    {"sinh", pypga3_sinh, METH_NOARGS, "sinh."},
    {"cosh", pypga3_cosh, METH_NOARGS, "cosh."},
    {"tanh", pypga3_tanh, METH_NOARGS, "tanh."},
    {"asinh", pypga3_asinh, METH_NOARGS, "asinh."},
    {"acosh", pypga3_acosh, METH_NOARGS, "acosh."},
    {"atanh", pypga3_atanh, METH_NOARGS, "atanh."},
    {"re", pypga3_re, METH_NOARGS, "re."},
    {"im", pypga3_im, METH_NOARGS, "im."},
    {"arcsin", pypga3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pypga3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pypga3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pypga3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pypga3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pypga3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pypga3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pypga3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pypga3_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pypga3_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pypga3_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pypga3_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"basis", pypga3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pypga3_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pypga3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pypga3_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pypga3_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pypga3_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pypga3_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pypga3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pypga3_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pypga3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pypga3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pypga3_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pypga3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pypga3_exp, METH_NOARGS, "exp."},
    {"log", pypga3_log, METH_NOARGS, "log."},
    {"sqrt", pypga3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pypga3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pypga3_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pypga3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyPga3 *)self)->v.data[i]);
}

static PyObject *
pypga3_subscript(PyObject *self, PyObject *key)
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
        i += 16;
    }
    return pypga3_item(self, i);
}

static int
pypga3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pypga3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pypga3_slots[] = {
    {Py_tp_dealloc, pypga3_dealloc},
    {Py_tp_traverse, pypga3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pypga3_repr},
    {Py_tp_str, pypga3_repr},
    {Py_tp_hash, pypga3_hash},
    {Py_tp_doc, (void *)"Projective geometric algebra in 3D with 16 float64 coefficients"},
    {Py_tp_richcompare, pypga3_richcompare},
    {Py_tp_methods, pypga3_methods},
    {Py_tp_getset, pypga3_getset},
    {Py_sq_length, pypga3_length},
    {Py_sq_item, pypga3_item},
    {Py_mp_subscript, pypga3_subscript},
    {Py_tp_new, pypga3_new},
    {Py_nb_add, pypga3_add},
    {Py_nb_subtract, pypga3_sub},
    {Py_nb_multiply, pypga3_mul},
    {Py_nb_true_divide, pypga3_div},
    {Py_nb_power, (void *)pypga3_pow},
    {Py_nb_positive, pypga3_pos},
    {Py_nb_negative, pypga3_neg},
    {Py_nb_or, pypga3_join},
    {Py_nb_and, pypga3_meet},
    {Py_nb_invert, pypga3_invert},
    {Py_nb_bool, pypga3_bool},
    {0, NULL}
};

static PyType_Spec pypga3_spec = {
    .name = "algpy.algebras.clifford.pga3",
    .basicsize = sizeof(PyPga3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pypga3_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npypga3_getitem(void *data, void *arr)
{
    pga3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_PGA3(v);
}

static int
npypga3_setitem(PyObject *item, void *data, void *arr)
{
    pga3 v;
    int r = pga3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to pga3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "pga3 needs exactly 16 components");
            return -1;
        }
        for (int i = 0; i < 16; i++) {
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
pypga3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(pga3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pypga3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pypga3_copyswapn(dst, sizeof(pga3), src, sizeof(pga3), 1,
                         swap, arr);
}

static npy_bool
pypga3_nonzero(void *data, void *arr)
{
    pga3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pypga3_arrfuncs;
typedef struct { char c; pga3 v; } pypga3_align_test;

static PyArray_DescrProto pypga3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(pga3), offsetof(pypga3_align_test, v), 0, 0, 0,
    &pypga3_arrfuncs,
};

#define PGA3_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pypga3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        pga3 *dst = (pga3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = pga3_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define PGA3_REAL_TYPES(X) \
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

PGA3_REAL_TYPES(PGA3_CAST_FROM)

static int
pypga3_register_casts(void)
{
#define REGISTER_PGA3_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_pga3,                    \
                                     pypga3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_pga3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    PGA3_REAL_TYPES(REGISTER_PGA3_CAST)
#undef REGISTER_PGA3_CAST
    return 0;
}

typedef pga3 (*pga3_unary_op)(pga3);
typedef pga3 (*pga3_binary_op)(pga3, pga3);
typedef int (*pga3_predicate)(pga3, pga3);

static void
loop_pga3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    pga3_unary_op op = (pga3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_pga3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    pga3_binary_op op = (pga3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga3 a, b, r;
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
loop_pga3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    pga3_predicate op = (pga3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pypga3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_pga3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_pga3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pypga3_register_ufuncs(PyObject *numpy)
{
    return
        pypga3_register_loop(numpy, "add", loop_pga3_binary, 3, 0, pga3_add) ||
        pypga3_register_loop(numpy, "subtract", loop_pga3_binary, 3, 0, pga3_sub) ||
        pypga3_register_loop(numpy, "multiply", loop_pga3_binary, 3, 0, pga3_mul) ||
        pypga3_register_loop(numpy, "divide", loop_pga3_binary, 3, 0, pga3_div) ||
        pypga3_register_loop(numpy, "power", loop_pga3_binary, 3, 0, pga3_rpow) ||
        pypga3_register_loop(numpy, "positive", loop_pga3_unary, 2, 0, pga3_pos) ||
        pypga3_register_loop(numpy, "negative", loop_pga3_unary, 2, 0, pga3_neg) ||
        pypga3_register_loop(numpy, "conjugate", loop_pga3_unary, 2, 0, pga3_conjugate) ||
        pypga3_register_loop(numpy, "reciprocal", loop_pga3_unary, 2, 0, pga3_inv) ||
        pypga3_register_loop(numpy, "square", loop_pga3_unary, 2, 0, pga3_square) ||
        pypga3_register_loop(numpy, "exp", loop_pga3_unary, 2, 0, pga3_exp) ||
        pypga3_register_loop(numpy, "expm1", loop_pga3_unary, 2, 0, pga3_expm1) ||
        pypga3_register_loop(numpy, "log", loop_pga3_unary, 2, 0, pga3_log) ||
        pypga3_register_loop(numpy, "log1p", loop_pga3_unary, 2, 0, pga3_log1p) ||
        pypga3_register_loop(numpy, "sqrt", loop_pga3_unary, 2, 0, pga3_sqrt) ||
        pypga3_register_loop(numpy, "sin", loop_pga3_unary, 2, 0, pga3_sin) ||
        pypga3_register_loop(numpy, "cos", loop_pga3_unary, 2, 0, pga3_cos) ||
        pypga3_register_loop(numpy, "tan", loop_pga3_unary, 2, 0, pga3_tan) ||
        pypga3_register_loop(numpy, "arcsin", loop_pga3_unary, 2, 0, pga3_asin) ||
        pypga3_register_loop(numpy, "arccos", loop_pga3_unary, 2, 0, pga3_acos) ||
        pypga3_register_loop(numpy, "arctan", loop_pga3_unary, 2, 0, pga3_atan) ||
        pypga3_register_loop(numpy, "sinh", loop_pga3_unary, 2, 0, pga3_sinh) ||
        pypga3_register_loop(numpy, "cosh", loop_pga3_unary, 2, 0, pga3_cosh) ||
        pypga3_register_loop(numpy, "tanh", loop_pga3_unary, 2, 0, pga3_tanh) ||
        pypga3_register_loop(numpy, "arcsinh", loop_pga3_unary, 2, 0, pga3_asinh) ||
        pypga3_register_loop(numpy, "arccosh", loop_pga3_unary, 2, 0, pga3_acosh) ||
        pypga3_register_loop(numpy, "arctanh", loop_pga3_unary, 2, 0, pga3_atanh) ||
        pypga3_register_loop(numpy, "equal", loop_pga3_predicate, 3, 1, pga3_equal) ||
        pypga3_register_loop(numpy, "bitwise_or", loop_pga3_binary, 3, 0, pga3_join) ||
        pypga3_register_loop(numpy, "bitwise_and", loop_pga3_binary, 3, 0, pga3_meet) ||
        pypga3_register_loop(numpy, "not_equal", loop_pga3_predicate, 3, 1, pga3_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pypga3_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pypga3_register_loop(ufuncs, "ldual", loop_pga3_unary, 2, 0, pga3_ldual) ||
        pypga3_register_loop(ufuncs, "rdual", loop_pga3_unary, 2, 0, pga3_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_pga3", "pga3 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__pga3(void)
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
    PyPga3_Type = PyType_FromModuleAndSpec(
        m, &pypga3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyPga3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pypga3_arrfuncs);
    pypga3_arrfuncs.getitem = npypga3_getitem;
    pypga3_arrfuncs.setitem = npypga3_setitem;
    pypga3_arrfuncs.copyswapn = pypga3_copyswapn;
    pypga3_arrfuncs.copyswap = pypga3_copyswap;
    pypga3_arrfuncs.nonzero = pypga3_nonzero;
    pypga3_descr_proto.typeobj = (PyTypeObject *)PyPga3_Type;
    Py_SET_TYPE((PyObject *)&pypga3_descr_proto, &PyArrayDescr_Type);
    npy_pga3 = PyArray_RegisterDataType(&pypga3_descr_proto);
    if (npy_pga3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_pga3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyPga3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyPga3_Type, 16) < 0) {
        goto fail;
    }
    if (pypga3_register_casts() < 0 ||
            pypga3_register_ufuncs(numpy) < 0 ||
            pypga3_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "pga3", PyPga3_Type) < 0) {
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
