/*
 * cga2 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in cga2.h.
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

#include "cga2.h"

typedef struct {
    PyObject_HEAD
    cga2 v;
} PyCga2;

static PyObject *PyCga2_Type = NULL;
static int npy_cga2 = -1;

static PyObject *
PyCga2_FromValue(PyTypeObject *type, cga2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCga2 *p = (PyCga2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_CGA2(v) \
    PyCga2_FromValue((PyTypeObject *)PyCga2_Type, (v))

static int
cga2_from_object(PyObject *o, cga2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyCga2_Type)) {
        *out = ((PyCga2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cga2_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycga2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PyCga2_FromValue(type, cga2_make(d));
}

static PyObject *
pycga2_repr(PyObject *self)
{
    cga2 v = ((PyCga2 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cga2(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycga2_hash(PyObject *self)
{
    cga2 v = ((PyCga2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CGA2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cga2_is_real(v)) {
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
pycga2_richcompare(PyObject *a, PyObject *b, int op)
{
    cga2 x, y;
    int r;
    if ((r = cga2_from_object(a, &x)) <= 0 ||
            (r = cga2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(cga2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!cga2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define CGA2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycga2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        cga2 x, y;                                                       \
        int r;                                                               \
        if ((r = cga2_from_object(a, &x)) <= 0 ||                        \
                (r = cga2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_CGA2(expr(x, y));                                     \
    }

CGA2_BINARY(add, cga2_add)
CGA2_BINARY(sub, cga2_sub)
CGA2_BINARY(mul, cga2_mul)
CGA2_BINARY(div, cga2_div)
CGA2_BINARY(join, cga2_join)
CGA2_BINARY(meet, cga2_meet)

static PyObject *
pycga2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cga2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cga2 does not support modular power");
        return NULL;
    }
    if ((r = cga2_from_object(a, &x)) <= 0 ||
            (r = cga2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_CGA2(cga2_rpow(x, y));
}

#define CGA2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycga2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_CGA2(expr(((PyCga2 *)self)->v));                  \
    }

CGA2_UNARY(pos, cga2_pos)
CGA2_UNARY(neg, cga2_neg)

/* ldual and rdual differ in even dimensions, so ~ is deliberately undefined. */
static PyObject *
pycga2_invert(PyObject *Py_UNUSED(self))
{
    PyErr_SetString(PyExc_TypeError,
                    "~ is undefined for cga2: ldual and rdual differ in "
                    "even dimensions; use ldual() or rdual()");
    return NULL;
}

static PyObject *
pycga2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue((PyTypeObject *)cls, cga2_zero());
}

static PyObject *
pycga2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue((PyTypeObject *)cls, cga2_one());
}

static PyObject *
pycga2_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue((PyTypeObject *)cls, cga2_bottom());
}

static PyObject *
pycga2_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue((PyTypeObject *)cls, cga2_top());
}

static PyObject *
pycga2_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_ldual(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_rdual(((PyCga2 *)self)->v));
}

static int
pycga2_bool(PyObject *self)
{
    cga2 v = ((PyCga2 *)self)->v;
    return vec16_any_nonzero(v);
}

#define CGA2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycga2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_CGA2(op(((PyCga2 *)self)->v));                   \
    }

CGA2_METHOD(conjugate, cga2_conjugate)
CGA2_METHOD(exp, cga2_exp)
CGA2_METHOD(log, cga2_log)
CGA2_METHOD(sqrt, cga2_sqrt)

static PyObject *
pycga2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = cga2_norm(((PyCga2 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycga2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_CGA2(cga2_inv(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga2_lmat(((PyCga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga2_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga2_rmat(((PyCga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga2_minpoly(((PyCga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    cga2_eivals(((PyCga2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycga2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga2);
    cga2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga2 *)PyArray_DATA(arr);
    cga2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga2_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga2);
    cga2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga2 *)PyArray_DATA(arr);
    cga2_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {16, 16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_cga2);
    cga2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (cga2 *)PyArray_DATA(arr);
    cga2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycga2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    cga2 v = ((PyCga2 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pycga2_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_square(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_expm1(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_log1p(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_sin(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_cos(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_tan(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_asin(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_acos(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_atan(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_sinh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_cosh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_tanh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_asinh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_acosh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_atanh(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_re(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga2_FromValue(Py_TYPE(self), cga2_im(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_real(PyObject *self, void *closure)
{
    return pycga2_re(self, NULL);
}

static PyObject *
pycga2_imag(PyObject *self, void *closure)
{
    return pycga2_im(self, NULL);
}

static PyObject *pycga2_grades(PyObject *self, void *closure);

static PyGetSetDef pycga2_getset[] = {
    {"grade", pycga2_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycga2_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycga2_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycga2_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cga2_is_real(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cga2_to_re(((PyCga2 *)self)->v));
}

static PyObject *
pycga2_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyCga2_FromValue((PyTypeObject *)cls, cga2_from_re(a));
}

static PyObject *
pycga2_lpow(PyObject *self, PyObject *other)
{
    cga2 b;
    int r;
    r = cga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga2_FromValue(Py_TYPE(self),
        cga2_lpow(((PyCga2 *)self)->v, b));
}

static PyObject *
pycga2_rpow(PyObject *self, PyObject *other)
{
    cga2 b;
    int r;
    r = cga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga2_FromValue(Py_TYPE(self),
        cga2_rpow(((PyCga2 *)self)->v, b));
}

static PyObject *
pycga2_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cga2 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga2_FromValue(Py_TYPE(self),
        cga2_mix(((PyCga2 *)self)->v, b, t));
}

static PyObject *
pycga2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_cga2, 16);
}

static PyObject *
pycga2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_cga2, 16);
}

static PyObject *
pycga2_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_CGA2(cga2_ipow(((PyCga2 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycga2_grade(PyObject *self, Py_ssize_t k)
{
    cga2 v = ((PyCga2 *)self)->v;
    switch (k) {
        case 0:
            return PyCga2_FromValue(Py_TYPE(self), cga2_grade0(v));
        case 1:
            return PyCga2_FromValue(Py_TYPE(self), cga2_grade1(v));
        case 2:
            return PyCga2_FromValue(Py_TYPE(self), cga2_grade2(v));
        case 3:
            return PyCga2_FromValue(Py_TYPE(self), cga2_grade3(v));
        case 4:
            return PyCga2_FromValue(Py_TYPE(self), cga2_grade4(v));
        default:
            return PyCga2_FromValue(Py_TYPE(self), cga2_zero());
    }
}

static PyObject *
pycga2_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycga2_grade, 5);
}

static PyMethodDef pycga2_methods[] = {
    {"to_float", (PyCFunction)pycga2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycga2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycga2_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycga2_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycga2_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycga2_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycga2_square, METH_NOARGS, "square."},
    {"expm1", pycga2_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycga2_log1p, METH_NOARGS, "log1p."},
    {"sin", pycga2_sin, METH_NOARGS, "sin."},
    {"cos", pycga2_cos, METH_NOARGS, "cos."},
    {"tan", pycga2_tan, METH_NOARGS, "tan."},
    {"asin", pycga2_asin, METH_NOARGS, "asin."},
    {"acos", pycga2_acos, METH_NOARGS, "acos."},
    {"atan", pycga2_atan, METH_NOARGS, "atan."},
    {"sinh", pycga2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycga2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycga2_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycga2_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycga2_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycga2_atanh, METH_NOARGS, "atanh."},
    {"re", pycga2_re, METH_NOARGS, "re."},
    {"im", pycga2_im, METH_NOARGS, "im."},
    {"arcsin", pycga2_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycga2_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycga2_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycga2_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycga2_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycga2_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycga2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycga2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pycga2_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pycga2_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pycga2_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pycga2_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"basis", pycga2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycga2_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pycga2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycga2_norm, METH_NOARGS, "Norm as a Python complex."},
    {"lpow", pycga2_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycga2_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycga2_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycga2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycga2_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pycga2_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycga2_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycga2_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"eivals", pycga2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycga2_exp, METH_NOARGS, "exp."},
    {"log", pycga2_log, METH_NOARGS, "log."},
    {"sqrt", pycga2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycga2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycga2_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pycga2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCga2 *)self)->v.data[i]);
}

static PyObject *
pycga2_subscript(PyObject *self, PyObject *key)
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
    return pycga2_item(self, i);
}

static int
pycga2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycga2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycga2_slots[] = {
    {Py_tp_dealloc, pycga2_dealloc},
    {Py_tp_traverse, pycga2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycga2_repr},
    {Py_tp_str, pycga2_repr},
    {Py_tp_hash, pycga2_hash},
    {Py_tp_doc, (void *)"Conformal geometric algebra in 2D with 16 float64 coefficients"},
    {Py_tp_richcompare, pycga2_richcompare},
    {Py_tp_methods, pycga2_methods},
    {Py_tp_getset, pycga2_getset},
    {Py_sq_length, pycga2_length},
    {Py_sq_item, pycga2_item},
    {Py_mp_subscript, pycga2_subscript},
    {Py_tp_new, pycga2_new},
    {Py_nb_add, pycga2_add},
    {Py_nb_subtract, pycga2_sub},
    {Py_nb_multiply, pycga2_mul},
    {Py_nb_true_divide, pycga2_div},
    {Py_nb_power, (void *)pycga2_pow},
    {Py_nb_positive, pycga2_pos},
    {Py_nb_negative, pycga2_neg},
    {Py_nb_or, pycga2_join},
    {Py_nb_and, pycga2_meet},
    {Py_nb_invert, pycga2_invert},
    {Py_nb_bool, pycga2_bool},
    {0, NULL}
};

static PyType_Spec pycga2_spec = {
    .name = "algpy.algebras.clifford.cga2",
    .basicsize = sizeof(PyCga2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycga2_slots,
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
npycga2_getitem(void *data, void *arr)
{
    cga2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_CGA2(v);
}

static int
npycga2_setitem(PyObject *item, void *data, void *arr)
{
    cga2 v;
    int r = cga2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cga2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cga2 needs exactly 16 components");
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
pycga2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cga2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycga2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycga2_copyswapn(dst, sizeof(cga2), src, sizeof(cga2), 1,
                         swap, arr);
}

static npy_bool
pycga2_nonzero(void *data, void *arr)
{
    cga2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycga2_arrfuncs;
typedef struct { char c; cga2 v; } pycga2_align_test;

static PyArray_DescrProto pycga2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(cga2), offsetof(pycga2_align_test, v), 0, 0, 0,
    &pycga2_arrfuncs,
};

#define CGA2_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycga2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        cga2 *dst = (cga2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = cga2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define CGA2_REAL_TYPES(X) \
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

CGA2_REAL_TYPES(CGA2_CAST_FROM)

static int
pycga2_register_casts(void)
{
#define REGISTER_CGA2_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_cga2,                    \
                                     pycga2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_cga2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    CGA2_REAL_TYPES(REGISTER_CGA2_CAST)
#undef REGISTER_CGA2_CAST
    return 0;
}

typedef cga2 (*cga2_unary_op)(cga2);
typedef cga2 (*cga2_binary_op)(cga2, cga2);
typedef int (*cga2_predicate)(cga2, cga2);

static void
loop_cga2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    cga2_unary_op op = (cga2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_cga2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    cga2_binary_op op = (cga2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga2 a, b, r;
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
loop_cga2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    cga2_predicate op = (cga2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycga2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_cga2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_cga2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycga2_register_ufuncs(PyObject *numpy)
{
    return
        pycga2_register_loop(numpy, "add", loop_cga2_binary, 3, 0, cga2_add) ||
        pycga2_register_loop(numpy, "subtract", loop_cga2_binary, 3, 0, cga2_sub) ||
        pycga2_register_loop(numpy, "multiply", loop_cga2_binary, 3, 0, cga2_mul) ||
        pycga2_register_loop(numpy, "divide", loop_cga2_binary, 3, 0, cga2_div) ||
        pycga2_register_loop(numpy, "power", loop_cga2_binary, 3, 0, cga2_rpow) ||
        pycga2_register_loop(numpy, "positive", loop_cga2_unary, 2, 0, cga2_pos) ||
        pycga2_register_loop(numpy, "negative", loop_cga2_unary, 2, 0, cga2_neg) ||
        pycga2_register_loop(numpy, "conjugate", loop_cga2_unary, 2, 0, cga2_conjugate) ||
        pycga2_register_loop(numpy, "reciprocal", loop_cga2_unary, 2, 0, cga2_inv) ||
        pycga2_register_loop(numpy, "square", loop_cga2_unary, 2, 0, cga2_square) ||
        pycga2_register_loop(numpy, "exp", loop_cga2_unary, 2, 0, cga2_exp) ||
        pycga2_register_loop(numpy, "expm1", loop_cga2_unary, 2, 0, cga2_expm1) ||
        pycga2_register_loop(numpy, "log", loop_cga2_unary, 2, 0, cga2_log) ||
        pycga2_register_loop(numpy, "log1p", loop_cga2_unary, 2, 0, cga2_log1p) ||
        pycga2_register_loop(numpy, "sqrt", loop_cga2_unary, 2, 0, cga2_sqrt) ||
        pycga2_register_loop(numpy, "sin", loop_cga2_unary, 2, 0, cga2_sin) ||
        pycga2_register_loop(numpy, "cos", loop_cga2_unary, 2, 0, cga2_cos) ||
        pycga2_register_loop(numpy, "tan", loop_cga2_unary, 2, 0, cga2_tan) ||
        pycga2_register_loop(numpy, "arcsin", loop_cga2_unary, 2, 0, cga2_asin) ||
        pycga2_register_loop(numpy, "arccos", loop_cga2_unary, 2, 0, cga2_acos) ||
        pycga2_register_loop(numpy, "arctan", loop_cga2_unary, 2, 0, cga2_atan) ||
        pycga2_register_loop(numpy, "sinh", loop_cga2_unary, 2, 0, cga2_sinh) ||
        pycga2_register_loop(numpy, "cosh", loop_cga2_unary, 2, 0, cga2_cosh) ||
        pycga2_register_loop(numpy, "tanh", loop_cga2_unary, 2, 0, cga2_tanh) ||
        pycga2_register_loop(numpy, "arcsinh", loop_cga2_unary, 2, 0, cga2_asinh) ||
        pycga2_register_loop(numpy, "arccosh", loop_cga2_unary, 2, 0, cga2_acosh) ||
        pycga2_register_loop(numpy, "arctanh", loop_cga2_unary, 2, 0, cga2_atanh) ||
        pycga2_register_loop(numpy, "equal", loop_cga2_predicate, 3, 1, cga2_equal) ||
        pycga2_register_loop(numpy, "bitwise_or", loop_cga2_binary, 3, 0, cga2_join) ||
        pycga2_register_loop(numpy, "bitwise_and", loop_cga2_binary, 3, 0, cga2_meet) ||
        pycga2_register_loop(numpy, "not_equal", loop_cga2_predicate, 3, 1, cga2_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pycga2_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pycga2_register_loop(ufuncs, "ldual", loop_cga2_unary, 2, 0, cga2_ldual) ||
        pycga2_register_loop(ufuncs, "rdual", loop_cga2_unary, 2, 0, cga2_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cga2", "cga2 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cga2(void)
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
    PyCga2_Type = PyType_FromModuleAndSpec(
        m, &pycga2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyCga2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycga2_arrfuncs);
    pycga2_arrfuncs.getitem = npycga2_getitem;
    pycga2_arrfuncs.setitem = npycga2_setitem;
    pycga2_arrfuncs.copyswapn = pycga2_copyswapn;
    pycga2_arrfuncs.copyswap = pycga2_copyswap;
    pycga2_arrfuncs.nonzero = pycga2_nonzero;
    pycga2_descr_proto.typeobj = (PyTypeObject *)PyCga2_Type;
    Py_SET_TYPE((PyObject *)&pycga2_descr_proto, &PyArrayDescr_Type);
    npy_cga2 = PyArray_RegisterDataType(&pycga2_descr_proto);
    if (npy_cga2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_cga2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyCga2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyCga2_Type, 16) < 0) {
        goto fail;
    }
    if (pycga2_register_casts() < 0 ||
            pycga2_register_ufuncs(numpy) < 0 ||
            pycga2_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "cga2", PyCga2_Type) < 0) {
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
