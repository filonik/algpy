/*
 * Cayley-Dickson numbers x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 +
 * The 16 coefficients of e0 through e3, with e_i*e_i = k_i, exposed to Python
 * as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `cayleydickson4(k0, k1, k2, k3)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `cayleydickson4.h`.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "../../../_conversions.h"

#include "cayleydickson4.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec4 k;
} Info;

typedef struct {
    PyObject_HEAD
    cayleydickson4 v;
} PyCD4;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.cayleydickson4.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cd4_from_value(PyTypeObject *type, cayleydickson4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCD4 *p = (PyCD4 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a cayleydickson4 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cd4_from_object(Info *info, PyObject *o, cayleydickson4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCD4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cayleydickson4_make(
            d, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycd4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12",
                             "x012", "x3", "x03", "x13", "x013", "x23",
                             "x023", "x123", "x0123", NULL};
    double c[16] = {0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
                                     c, c + 1, c + 2, c + 3, c + 4, c + 5,
                                     c + 6, c + 7, c + 8, c + 9, c + 10,
                                     c + 11, c + 12, c + 13, c + 14, c + 15)) {
        return NULL;
    }
    return cd4_from_value(type, cayleydickson4_make(
        c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
        c[8], c[9], c[10], c[11], c[12], c[13], c[14], c[15]));
}

static PyObject *
pycd4_repr(PyObject *self)
{
    cayleydickson4 v = ((PyCD4 *)self)->v;
    char *s[16] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cayleydickson4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7],
                                  s[8], s[9], s[10], s[11], s[12], s[13], s[14],
                                  s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycd4_hash(PyObject *self)
{
    cayleydickson4 v = ((PyCD4 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CAYLEYDICKSON4_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cayleydickson4_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddd)",
                            v.data[0], v.data[1], v.data[2], v.data[3],
                            v.data[4], v.data[5], v.data[6], v.data[7],
                            v.data[8], v.data[9], v.data[10], v.data[11],
                            v.data[12], v.data[13], v.data[14], v.data[15]);
    }
    if (obj == NULL) {
        return -1;
    }
    h = PyObject_Hash(obj);
    Py_DECREF(obj);
    return h;
}

/* Resolve the type info of whichever operand is a cayleydickson4. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyTypeObject *type;
    /* A subclass may override __new__; recognize its algebra base type. */
    for (type = Py_TYPE(a); type != NULL; type = type->tp_base) {
        if (PyType_GetSlot(type, Py_tp_new) == (void *)pycd4_new) {
            return info_of_type(Py_TYPE(a));
        }
    }
    return info_of_type(Py_TYPE(b));
}

static PyObject *
pycd4_richcompare(PyObject *a, PyObject *b, int op)
{
    cayleydickson4 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cd4_from_object(info, a, &x)) <= 0 ||
            (r = cd4_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(cayleydickson4_equal(x, y) == (op == Py_EQ));
}

#define CD4_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycd4_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        cayleydickson4 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cd4_from_object(info, a, &x)) <= 0 ||                     \
                (r = cd4_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cd4_from_value((PyTypeObject *)info->type, expr);           \
    }

CD4_BINARY(add, cayleydickson4_add(x, y))
CD4_BINARY(sub, cayleydickson4_sub(x, y))
CD4_BINARY(mul, cayleydickson4_mul(info->k, x, y))
CD4_BINARY(div, cayleydickson4_div(info->k, x, y))

static PyObject *
pycd4_pos(PyObject *self)
{
    return cd4_from_value(Py_TYPE(self), cayleydickson4_pos(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_neg(PyObject *self)
{
    return cd4_from_value(Py_TYPE(self), cayleydickson4_neg(((PyCD4 *)self)->v));
}

static int
pycd4_bool(PyObject *self)
{
    return vec16_any_nonzero(((PyCD4 *)self)->v);
}

static PyObject *
pycd4_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd4_from_value(Py_TYPE(self),
                          cayleydickson4_conjugate(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd4_from_value((PyTypeObject *)cls, cayleydickson4_zero());
}

static PyObject *
pycd4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd4_from_value((PyTypeObject *)cls, cayleydickson4_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycd4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {16, 16};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_lmat(info->k, ((PyCD4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {16, 16};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_rmat(info->k, ((PyCD4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {3};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_minpoly(info->k, ((PyCD4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {2};
    complex ev[2];
    npy_cdouble *d;
    PyObject *arr;
    int i;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_eivals(info->k, ((PyCD4 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pycd4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    cayleydickson4 v = ((PyCD4 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(dddddddddddddddddddd))", rebuild_func, info->k.data[0],
                         info->k.data[1], info->k.data[2], info->k.data[3], v.data[0],
                         v.data[1], v.data[2], v.data[3], v.data[4],
                         v.data[5], v.data[6], v.data[7], v.data[8],
                         v.data[9], v.data[10], v.data[11], v.data[12],
                         v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pycd4_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CAYLEYDICKSON4_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_basis((cayleydickson4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd4_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CAYLEYDICKSON4_N, ALGPY_CAYLEYDICKSON4_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson4_coeffs(
        info->k, (cayleydickson4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = cayleydickson4_norm(info->k, ((PyCD4 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycd4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_inv(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_square(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_exp(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_expm1(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_log(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_log1p(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_sqrt(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_sin(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_cos(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_tan(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_asin(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_acos(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_atan(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_sinh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_cosh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_tanh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_asinh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_acosh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_atanh(info->k, ((PyCD4 *)self)->v));
}

static PyObject *
pycd4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd4_from_value(Py_TYPE(self), cayleydickson4_re(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd4_from_value(Py_TYPE(self), cayleydickson4_im(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_real(PyObject *self, void *closure)
{
    return pycd4_re(self, NULL);
}

static PyObject *
pycd4_imag(PyObject *self, void *closure)
{
    return pycd4_im(self, NULL);
}

static PyGetSetDef pycd4_getset[] = {
    {"real", pycd4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycd4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycd4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cayleydickson4_is_real(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cayleydickson4_to_re(((PyCD4 *)self)->v));
}

static PyObject *
pycd4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cd4_from_value((PyTypeObject *)cls, cayleydickson4_from_re(a));
}

static PyObject *
pycd4_lpow(PyObject *self, PyObject *other)
{
    cayleydickson4 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd4_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self),
        cayleydickson4_lpow(info->k, ((PyCD4 *)self)->v, b));
}

static PyObject *
pycd4_rpow(PyObject *self, PyObject *other)
{
    cayleydickson4 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd4_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self),
        cayleydickson4_rpow(info->k, ((PyCD4 *)self)->v, b));
}

static PyObject *
pycd4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cayleydickson4 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cd4_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self),
        cayleydickson4_mix(info->k, ((PyCD4 *)self)->v, b, t));
}

static PyObject *
pycd4_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 16);
}

static PyObject *
pycd4_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 16);
}

static PyObject *
pycd4_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd4_from_value(Py_TYPE(self), cayleydickson4_ipow(info->k, ((PyCD4 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pycd4_methods[] = {
    {"to_float", (PyCFunction)pycd4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycd4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycd4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycd4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycd4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycd4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pycd4_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycd4_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycd4_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycd4_inv, METH_NOARGS, "inv."},
    {"square", pycd4_square, METH_NOARGS, "square."},
    {"exp", pycd4_exp, METH_NOARGS, "exp."},
    {"expm1", pycd4_expm1, METH_NOARGS, "expm1."},
    {"log", pycd4_log, METH_NOARGS, "log."},
    {"log1p", pycd4_log1p, METH_NOARGS, "log1p."},
    {"sqrt", pycd4_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pycd4_sin, METH_NOARGS, "sin."},
    {"cos", pycd4_cos, METH_NOARGS, "cos."},
    {"tan", pycd4_tan, METH_NOARGS, "tan."},
    {"asin", pycd4_asin, METH_NOARGS, "asin."},
    {"acos", pycd4_acos, METH_NOARGS, "acos."},
    {"atan", pycd4_atan, METH_NOARGS, "atan."},
    {"sinh", pycd4_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycd4_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycd4_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycd4_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycd4_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycd4_atanh, METH_NOARGS, "atanh."},
    {"re", pycd4_re, METH_NOARGS, "re."},
    {"im", pycd4_im, METH_NOARGS, "im."},
    {"arcsin", pycd4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycd4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycd4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycd4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycd4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycd4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycd4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycd4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycd4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycd4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycd4_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycd4_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pycd4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycd4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycd4_minpoly, METH_NOARGS, "Coefficients of the generic minimal polynomial."},
    {"eivals", pycd4_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycd4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycd4_length(PyObject *self)
{
    return 16;
}

static PyObject *
pycd4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCD4 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycd4_subscript(PyObject *self, PyObject *key)
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
    return pycd4_item(self, i);
}

static int
pycd4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycd4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycd4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cayleydickson4 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cayleydickson4 does not support modular power");
        return NULL;
    }
    if ((r = cd4_from_object(info, a, &x)) <= 0 ||
            (r = cd4_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cd4_from_value((PyTypeObject *)info->type, cayleydickson4_rpow(info->k, x, y));
}

static PyType_Slot pycd4_slots[] = {
    {Py_tp_dealloc, pycd4_dealloc},
    {Py_tp_traverse, pycd4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycd4_repr},
    {Py_tp_str, pycd4_repr},
    {Py_tp_hash, pycd4_hash},
    {Py_tp_doc, (void *)"Cayley-Dickson number with float64 components"},
    {Py_tp_richcompare, pycd4_richcompare},
    {Py_tp_methods, pycd4_methods},
    {Py_tp_getset, pycd4_getset},
    {Py_sq_length, pycd4_length},
    {Py_sq_item, pycd4_item},
    {Py_mp_subscript, pycd4_subscript},
    {Py_tp_new, pycd4_new},
    {Py_nb_add, pycd4_add},
    {Py_nb_subtract, pycd4_sub},
    {Py_nb_multiply, pycd4_mul},
    {Py_nb_true_divide, pycd4_div},
    {Py_nb_power, (void *)pycd4_pow},
    {Py_nb_positive, pycd4_pos},
    {Py_nb_negative, pycd4_neg},
    {Py_nb_bool, pycd4_bool},
    {0, NULL}
};

static PyType_Spec pycd4_spec = {
    .name = "algpy.algebras.cayleydickson4",
    .basicsize = sizeof(PyCD4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycd4_slots,
};

/* NumPy dtype support */

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
npycd4_getitem(void *data, void *arr)
{
    cayleydickson4 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return cd4_from_value(descr->typeobj, v);
}

static int
npycd4_setitem(PyObject *item, void *data, void *arr)
{
    cayleydickson4 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cd4_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cayleydickson4");
        double d[16];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cayleydickson4 needs exactly 16 components");
            return -1;
        }
        for (i = 0; i < 16; i++) {
            d[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (d[i] == -1.0 && PyErr_Occurred()) {
                Py_DECREF(seq);
                return -1;
            }
        }
        Py_DECREF(seq);
        v = cayleydickson4_make(
            d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7],
            d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
npycd4_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cayleydickson4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycd4_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(cayleydickson4));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycd4_nonzero(void *data, void *arr)
{
    cayleydickson4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycd4_arrfuncs;

typedef struct { char c; cayleydickson4 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CD4_CAST_FROM(name, type, convert)                                          \
    static void                                                            \
    npycast_##name##_to_cd4(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        cayleydickson4 *dst = (cayleydickson4 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = cayleydickson4_make(                                 \
                convert(src[i]), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); \
        }                                                                  \
    }

#define CD4_REAL_TYPES(X) \
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

CD4_REAL_TYPES(CD4_CAST_FROM)

static int
register_casts(int typenum)
{
#define REGISTER_CAST(name, type, convert)                                          \
    {                                                                      \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);           \
        if (from == NULL) {                                                \
            return -1;                                                     \
        }                                                                  \
        if (PyArray_RegisterCastFunc(from, typenum,                        \
                npycast_##name##_to_cd4) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CD4_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec4 *k;
    cayleydickson4 (*unary)(cayleydickson4);
    cayleydickson4 (*unary_k)(vec4, cayleydickson4);
    cayleydickson4 (*binary)(vec4, cayleydickson4, cayleydickson4);
    int (*pred)(cayleydickson4, cayleydickson4);
} Op;

static cayleydickson4 op_add(vec4 k, cayleydickson4 a, cayleydickson4 b)
{ return cayleydickson4_add(a, b); }
static cayleydickson4 op_sub(vec4 k, cayleydickson4 a, cayleydickson4 b)
{ return cayleydickson4_sub(a, b); }
static int op_equal(cayleydickson4 a, cayleydickson4 b)
{ return cayleydickson4_equal(a, b); }
static int op_not_equal(cayleydickson4 a, cayleydickson4 b)
{ return !cayleydickson4_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        cayleydickson4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op->unary(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_unary_k(char **args, npy_intp const *dims, npy_intp const *steps, void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cayleydickson4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op->unary_k(*op->k, a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_binary(char **args, npy_intp const *dims, npy_intp const *steps,
            void *data)
{
    Op *op = (Op *)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        cayleydickson4 a, b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = op->binary(*op->k, a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_predicate(char **args, npy_intp const *dims, npy_intp const *steps,
               void *data)
{
    Op *op = (Op *)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        cayleydickson4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op->pred(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
register_loop(Info *info, PyObject *numpy, const char *name,
              PyUFuncGenericFunction loop, int nargs, int out_bool, Op op)
{
    int types[3];
    int i, res;
    Op *data;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    data = (Op *)malloc(sizeof(Op));
    if (data == NULL) {
        Py_DECREF(ufunc);
        PyErr_NoMemory();
        return -1;
    }
    *data = op;
    data->k = &info->k;  /* read at call time: K is mutable */
    for (i = 0; i < nargs; i++) {
        types[i] = info->typenum;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    res = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, info->typenum,
                                      loop, types, data);
    Py_DECREF(ufunc);
    return res;
}

static int
register_ufuncs(Info *info, PyObject *numpy)
{
    Op none = {NULL, NULL, NULL, NULL, NULL};
    Op o;
#define BIN(name, f) (o = none, o.binary = f, register_loop(info, numpy, name, loop_binary, 3, 0, o))
#define UN(name, f) (o = none, o.unary = f, register_loop(info, numpy, name, loop_unary, 2, 0, o))
#define UNK(name, f) (o = none, o.unary_k = f, register_loop(info, numpy, name, loop_unary_k, 2, 0, o))
#define PRED(name, f) (o = none, o.pred = f, register_loop(info, numpy, name, loop_predicate, 3, 1, o))
    return
        BIN("add", op_add) ||
        BIN("subtract", op_sub) ||
        BIN("multiply", cayleydickson4_mul) ||
        BIN("divide", cayleydickson4_div) ||
        BIN("power", cayleydickson4_rpow) ||
        UN("positive", cayleydickson4_pos) ||
        UN("negative", cayleydickson4_neg) ||
        UN("conjugate", cayleydickson4_conjugate) ||
        UNK("reciprocal", cayleydickson4_inv) ||
        UNK("square", cayleydickson4_square) ||
        UNK("exp", cayleydickson4_exp) ||
        UNK("expm1", cayleydickson4_expm1) ||
        UNK("log", cayleydickson4_log) ||
        UNK("log1p", cayleydickson4_log1p) ||
        UNK("sqrt", cayleydickson4_sqrt) ||
        UNK("sin", cayleydickson4_sin) ||
        UNK("cos", cayleydickson4_cos) ||
        UNK("tan", cayleydickson4_tan) ||
        UNK("arcsin", cayleydickson4_asin) ||
        UNK("arccos", cayleydickson4_acos) ||
        UNK("arctan", cayleydickson4_atan) ||
        UNK("sinh", cayleydickson4_sinh) ||
        UNK("cosh", cayleydickson4_cosh) ||
        UNK("tanh", cayleydickson4_tanh) ||
        UNK("arcsinh", cayleydickson4_asinh) ||
        UNK("arccosh", cayleydickson4_acosh) ||
        UNK("arctanh", cayleydickson4_atanh) ||
        PRED("equal", op_equal) ||
        PRED("not_equal", op_not_equal)
        ? -1 : 0;
#undef BIN
#undef UN
#undef UNK
#undef PRED
}

/*
 * Metatype giving every scalar type a mutable parameter sequence `K`.
 * Info is created before the type becomes visible, so `info_of_type` works.
 */

static PyObject *
meta_get_K(PyObject *self, void *closure)
{
    Info *info = info_of_type((PyTypeObject *)self);
    PyObject *list;
    int i;
    if (info == NULL) {
        return NULL;
    }
    list = PyList_New(4);
    for (i = 0; list != NULL && i < 4; i++) {
        PyObject *v = PyFloat_FromDouble(info->k.data[i]);
        if (v == NULL) {
            Py_CLEAR(list);
            break;
        }
        PyList_SET_ITEM(list, i, v);
    }
    return list;
}

static int
meta_set_K(PyObject *self, PyObject *value, void *closure)
{
    Info *info = info_of_type((PyTypeObject *)self);
    PyObject *seq;
    vec4 k;
    int i;
    if (info == NULL) {
        return -1;
    }
    if (value == NULL) {
        PyErr_SetString(PyExc_AttributeError, "cannot delete K");
        return -1;
    }
    seq = PySequence_Fast(value, "K must be a sequence of 4 numbers");
    if (seq == NULL) {
        return -1;
    }
    if (PySequence_Fast_GET_SIZE(seq) != 4) {
        PyErr_SetString(PyExc_ValueError, "K must have exactly 4 element(s)");
        Py_DECREF(seq);
        return -1;
    }
    for (i = 0; i < 4; i++) {
        k.data[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        if (k.data[i] == -1.0 && PyErr_Occurred()) {
            Py_DECREF(seq);
            return -1;
        }
    }
    Py_DECREF(seq);
    info->k = k;
    return 0;
}

static PyGetSetDef meta_getset[] = {
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 4 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.cayleydickson4.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = meta_getset,
};

/* Type factory */

static PyObject *
create_type(double k0, double k1, double k2, double k3)
{
    Info *info = NULL;
    PyObject *type = NULL, *cap = NULL, *numpy = NULL, *descr_obj = NULL;
    PyArray_DescrProto *proto = NULL;
    PyArray_Descr *descr;

    info = (Info *)calloc(1, sizeof(Info));
    proto = (PyArray_DescrProto *)calloc(1, sizeof(PyArray_DescrProto));
    if (info == NULL || proto == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    info->k.data[0] = k0;
    info->k.data[1] = k1;
    info->k.data[2] = k2;
    info->k.data[3] = k3;

    type = PyType_FromSpecWithBases(&pycd4_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.cayleydickson4.info", NULL);
    if (cap == NULL ||
            PyObject_SetAttrString(type, "_info", cap) < 0) {
        goto fail;
    }
    Py_CLEAR(cap);

    proto->typeobj = (PyTypeObject *)type;
    proto->kind = 'V';
    proto->type = 'q';
    proto->byteorder = '=';
    proto->flags = NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM;
    proto->elsize = sizeof(cayleydickson4);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycd4_arrfuncs;
    Py_SET_TYPE((PyObject *)proto, &PyArrayDescr_Type);

    info->typenum = PyArray_RegisterDataType(proto);
    if (info->typenum < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(info->typenum);
    if (descr == NULL) {
        goto fail;
    }
    descr_obj = (PyObject *)descr;
    if (PyObject_SetAttrString(type, "dtype", descr_obj) < 0) {
        Py_DECREF(descr_obj);
        goto fail;
    }
    Py_DECREF(descr_obj);

    if (algpy_set_dimension(type, 16) < 0) {
        goto fail;
    }

    numpy = PyImport_ImportModule("numpy");
    if (numpy == NULL) {
        goto fail;
    }
    if (register_casts(info->typenum) < 0 ||
            register_ufuncs(info, numpy) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return type;

fail:
    /* A registered dtype cannot be unregistered; only free if not registered. */
    Py_XDECREF(cap);
    Py_XDECREF(numpy);
    return NULL;
}

static PyObject *
cd4_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3;
    if (!PyArg_ParseTuple(args, "dddd", &k0, &k1, &k2, &k3)) {
        return NULL;
    }
    return create_type(k0, k1, k2, k3);
}

static PyObject *
cd4_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, c[16];
    PyObject *type, *res;
    if (!PyArg_ParseTuple(args, "dddddddddddddddddddd", &k0, &k1, &k2, &k3, c, c + 1, c + 2, c + 3,
                          c + 4, c + 5, c + 6, c + 7, c + 8, c + 9, c + 10, c + 11,
                          c + 12, c + 13, c + 14, c + 15)) {
        return NULL;
    }
    type = create_type(k0, k1, k2, k3);
    if (type == NULL) {
        return NULL;
    }
    res = cd4_from_value((PyTypeObject *)type, cayleydickson4_make(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13], c[14], c[15]));
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"cayleydickson4", cd4_factory, METH_VARARGS,
     "cayleydickson4(k0, k1, k2, k3): scalar type / dtype with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2, e3*e3 = k3."},
    {"_rebuild", cd4_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cayleydickson4",
    "Cayley-Dickson number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cayleydickson4(void)
{
    PyObject *m;

    import_array();
    import_umath();

    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }
    MetaType.tp_base = &PyType_Type;  /* not a constant initializer on MSVC */
    if (PyType_Ready(&MetaType) < 0) {
        Py_DECREF(m);
        return NULL;
    }
    PyArray_InitArrFuncs(&npycd4_arrfuncs);
    npycd4_arrfuncs.getitem = npycd4_getitem;
    npycd4_arrfuncs.setitem = npycd4_setitem;
    npycd4_arrfuncs.copyswapn = npycd4_copyswapn;
    npycd4_arrfuncs.copyswap = npycd4_copyswap;
    npycd4_arrfuncs.nonzero = npycd4_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
