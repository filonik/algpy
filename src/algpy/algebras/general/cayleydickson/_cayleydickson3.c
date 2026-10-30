/*
 * Cayley-Dickson numbers x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 +
 * x12*e1*e2 + x012*e0*e1*e2 (e0*e0 = k0, e1*e1 = k1, e2*e2 = k2) exposed to
 * Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `cayleydickson3(k0, k1, k2)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `cayleydickson3.h`.
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

#include "cayleydickson3.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec3 k;
} Info;

typedef struct {
    PyObject_HEAD
    cayleydickson3 v;
} PyCD3;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.cayleydickson3.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cd3_from_value(PyTypeObject *type, cayleydickson3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCD3 *p = (PyCD3 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a cayleydickson3 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cd3_from_object(Info *info, PyObject *o, cayleydickson3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCD3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cayleydickson3_make(d, 0, 0, 0, 0, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycd3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12",
                             "x012", NULL};
    double c[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
                                     c, c + 1, c + 2, c + 3, c + 4, c + 5,
                                     c + 6, c + 7)) {
        return NULL;
    }
    return cd3_from_value(type, cayleydickson3_make(c[0], c[1], c[2], c[3],
                                                   c[4], c[5], c[6], c[7]));
}

static PyObject *
pycd3_repr(PyObject *self)
{
    cayleydickson3 v = ((PyCD3 *)self)->v;
    char *s[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cayleydickson3(%s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6],
                                  s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycd3_hash(PyObject *self)
{
    cayleydickson3 v = ((PyCD3 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CAYLEYDICKSON3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cayleydickson3_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddd)", v.data[0], v.data[1], v.data[2],
                            v.data[3], v.data[4], v.data[5], v.data[6],
                            v.data[7]);
    }
    if (obj == NULL) {
        return -1;
    }
    h = PyObject_Hash(obj);
    Py_DECREF(obj);
    return h;
}

/* Resolve the type info of whichever operand is a cayleydickson3. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyTypeObject *type;
    /* A subclass may override __new__; recognize its algebra base type. */
    for (type = Py_TYPE(a); type != NULL; type = type->tp_base) {
        if (PyType_GetSlot(type, Py_tp_new) == (void *)pycd3_new) {
            return info_of_type(Py_TYPE(a));
        }
    }
    return info_of_type(Py_TYPE(b));
}

static PyObject *
pycd3_richcompare(PyObject *a, PyObject *b, int op)
{
    cayleydickson3 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cd3_from_object(info, a, &x)) <= 0 ||
            (r = cd3_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(cayleydickson3_equal(x, y) == (op == Py_EQ));
}

#define CD3_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycd3_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        cayleydickson3 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cd3_from_object(info, a, &x)) <= 0 ||                     \
                (r = cd3_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cd3_from_value((PyTypeObject *)info->type, expr);           \
    }

CD3_BINARY(add, cayleydickson3_add(x, y))
CD3_BINARY(sub, cayleydickson3_sub(x, y))
CD3_BINARY(mul, cayleydickson3_mul(info->k, x, y))
CD3_BINARY(div, cayleydickson3_div(info->k, x, y))

static PyObject *
pycd3_pos(PyObject *self)
{
    return cd3_from_value(Py_TYPE(self), cayleydickson3_pos(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_neg(PyObject *self)
{
    return cd3_from_value(Py_TYPE(self), cayleydickson3_neg(((PyCD3 *)self)->v));
}

static int
pycd3_bool(PyObject *self)
{
    return vec8_any_nonzero(((PyCD3 *)self)->v);
}

static PyObject *
pycd3_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd3_from_value(Py_TYPE(self),
                          cayleydickson3_conjugate(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd3_from_value((PyTypeObject *)cls, cayleydickson3_zero());
}

static PyObject *
pycd3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd3_from_value((PyTypeObject *)cls, cayleydickson3_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycd3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {8, 8};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson3_lmat(info->k, ((PyCD3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {8, 8};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson3_rmat(info->k, ((PyCD3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson3_minpoly(info->k, ((PyCD3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson3_eivals(info->k, ((PyCD3 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pycd3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    cayleydickson3 v = ((PyCD3 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(ddddddddddd))", rebuild_func, info->k.data[0],
                         info->k.data[1], info->k.data[2], v.data[0],
                         v.data[1], v.data[2], v.data[3], v.data[4],
                         v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pycd3_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CAYLEYDICKSON3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson3_basis((cayleydickson3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd3_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CAYLEYDICKSON3_N, ALGPY_CAYLEYDICKSON3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson3_coeffs(
        info->k, (cayleydickson3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = cayleydickson3_norm(info->k, ((PyCD3 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycd3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_inv(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_square(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_exp(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_expm1(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_log(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_log1p(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_sqrt(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_sin(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_cos(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_tan(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_asin(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_acos(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_atan(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_sinh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_cosh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_tanh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_asinh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_acosh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_atanh(info->k, ((PyCD3 *)self)->v));
}

static PyObject *
pycd3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd3_from_value(Py_TYPE(self), cayleydickson3_re(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd3_from_value(Py_TYPE(self), cayleydickson3_im(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_real(PyObject *self, void *closure)
{
    return pycd3_re(self, NULL);
}

static PyObject *
pycd3_imag(PyObject *self, void *closure)
{
    return pycd3_im(self, NULL);
}

static PyGetSetDef pycd3_getset[] = {
    {"real", pycd3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycd3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycd3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cayleydickson3_is_real(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cayleydickson3_to_re(((PyCD3 *)self)->v));
}

static PyObject *
pycd3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cd3_from_value((PyTypeObject *)cls, cayleydickson3_from_re(a));
}

static PyObject *
pycd3_lpow(PyObject *self, PyObject *other)
{
    cayleydickson3 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd3_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self),
        cayleydickson3_lpow(info->k, ((PyCD3 *)self)->v, b));
}

static PyObject *
pycd3_rpow(PyObject *self, PyObject *other)
{
    cayleydickson3 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd3_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self),
        cayleydickson3_rpow(info->k, ((PyCD3 *)self)->v, b));
}

static PyObject *
pycd3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cayleydickson3 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cd3_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self),
        cayleydickson3_mix(info->k, ((PyCD3 *)self)->v, b, t));
}

static PyObject *
pycd3_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 8);
}

static PyObject *
pycd3_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 8);
}

static PyObject *
pycd3_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd3_from_value(Py_TYPE(self), cayleydickson3_ipow(info->k, ((PyCD3 *)self)->v, (int64_t)exponent));
}

static PyMethodDef pycd3_methods[] = {
    {"to_float", (PyCFunction)pycd3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycd3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycd3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycd3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycd3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycd3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pycd3_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycd3_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycd3_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycd3_inv, METH_NOARGS, "inv."},
    {"square", pycd3_square, METH_NOARGS, "square."},
    {"exp", pycd3_exp, METH_NOARGS, "exp."},
    {"expm1", pycd3_expm1, METH_NOARGS, "expm1."},
    {"log", pycd3_log, METH_NOARGS, "log."},
    {"log1p", pycd3_log1p, METH_NOARGS, "log1p."},
    {"sqrt", pycd3_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pycd3_sin, METH_NOARGS, "sin."},
    {"cos", pycd3_cos, METH_NOARGS, "cos."},
    {"tan", pycd3_tan, METH_NOARGS, "tan."},
    {"asin", pycd3_asin, METH_NOARGS, "asin."},
    {"acos", pycd3_acos, METH_NOARGS, "acos."},
    {"atan", pycd3_atan, METH_NOARGS, "atan."},
    {"sinh", pycd3_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycd3_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycd3_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycd3_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycd3_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycd3_atanh, METH_NOARGS, "atanh."},
    {"re", pycd3_re, METH_NOARGS, "re."},
    {"im", pycd3_im, METH_NOARGS, "im."},
    {"arcsin", pycd3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycd3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycd3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycd3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycd3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycd3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycd3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycd3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycd3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycd3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycd3_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycd3_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pycd3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycd3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycd3_minpoly, METH_NOARGS, "Coefficients of the generic minimal polynomial."},
    {"eivals", pycd3_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycd3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycd3_length(PyObject *self)
{
    return 8;
}

static PyObject *
pycd3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCD3 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycd3_subscript(PyObject *self, PyObject *key)
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
        i += 8;
    }
    return pycd3_item(self, i);
}

static int
pycd3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycd3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycd3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cayleydickson3 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cayleydickson3 does not support modular power");
        return NULL;
    }
    if ((r = cd3_from_object(info, a, &x)) <= 0 ||
            (r = cd3_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cd3_from_value((PyTypeObject *)info->type, cayleydickson3_rpow(info->k, x, y));
}

static PyType_Slot pycd3_slots[] = {
    {Py_tp_dealloc, pycd3_dealloc},
    {Py_tp_traverse, pycd3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycd3_repr},
    {Py_tp_str, pycd3_repr},
    {Py_tp_hash, pycd3_hash},
    {Py_tp_doc, (void *)"Cayley-Dickson number with float64 components"},
    {Py_tp_richcompare, pycd3_richcompare},
    {Py_tp_methods, pycd3_methods},
    {Py_tp_getset, pycd3_getset},
    {Py_sq_length, pycd3_length},
    {Py_sq_item, pycd3_item},
    {Py_mp_subscript, pycd3_subscript},
    {Py_tp_new, pycd3_new},
    {Py_nb_add, pycd3_add},
    {Py_nb_subtract, pycd3_sub},
    {Py_nb_multiply, pycd3_mul},
    {Py_nb_true_divide, pycd3_div},
    {Py_nb_power, (void *)pycd3_pow},
    {Py_nb_positive, pycd3_pos},
    {Py_nb_negative, pycd3_neg},
    {Py_nb_bool, pycd3_bool},
    {0, NULL}
};

static PyType_Spec pycd3_spec = {
    .name = "algpy.algebras.cayleydickson3",
    .basicsize = sizeof(PyCD3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycd3_slots,
};

/* NumPy dtype support */

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npycd3_getitem(void *data, void *arr)
{
    cayleydickson3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return cd3_from_value(descr->typeobj, v);
}

static int
npycd3_setitem(PyObject *item, void *data, void *arr)
{
    cayleydickson3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cd3_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cayleydickson3");
        double d[8];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cayleydickson3 needs exactly 8 components");
            return -1;
        }
        for (i = 0; i < 8; i++) {
            d[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (d[i] == -1.0 && PyErr_Occurred()) {
                Py_DECREF(seq);
                return -1;
            }
        }
        Py_DECREF(seq);
        v = cayleydickson3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
npycd3_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cayleydickson3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycd3_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(cayleydickson3));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycd3_nonzero(void *data, void *arr)
{
    cayleydickson3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycd3_arrfuncs;

typedef struct { char c; cayleydickson3 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CD3_CAST_FROM(name, type, convert)                                          \
    static void                                                            \
    npycast_##name##_to_cd3(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        cayleydickson3 *dst = (cayleydickson3 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = cayleydickson3_make(convert(src[i]), 0, 0, 0, 0, 0, 0, 0);                \
        }                                                                  \
    }

#define CD3_REAL_TYPES(X) \
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

CD3_REAL_TYPES(CD3_CAST_FROM)

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
                npycast_##name##_to_cd3) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CD3_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec3 *k;
    cayleydickson3 (*unary)(cayleydickson3);
    cayleydickson3 (*unary_k)(vec3, cayleydickson3);
    cayleydickson3 (*binary)(vec3, cayleydickson3, cayleydickson3);
    int (*pred)(cayleydickson3, cayleydickson3);
} Op;

static cayleydickson3 op_add(vec3 k, cayleydickson3 a, cayleydickson3 b)
{ return cayleydickson3_add(a, b); }
static cayleydickson3 op_sub(vec3 k, cayleydickson3 a, cayleydickson3 b)
{ return cayleydickson3_sub(a, b); }
static int op_equal(cayleydickson3 a, cayleydickson3 b)
{ return cayleydickson3_equal(a, b); }
static int op_not_equal(cayleydickson3 a, cayleydickson3 b)
{ return !cayleydickson3_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        cayleydickson3 a, r;
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
        cayleydickson3 a, r;
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
        cayleydickson3 a, b, r;
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
        cayleydickson3 a, b;
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
        BIN("multiply", cayleydickson3_mul) ||
        BIN("divide", cayleydickson3_div) ||
        BIN("power", cayleydickson3_rpow) ||
        UN("positive", cayleydickson3_pos) ||
        UN("negative", cayleydickson3_neg) ||
        UN("conjugate", cayleydickson3_conjugate) ||
        UNK("reciprocal", cayleydickson3_inv) ||
        UNK("square", cayleydickson3_square) ||
        UNK("exp", cayleydickson3_exp) ||
        UNK("expm1", cayleydickson3_expm1) ||
        UNK("log", cayleydickson3_log) ||
        UNK("log1p", cayleydickson3_log1p) ||
        UNK("sqrt", cayleydickson3_sqrt) ||
        UNK("sin", cayleydickson3_sin) ||
        UNK("cos", cayleydickson3_cos) ||
        UNK("tan", cayleydickson3_tan) ||
        UNK("arcsin", cayleydickson3_asin) ||
        UNK("arccos", cayleydickson3_acos) ||
        UNK("arctan", cayleydickson3_atan) ||
        UNK("sinh", cayleydickson3_sinh) ||
        UNK("cosh", cayleydickson3_cosh) ||
        UNK("tanh", cayleydickson3_tanh) ||
        UNK("arcsinh", cayleydickson3_asinh) ||
        UNK("arccosh", cayleydickson3_acosh) ||
        UNK("arctanh", cayleydickson3_atanh) ||
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
    list = PyList_New(3);
    for (i = 0; list != NULL && i < 3; i++) {
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
    vec3 k;
    int i;
    if (info == NULL) {
        return -1;
    }
    if (value == NULL) {
        PyErr_SetString(PyExc_AttributeError, "cannot delete K");
        return -1;
    }
    seq = PySequence_Fast(value, "K must be a sequence of 3 numbers");
    if (seq == NULL) {
        return -1;
    }
    if (PySequence_Fast_GET_SIZE(seq) != 3) {
        PyErr_SetString(PyExc_ValueError, "K must have exactly 3 element(s)");
        Py_DECREF(seq);
        return -1;
    }
    for (i = 0; i < 3; i++) {
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
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 3 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.cayleydickson3.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = meta_getset,
};

/* Type factory */

static PyObject *
create_type(double k0, double k1, double k2)
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

    type = PyType_FromSpecWithBases(&pycd3_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.cayleydickson3.info", NULL);
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
    proto->elsize = sizeof(cayleydickson3);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycd3_arrfuncs;
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

    if (algpy_set_dimension(type, 8) < 0) {
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
cd3_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2;
    if (!PyArg_ParseTuple(args, "ddd", &k0, &k1, &k2)) {
        return NULL;
    }
    return create_type(k0, k1, k2);
}

static PyObject *
cd3_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, c[8];
    PyObject *type, *res;
    if (!PyArg_ParseTuple(args, "ddddddddddd", &k0, &k1, &k2, c, c + 1, c + 2, c + 3,
                          c + 4, c + 5, c + 6, c + 7)) {
        return NULL;
    }
    type = create_type(k0, k1, k2);
    if (type == NULL) {
        return NULL;
    }
    res = cd3_from_value((PyTypeObject *)type, cayleydickson3_make(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]));
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"cayleydickson3", cd3_factory, METH_VARARGS,
     "cayleydickson3(k0, k1, k2): scalar type / dtype of numbers with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2."},
    {"_rebuild", cd3_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cayleydickson3",
    "Cayley-Dickson number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cayleydickson3(void)
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
    PyArray_InitArrFuncs(&npycd3_arrfuncs);
    npycd3_arrfuncs.getitem = npycd3_getitem;
    npycd3_arrfuncs.setitem = npycd3_setitem;
    npycd3_arrfuncs.copyswapn = npycd3_copyswapn;
    npycd3_arrfuncs.copyswap = npycd3_copyswap;
    npycd3_arrfuncs.nonzero = npycd3_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
