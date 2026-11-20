/*
 * Multicomplex numbers with 32 coefficients over the basis e_S (S a subset of
 * {0,1,2,3,4}, e_i*e_i = k_i), exposed to Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `multicomplex5(k0, k1, k2, k3, k4)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `multicomplex5.h`.
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
#include "../../../_grades.h"

#include "multicomplex5.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec5 k;
} Info;

typedef struct {
    PyObject_HEAD
    multicomplex5 v;
} PyMC5;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.multicomplex5.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
mc5_from_value(PyTypeObject *type, multicomplex5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyMC5 *p = (PyMC5 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a multicomplex5 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
mc5_from_object(Info *info, PyObject *o, multicomplex5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyMC5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = multicomplex5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pymc5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return mc5_from_value(type, multicomplex5_make(d));
}

static PyObject *
pymc5_repr(PyObject *self)
{
    multicomplex5 v = ((PyMC5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("multicomplex5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pymc5_hash(PyObject *self)
{
    multicomplex5 v = ((PyMC5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_MULTICOMPLEX5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (multicomplex5_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddddddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

/* Resolve the type info of whichever operand is a multicomplex5. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyTypeObject *type;
    /* A subclass may override __new__; recognize its algebra base type. */
    for (type = Py_TYPE(a); type != NULL; type = type->tp_base) {
        if (PyType_GetSlot(type, Py_tp_new) == (void *)pymc5_new) {
            return info_of_type(Py_TYPE(a));
        }
    }
    return info_of_type(Py_TYPE(b));
}

static PyObject *
pymc5_richcompare(PyObject *a, PyObject *b, int op)
{
    multicomplex5 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = mc5_from_object(info, a, &x)) <= 0 ||
            (r = mc5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(multicomplex5_equal(x, y) == (op == Py_EQ));
}

#define MC5_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pymc5_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        multicomplex5 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = mc5_from_object(info, a, &x)) <= 0 ||                     \
                (r = mc5_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return mc5_from_value((PyTypeObject *)info->type, expr);           \
    }

MC5_BINARY(add, multicomplex5_add(x, y))
MC5_BINARY(sub, multicomplex5_sub(x, y))
MC5_BINARY(mul, multicomplex5_mul(info->k, x, y))
MC5_BINARY(div, multicomplex5_div(info->k, x, y))

static PyObject *
pymc5_pos(PyObject *self)
{
    return mc5_from_value(Py_TYPE(self), multicomplex5_pos(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_neg(PyObject *self)
{
    return mc5_from_value(Py_TYPE(self), multicomplex5_neg(((PyMC5 *)self)->v));
}

static int
pymc5_bool(PyObject *self)
{
    return vec32_any_nonzero(((PyMC5 *)self)->v);
}

static PyObject *
pymc5_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc5_from_value(Py_TYPE(self),
                          multicomplex5_conjugate(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc5_from_value((PyTypeObject *)cls, multicomplex5_zero());
}

static PyObject *
pymc5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc5_from_value((PyTypeObject *)cls, multicomplex5_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pymc5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {32, 32};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_lmat(info->k, ((PyMC5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {32, 32};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_rmat(info->k, ((PyMC5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {33};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_minpoly(info->k, ((PyMC5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {32};
    complex ev[32];
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
    multicomplex5_eivals(info->k, ((PyMC5 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 32; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pymc5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    multicomplex5 v = ((PyMC5 *)self)->v;
    if (info == NULL) return NULL;
    return Py_BuildValue("(O(ddddddddddddddddddddddddddddddddddddd))", rebuild_func,
        info->k.data[0], info->k.data[1], info->k.data[2], info->k.data[3], info->k.data[4],
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pymc5_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_MULTICOMPLEX5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_basis((multicomplex5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc5_gens(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {5};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_gens((multicomplex5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc5_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_MULTICOMPLEX5_N, ALGPY_MULTICOMPLEX5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex5_coeffs(
        info->k, (multicomplex5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = multicomplex5_norm(info->k, ((PyMC5 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pymc5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_inv(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_square(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_exp(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_expm1(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_log(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_log1p(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_sqrt(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_sin(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_cos(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_tan(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_asin(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_acos(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_atan(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_sinh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_cosh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_tanh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_asinh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_acosh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_atanh(info->k, ((PyMC5 *)self)->v));
}

static PyObject *
pymc5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc5_from_value(Py_TYPE(self), multicomplex5_re(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc5_from_value(Py_TYPE(self), multicomplex5_im(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_real(PyObject *self, void *closure)
{
    return pymc5_re(self, NULL);
}

static PyObject *
pymc5_imag(PyObject *self, void *closure)
{
    return pymc5_im(self, NULL);
}

static PyObject *pymc5_grades(PyObject *self, void *closure);

static PyGetSetDef pymc5_getset[] = {
    {"grade", pymc5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pymc5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pymc5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pymc5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(multicomplex5_is_real(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(multicomplex5_to_re(((PyMC5 *)self)->v));
}

static PyObject *
pymc5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return mc5_from_value((PyTypeObject *)cls, multicomplex5_from_re(a));
}

static PyObject *
pymc5_lpow(PyObject *self, PyObject *other)
{
    multicomplex5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = mc5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self),
        multicomplex5_lpow(info->k, ((PyMC5 *)self)->v, b));
}

static PyObject *
pymc5_rpow(PyObject *self, PyObject *other)
{
    multicomplex5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = mc5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self),
        multicomplex5_rpow(info->k, ((PyMC5 *)self)->v, b));
}

static PyObject *
pymc5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    multicomplex5 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = mc5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self),
        multicomplex5_mix(info->k, ((PyMC5 *)self)->v, b, t));
}

static PyObject *
pymc5_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pymc5_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pymc5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc5_from_value(Py_TYPE(self), multicomplex5_ipow(info->k, ((PyMC5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pymc5_grade(PyObject *self, Py_ssize_t k)
{
    multicomplex5 v = ((PyMC5 *)self)->v;
    switch (k) {
        case 0:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade0(v));
        case 1:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade1(v));
        case 2:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade2(v));
        case 3:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade3(v));
        case 4:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade4(v));
        case 5:
            return mc5_from_value(Py_TYPE(self), multicomplex5_grade5(v));
        default:
            return mc5_from_value(Py_TYPE(self), multicomplex5_zero());
    }
}

static PyObject *
pymc5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pymc5_grade, 6);
}

static PyMethodDef pymc5_methods[] = {
    {"to_float", (PyCFunction)pymc5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pymc5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pymc5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pymc5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pymc5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pymc5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pymc5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pymc5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pymc5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pymc5_inv, METH_NOARGS, "inv."},
    {"square", pymc5_square, METH_NOARGS, "square."},
    {"exp", pymc5_exp, METH_NOARGS, "exp."},
    {"expm1", pymc5_expm1, METH_NOARGS, "expm1."},
    {"log", pymc5_log, METH_NOARGS, "log."},
    {"log1p", pymc5_log1p, METH_NOARGS, "log1p."},
    {"sqrt", pymc5_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pymc5_sin, METH_NOARGS, "sin."},
    {"cos", pymc5_cos, METH_NOARGS, "cos."},
    {"tan", pymc5_tan, METH_NOARGS, "tan."},
    {"asin", pymc5_asin, METH_NOARGS, "asin."},
    {"acos", pymc5_acos, METH_NOARGS, "acos."},
    {"atan", pymc5_atan, METH_NOARGS, "atan."},
    {"sinh", pymc5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pymc5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pymc5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pymc5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pymc5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pymc5_atanh, METH_NOARGS, "atanh."},
    {"re", pymc5_re, METH_NOARGS, "re."},
    {"im", pymc5_im, METH_NOARGS, "im."},
    {"arcsin", pymc5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pymc5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pymc5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pymc5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pymc5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pymc5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pymc5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pymc5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pymc5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pymc5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pymc5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pymc5_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pymc5_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pymc5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pymc5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pymc5_minpoly, METH_NOARGS, "Coefficients of the generic minimal polynomial."},
    {"eivals", pymc5_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pymc5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pymc5_length(PyObject *self)
{
    return 32;
}

static PyObject *
pymc5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyMC5 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pymc5_subscript(PyObject *self, PyObject *key)
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
        i += 32;
    }
    return pymc5_item(self, i);
}

static int
pymc5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pymc5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pymc5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    multicomplex5 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "multicomplex5 does not support modular power");
        return NULL;
    }
    if ((r = mc5_from_object(info, a, &x)) <= 0 ||
            (r = mc5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return mc5_from_value((PyTypeObject *)info->type, multicomplex5_rpow(info->k, x, y));
}

static PyType_Slot pymc5_slots[] = {
    {Py_tp_dealloc, pymc5_dealloc},
    {Py_tp_traverse, pymc5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pymc5_repr},
    {Py_tp_str, pymc5_repr},
    {Py_tp_hash, pymc5_hash},
    {Py_tp_doc, (void *)"Multicomplex number with float64 components"},
    {Py_tp_richcompare, pymc5_richcompare},
    {Py_tp_methods, pymc5_methods},
    {Py_tp_getset, pymc5_getset},
    {Py_sq_length, pymc5_length},
    {Py_sq_item, pymc5_item},
    {Py_mp_subscript, pymc5_subscript},
    {Py_tp_new, pymc5_new},
    {Py_nb_add, pymc5_add},
    {Py_nb_subtract, pymc5_sub},
    {Py_nb_multiply, pymc5_mul},
    {Py_nb_true_divide, pymc5_div},
    {Py_nb_power, (void *)pymc5_pow},
    {Py_nb_positive, pymc5_pos},
    {Py_nb_negative, pymc5_neg},
    {Py_nb_bool, pymc5_bool},
    {0, NULL}
};

static PyType_Spec pymc5_spec = {
    .name = "algpy.algebras.multicomplex5",
    .basicsize = sizeof(PyMC5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pymc5_slots,
};

/* NumPy dtype support */

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npymc5_getitem(void *data, void *arr)
{
    multicomplex5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return mc5_from_value(descr->typeobj, v);
}

static int
npymc5_setitem(PyObject *item, void *data, void *arr)
{
    multicomplex5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = mc5_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to multicomplex5");
        double d[32];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "multicomplex5 needs exactly 32 components");
            return -1;
        }
        for (i = 0; i < 32; i++) {
            d[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (d[i] == -1.0 && PyErr_Occurred()) {
                Py_DECREF(seq);
                return -1;
            }
        }
        Py_DECREF(seq);
        v = multicomplex5_make(d);
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
npymc5_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(multicomplex5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npymc5_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(multicomplex5));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npymc5_nonzero(void *data, void *arr)
{
    multicomplex5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npymc5_arrfuncs;

typedef struct { char c; multicomplex5 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define MC5_CAST_FROM(name, type, convert)                                          \
    static void                                                            \
    npycast_##name##_to_mc5(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        multicomplex5 *dst = (multicomplex5 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = multicomplex5_from_re(convert(src[i])); \
        }                                                                  \
    }

#define MC5_REAL_TYPES(X) \
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

MC5_REAL_TYPES(MC5_CAST_FROM)

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
                npycast_##name##_to_mc5) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    MC5_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec5 *k;
    multicomplex5 (*unary)(multicomplex5);
    multicomplex5 (*unary_k)(vec5, multicomplex5);
    multicomplex5 (*binary)(vec5, multicomplex5, multicomplex5);
    int (*pred)(multicomplex5, multicomplex5);
} Op;

static multicomplex5 op_add(vec5 k, multicomplex5 a, multicomplex5 b)
{ return multicomplex5_add(a, b); }
static multicomplex5 op_sub(vec5 k, multicomplex5 a, multicomplex5 b)
{ return multicomplex5_sub(a, b); }
static int op_equal(multicomplex5 a, multicomplex5 b)
{ return multicomplex5_equal(a, b); }
static int op_not_equal(multicomplex5 a, multicomplex5 b)
{ return !multicomplex5_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        multicomplex5 a, r;
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
        multicomplex5 a, r;
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
        multicomplex5 a, b, r;
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
        multicomplex5 a, b;
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
        BIN("multiply", multicomplex5_mul) ||
        BIN("divide", multicomplex5_div) ||
        BIN("power", multicomplex5_rpow) ||
        UN("positive", multicomplex5_pos) ||
        UN("negative", multicomplex5_neg) ||
        UN("conjugate", multicomplex5_conjugate) ||
        UNK("reciprocal", multicomplex5_inv) ||
        UNK("square", multicomplex5_square) ||
        UNK("exp", multicomplex5_exp) ||
        UNK("expm1", multicomplex5_expm1) ||
        UNK("log", multicomplex5_log) ||
        UNK("log1p", multicomplex5_log1p) ||
        UNK("sqrt", multicomplex5_sqrt) ||
        UNK("sin", multicomplex5_sin) ||
        UNK("cos", multicomplex5_cos) ||
        UNK("tan", multicomplex5_tan) ||
        UNK("arcsin", multicomplex5_asin) ||
        UNK("arccos", multicomplex5_acos) ||
        UNK("arctan", multicomplex5_atan) ||
        UNK("sinh", multicomplex5_sinh) ||
        UNK("cosh", multicomplex5_cosh) ||
        UNK("tanh", multicomplex5_tanh) ||
        UNK("arcsinh", multicomplex5_asinh) ||
        UNK("arccosh", multicomplex5_acosh) ||
        UNK("arctanh", multicomplex5_atanh) ||
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
    list = PyList_New(5);
    for (i = 0; list != NULL && i < 5; i++) {
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
    vec5 k;
    int i;
    if (info == NULL) {
        return -1;
    }
    if (value == NULL) {
        PyErr_SetString(PyExc_AttributeError, "cannot delete K");
        return -1;
    }
    seq = PySequence_Fast(value, "K must be a sequence of 5 numbers");
    if (seq == NULL) {
        return -1;
    }
    if (PySequence_Fast_GET_SIZE(seq) != 5) {
        PyErr_SetString(PyExc_ValueError, "K must have exactly 5 element(s)");
        Py_DECREF(seq);
        return -1;
    }
    for (i = 0; i < 5; i++) {
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
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 5 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.multicomplex5.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = meta_getset,
};

/* Type factory */

static PyObject *
create_type(double k0, double k1, double k2, double k3, double k4)
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
    info->k.data[4] = k4;

    type = PyType_FromSpecWithBases(&pymc5_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.multicomplex5.info", NULL);
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
    proto->elsize = sizeof(multicomplex5);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npymc5_arrfuncs;
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

    if (algpy_set_dimension(type, 32) < 0) {
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
mc5_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4;
    if (!PyArg_ParseTuple(args, "ddddd", &k0, &k1, &k2, &k3, &k4)) return NULL;
    return create_type(k0, k1, k2, k3, k4);
}

static PyObject *
mc5_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4, c[32];
    PyObject *type, *result;
    if (!PyArg_ParseTuple(args, "ddddddddddddddddddddddddddddddddddddd", &k0, &k1, &k2, &k3, &k4,
        c + 0, c + 1, c + 2, c + 3, c + 4, c + 5, c + 6, c + 7, c + 8, c + 9, c + 10, c + 11, c + 12, c + 13, c + 14, c + 15, c + 16, c + 17, c + 18, c + 19, c + 20, c + 21, c + 22, c + 23, c + 24, c + 25, c + 26, c + 27, c + 28, c + 29, c + 30, c + 31)) return NULL;
    type = create_type(k0, k1, k2, k3, k4);
    if (type == NULL) return NULL;
    result = mc5_from_value((PyTypeObject *)type, multicomplex5_make(c));
    Py_DECREF(type);
    return result;
}

static PyMethodDef module_methods[] = {
    {"multicomplex5", mc5_factory, METH_VARARGS,
     "multicomplex5(k0, k1, k2, k3, k4): scalar type / dtype with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2, e3*e3 = k3, e4*e4 = k4."},
    {"_rebuild", mc5_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_multicomplex5",
    "Multicomplex number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__multicomplex5(void)
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
    PyArray_InitArrFuncs(&npymc5_arrfuncs);
    npymc5_arrfuncs.getitem = npymc5_getitem;
    npymc5_arrfuncs.setitem = npymc5_setitem;
    npymc5_arrfuncs.copyswapn = npymc5_copyswapn;
    npymc5_arrfuncs.copyswap = npymc5_copyswap;
    npymc5_arrfuncs.nonzero = npymc5_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
