/*
 * Cayley-Dickson numbers with 32 coefficients over the basis e_S (S a subset of
 * {0,1,2,3,4}, e_i*e_i = k_i), exposed to Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `cayleydickson5(k0, k1, k2, k3, k4)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `cayleydickson5.h`.
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

#include "cayleydickson5.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec5 k;
} Info;

typedef struct {
    PyObject_HEAD
    cayleydickson5 v;
} PyCD5;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.cayleydickson5.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cd5_from_value(PyTypeObject *type, cayleydickson5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCD5 *p = (PyCD5 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a cayleydickson5 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cd5_from_object(Info *info, PyObject *o, cayleydickson5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCD5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cayleydickson5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycd5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return cd5_from_value(type, cayleydickson5_make(d));
}

static PyObject *
pycd5_repr(PyObject *self)
{
    cayleydickson5 v = ((PyCD5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cayleydickson5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycd5_hash(PyObject *self)
{
    cayleydickson5 v = ((PyCD5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CAYLEYDICKSON5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cayleydickson5_is_real(v)) {
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

/* Resolve the type info of whichever operand is a cayleydickson5. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyTypeObject *type;
    /* A subclass may override __new__; recognize its algebra base type. */
    for (type = Py_TYPE(a); type != NULL; type = type->tp_base) {
        if (PyType_GetSlot(type, Py_tp_new) == (void *)pycd5_new) {
            return info_of_type(Py_TYPE(a));
        }
    }
    return info_of_type(Py_TYPE(b));
}

static PyObject *
pycd5_richcompare(PyObject *a, PyObject *b, int op)
{
    cayleydickson5 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cd5_from_object(info, a, &x)) <= 0 ||
            (r = cd5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(cayleydickson5_equal(x, y) == (op == Py_EQ));
}

#define CD5_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycd5_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        cayleydickson5 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cd5_from_object(info, a, &x)) <= 0 ||                     \
                (r = cd5_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cd5_from_value((PyTypeObject *)info->type, expr);           \
    }

CD5_BINARY(add, cayleydickson5_add(x, y))
CD5_BINARY(sub, cayleydickson5_sub(x, y))
CD5_BINARY(mul, cayleydickson5_mul(info->k, x, y))
CD5_BINARY(div, cayleydickson5_div(info->k, x, y))

static PyObject *
pycd5_pos(PyObject *self)
{
    return cd5_from_value(Py_TYPE(self), cayleydickson5_pos(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_neg(PyObject *self)
{
    return cd5_from_value(Py_TYPE(self), cayleydickson5_neg(((PyCD5 *)self)->v));
}

static int
pycd5_bool(PyObject *self)
{
    return vec32_any_nonzero(((PyCD5 *)self)->v);
}

static PyObject *
pycd5_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd5_from_value(Py_TYPE(self),
                          cayleydickson5_conjugate(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd5_from_value((PyTypeObject *)cls, cayleydickson5_zero());
}

static PyObject *
pycd5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cd5_from_value((PyTypeObject *)cls, cayleydickson5_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycd5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson5_lmat(info->k, ((PyCD5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson5_rmat(info->k, ((PyCD5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson5_minpoly(info->k, ((PyCD5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycd5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
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
    cayleydickson5_eivals(info->k, ((PyCD5 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pycd5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    cayleydickson5 v = ((PyCD5 *)self)->v;
    if (info == NULL) return NULL;
    return Py_BuildValue("(O(ddddddddddddddddddddddddddddddddddddd))", rebuild_func,
        info->k.data[0], info->k.data[1], info->k.data[2], info->k.data[3], info->k.data[4],
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pycd5_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CAYLEYDICKSON5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson5_basis((cayleydickson5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd5_gens(PyObject *cls, PyObject *Py_UNUSED(args))
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
    cayleydickson5_gens((cayleydickson5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd5_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CAYLEYDICKSON5_N, ALGPY_CAYLEYDICKSON5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    cayleydickson5_coeffs(
        info->k, (cayleydickson5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycd5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = cayleydickson5_norm(info->k, ((PyCD5 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycd5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_inv(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_square(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_exp(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_expm1(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_log(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_log1p(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_sqrt(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_sin(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_cos(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_tan(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_asin(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_acos(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_atan(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_sinh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_cosh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_tanh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_asinh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_acosh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_atanh(info->k, ((PyCD5 *)self)->v));
}

static PyObject *
pycd5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd5_from_value(Py_TYPE(self), cayleydickson5_re(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cd5_from_value(Py_TYPE(self), cayleydickson5_im(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_real(PyObject *self, void *closure)
{
    return pycd5_re(self, NULL);
}

static PyObject *
pycd5_imag(PyObject *self, void *closure)
{
    return pycd5_im(self, NULL);
}

static PyObject *pycd5_grades(PyObject *self, void *closure);

static PyGetSetDef pycd5_getset[] = {
    {"grade", pycd5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycd5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycd5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycd5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cayleydickson5_is_real(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cayleydickson5_to_re(((PyCD5 *)self)->v));
}

static PyObject *
pycd5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cd5_from_value((PyTypeObject *)cls, cayleydickson5_from_re(a));
}

static PyObject *
pycd5_lpow(PyObject *self, PyObject *other)
{
    cayleydickson5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self),
        cayleydickson5_lpow(info->k, ((PyCD5 *)self)->v, b));
}

static PyObject *
pycd5_rpow(PyObject *self, PyObject *other)
{
    cayleydickson5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cd5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self),
        cayleydickson5_rpow(info->k, ((PyCD5 *)self)->v, b));
}

static PyObject *
pycd5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cayleydickson5 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cd5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self),
        cayleydickson5_mix(info->k, ((PyCD5 *)self)->v, b, t));
}

static PyObject *
pycd5_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pycd5_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pycd5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cd5_from_value(Py_TYPE(self), cayleydickson5_ipow(info->k, ((PyCD5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycd5_grade(PyObject *self, Py_ssize_t k)
{
    cayleydickson5 v = ((PyCD5 *)self)->v;
    switch (k) {
        case 0:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade0(v));
        case 1:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade1(v));
        case 2:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade2(v));
        case 3:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade3(v));
        case 4:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade4(v));
        case 5:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_grade5(v));
        default:
            return cd5_from_value(Py_TYPE(self), cayleydickson5_zero());
    }
}

static PyObject *
pycd5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycd5_grade, 6);
}

static PyMethodDef pycd5_methods[] = {
    {"to_float", (PyCFunction)pycd5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycd5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycd5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycd5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycd5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycd5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pycd5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycd5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycd5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycd5_inv, METH_NOARGS, "inv."},
    {"square", pycd5_square, METH_NOARGS, "square."},
    {"exp", pycd5_exp, METH_NOARGS, "exp."},
    {"expm1", pycd5_expm1, METH_NOARGS, "expm1."},
    {"log", pycd5_log, METH_NOARGS, "log."},
    {"log1p", pycd5_log1p, METH_NOARGS, "log1p."},
    {"sqrt", pycd5_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pycd5_sin, METH_NOARGS, "sin."},
    {"cos", pycd5_cos, METH_NOARGS, "cos."},
    {"tan", pycd5_tan, METH_NOARGS, "tan."},
    {"asin", pycd5_asin, METH_NOARGS, "asin."},
    {"acos", pycd5_acos, METH_NOARGS, "acos."},
    {"atan", pycd5_atan, METH_NOARGS, "atan."},
    {"sinh", pycd5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycd5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycd5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycd5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycd5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycd5_atanh, METH_NOARGS, "atanh."},
    {"re", pycd5_re, METH_NOARGS, "re."},
    {"im", pycd5_im, METH_NOARGS, "im."},
    {"arcsin", pycd5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycd5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycd5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycd5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycd5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycd5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycd5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycd5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycd5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycd5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pycd5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycd5_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycd5_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pycd5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycd5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycd5_minpoly, METH_NOARGS, "Coefficients of the generic minimal polynomial."},
    {"eivals", pycd5_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycd5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycd5_length(PyObject *self)
{
    return 32;
}

static PyObject *
pycd5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCD5 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycd5_subscript(PyObject *self, PyObject *key)
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
    return pycd5_item(self, i);
}

static int
pycd5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycd5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycd5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cayleydickson5 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cayleydickson5 does not support modular power");
        return NULL;
    }
    if ((r = cd5_from_object(info, a, &x)) <= 0 ||
            (r = cd5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cd5_from_value((PyTypeObject *)info->type, cayleydickson5_rpow(info->k, x, y));
}

static PyType_Slot pycd5_slots[] = {
    {Py_tp_dealloc, pycd5_dealloc},
    {Py_tp_traverse, pycd5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycd5_repr},
    {Py_tp_str, pycd5_repr},
    {Py_tp_hash, pycd5_hash},
    {Py_tp_doc, (void *)"Cayley-Dickson number with float64 components"},
    {Py_tp_richcompare, pycd5_richcompare},
    {Py_tp_methods, pycd5_methods},
    {Py_tp_getset, pycd5_getset},
    {Py_sq_length, pycd5_length},
    {Py_sq_item, pycd5_item},
    {Py_mp_subscript, pycd5_subscript},
    {Py_tp_new, pycd5_new},
    {Py_nb_add, pycd5_add},
    {Py_nb_subtract, pycd5_sub},
    {Py_nb_multiply, pycd5_mul},
    {Py_nb_true_divide, pycd5_div},
    {Py_nb_power, (void *)pycd5_pow},
    {Py_nb_positive, pycd5_pos},
    {Py_nb_negative, pycd5_neg},
    {Py_nb_bool, pycd5_bool},
    {0, NULL}
};

static PyType_Spec pycd5_spec = {
    .name = "algpy.algebras.cayleydickson5",
    .basicsize = sizeof(PyCD5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycd5_slots,
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
npycd5_getitem(void *data, void *arr)
{
    cayleydickson5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return cd5_from_value(descr->typeobj, v);
}

static int
npycd5_setitem(PyObject *item, void *data, void *arr)
{
    cayleydickson5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cd5_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cayleydickson5");
        double d[32];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cayleydickson5 needs exactly 32 components");
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
        v = cayleydickson5_make(d);
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
npycd5_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cayleydickson5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycd5_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(cayleydickson5));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycd5_nonzero(void *data, void *arr)
{
    cayleydickson5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycd5_arrfuncs;

typedef struct { char c; cayleydickson5 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CD5_CAST_FROM(name, type, convert)                                          \
    static void                                                            \
    npycast_##name##_to_cd5(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        cayleydickson5 *dst = (cayleydickson5 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = cayleydickson5_from_re(convert(src[i])); \
        }                                                                  \
    }

#define CD5_REAL_TYPES(X) \
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

CD5_REAL_TYPES(CD5_CAST_FROM)

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
                npycast_##name##_to_cd5) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CD5_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec5 *k;
    cayleydickson5 (*unary)(cayleydickson5);
    cayleydickson5 (*unary_k)(vec5, cayleydickson5);
    cayleydickson5 (*binary)(vec5, cayleydickson5, cayleydickson5);
    int (*pred)(cayleydickson5, cayleydickson5);
} Op;

static cayleydickson5 op_add(vec5 k, cayleydickson5 a, cayleydickson5 b)
{ return cayleydickson5_add(a, b); }
static cayleydickson5 op_sub(vec5 k, cayleydickson5 a, cayleydickson5 b)
{ return cayleydickson5_sub(a, b); }
static int op_equal(cayleydickson5 a, cayleydickson5 b)
{ return cayleydickson5_equal(a, b); }
static int op_not_equal(cayleydickson5 a, cayleydickson5 b)
{ return !cayleydickson5_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        cayleydickson5 a, r;
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
        cayleydickson5 a, r;
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
        cayleydickson5 a, b, r;
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
        cayleydickson5 a, b;
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
        BIN("multiply", cayleydickson5_mul) ||
        BIN("divide", cayleydickson5_div) ||
        BIN("power", cayleydickson5_rpow) ||
        UN("positive", cayleydickson5_pos) ||
        UN("negative", cayleydickson5_neg) ||
        UN("conjugate", cayleydickson5_conjugate) ||
        UNK("reciprocal", cayleydickson5_inv) ||
        UNK("square", cayleydickson5_square) ||
        UNK("exp", cayleydickson5_exp) ||
        UNK("expm1", cayleydickson5_expm1) ||
        UNK("log", cayleydickson5_log) ||
        UNK("log1p", cayleydickson5_log1p) ||
        UNK("sqrt", cayleydickson5_sqrt) ||
        UNK("sin", cayleydickson5_sin) ||
        UNK("cos", cayleydickson5_cos) ||
        UNK("tan", cayleydickson5_tan) ||
        UNK("arcsin", cayleydickson5_asin) ||
        UNK("arccos", cayleydickson5_acos) ||
        UNK("arctan", cayleydickson5_atan) ||
        UNK("sinh", cayleydickson5_sinh) ||
        UNK("cosh", cayleydickson5_cosh) ||
        UNK("tanh", cayleydickson5_tanh) ||
        UNK("arcsinh", cayleydickson5_asinh) ||
        UNK("arccosh", cayleydickson5_acosh) ||
        UNK("arctanh", cayleydickson5_atanh) ||
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
    .tp_name = "algpy.cayleydickson5.type",
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

    type = PyType_FromSpecWithBases(&pycd5_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.cayleydickson5.info", NULL);
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
    proto->elsize = sizeof(cayleydickson5);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycd5_arrfuncs;
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
cd5_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4;
    if (!PyArg_ParseTuple(args, "ddddd", &k0, &k1, &k2, &k3, &k4)) return NULL;
    return create_type(k0, k1, k2, k3, k4);
}

static PyObject *
cd5_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4, c[32];
    PyObject *type, *result;
    if (!PyArg_ParseTuple(args, "ddddddddddddddddddddddddddddddddddddd", &k0, &k1, &k2, &k3, &k4,
        c + 0, c + 1, c + 2, c + 3, c + 4, c + 5, c + 6, c + 7, c + 8, c + 9, c + 10, c + 11, c + 12, c + 13, c + 14, c + 15, c + 16, c + 17, c + 18, c + 19, c + 20, c + 21, c + 22, c + 23, c + 24, c + 25, c + 26, c + 27, c + 28, c + 29, c + 30, c + 31)) return NULL;
    type = create_type(k0, k1, k2, k3, k4);
    if (type == NULL) return NULL;
    result = cd5_from_value((PyTypeObject *)type, cayleydickson5_make(c));
    Py_DECREF(type);
    return result;
}

static PyMethodDef module_methods[] = {
    {"cayleydickson5", cd5_factory, METH_VARARGS,
     "cayleydickson5(k0, k1, k2, k3, k4): scalar type / dtype with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2, e3*e3 = k3, e4*e4 = k4."},
    {"_rebuild", cd5_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cayleydickson5",
    "Cayley-Dickson number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cayleydickson5(void)
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
    PyArray_InitArrFuncs(&npycd5_arrfuncs);
    npycd5_arrfuncs.getitem = npycd5_getitem;
    npycd5_arrfuncs.setitem = npycd5_setitem;
    npycd5_arrfuncs.copyswapn = npycd5_copyswapn;
    npycd5_arrfuncs.copyswap = npycd5_copyswap;
    npycd5_arrfuncs.nonzero = npycd5_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
