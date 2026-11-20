/*
 * Clifford numbers with 32 coefficients over the basis e_S (S a subset of
 * {0,1,2,3,4}, e_i*e_i = k_i), exposed to Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `clifford5(k0, k1, k2, k3, k4)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `clifford5.h`.
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

#include "clifford5.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec5 k;
    PyObject *center;  /* clifford1 type of the centre, created on first use */
} Info;

typedef struct {
    PyObject_HEAD
    clifford5 v;
} PyCL5;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.clifford5.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cl5_from_value(PyTypeObject *type, clifford5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCL5 *p = (PyCL5 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a clifford5 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cl5_from_object(Info *info, PyObject *o, clifford5 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCL5 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = clifford5_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycl5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return cl5_from_value(type, clifford5_make(d));
}

static PyObject *
pycl5_repr(PyObject *self)
{
    clifford5 v = ((PyCL5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("clifford5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycl5_hash(PyObject *self)
{
    clifford5 v = ((PyCL5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CLIFFORD5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (clifford5_is_real(v)) {
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

/* Resolve the type info of whichever operand is a clifford5. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyTypeObject *type;
    /* A subclass may override __new__; recognize its algebra base type. */
    for (type = Py_TYPE(a); type != NULL; type = type->tp_base) {
        if (PyType_GetSlot(type, Py_tp_new) == (void *)pycl5_new) {
            return info_of_type(Py_TYPE(a));
        }
    }
    return info_of_type(Py_TYPE(b));
}

static PyObject *
pycl5_richcompare(PyObject *a, PyObject *b, int op)
{
    clifford5 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cl5_from_object(info, a, &x)) <= 0 ||
            (r = cl5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(clifford5_equal(x, y) == (op == Py_EQ));
}

#define CL5_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycl5_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        clifford5 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cl5_from_object(info, a, &x)) <= 0 ||                     \
                (r = cl5_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cl5_from_value((PyTypeObject *)info->type, expr);           \
    }

CL5_BINARY(add, clifford5_add(x, y))
CL5_BINARY(sub, clifford5_sub(x, y))
CL5_BINARY(mul, clifford5_mul(info->k, x, y))
CL5_BINARY(div, clifford5_div(info->k, x, y))

static PyObject *
pycl5_pos(PyObject *self)
{
    return cl5_from_value(Py_TYPE(self), clifford5_pos(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_neg(PyObject *self)
{
    return cl5_from_value(Py_TYPE(self), clifford5_neg(((PyCL5 *)self)->v));
}

static int
pycl5_bool(PyObject *self)
{
    return vec32_any_nonzero(((PyCL5 *)self)->v);
}

static PyObject *
pycl5_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl5_from_value(Py_TYPE(self),
                          clifford5_conjugate(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl5_from_value((PyTypeObject *)cls, clifford5_zero());
}

static PyObject *
pycl5_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl5_from_value((PyTypeObject *)cls, clifford5_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycl5_lmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford5_lmat(info->k, ((PyCL5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl5_rmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford5_rmat(info->k, ((PyCL5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

/* The clifford1 type of the centre span{1, w}, w the pseudoscalar, with its
 * K set to [w*w] for the current parameters. */
static PyObject *
center_type(Info *info)
{
    double s = clifford5_pseudoscalar_square(info->k) + 0.0;  /* no -0.0 */
    PyObject *k;
    if (info->center == NULL) {
        PyObject *factory = algpy_import_attribute(
            "algpy.algebras.general.clifford._clifford1", "clifford1");
        if (factory == NULL) {
            return NULL;
        }
        info->center = PyObject_CallFunction(factory, "d", s);
        Py_DECREF(factory);
        if (info->center == NULL) {
            return NULL;
        }
    }
    k = Py_BuildValue("[d]", s);
    if (k == NULL || PyObject_SetAttrString(info->center, "K", k) < 0) {
        Py_XDECREF(k);
        return NULL;
    }
    Py_DECREF(k);
    return Py_NewRef(info->center);
}

static PyObject *
pycl5_center_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    clifford1 p[5];
    PyObject *center, *result;
    if (info == NULL || (center = center_type(info)) == NULL) {
        return NULL;
    }
    clifford5_center_minpoly(info->k, ((PyCL5 *)self)->v, p);
    result = algpy_array_of_type(center, p, 5, sizeof(clifford1));
    Py_DECREF(center);
    return result;
}

static PyObject *
pycl5_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {9};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    clifford5_minpoly(info->k, ((PyCL5 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl5_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {8};
    complex ev[8];
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
    clifford5_eivals(info->k, ((PyCL5 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pycl5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    clifford5 v = ((PyCL5 *)self)->v;
    if (info == NULL) return NULL;
    return Py_BuildValue("(O(ddddddddddddddddddddddddddddddddddddd))", rebuild_func,
        info->k.data[0], info->k.data[1], info->k.data[2], info->k.data[3], info->k.data[4],
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pycl5_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CLIFFORD5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford5_basis((clifford5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl5_gens(PyObject *cls, PyObject *Py_UNUSED(args))
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
    clifford5_gens((clifford5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl5_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CLIFFORD5_N, ALGPY_CLIFFORD5_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford5_coeffs(
        info->k, (clifford5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl5_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = clifford5_norm(info->k, ((PyCL5 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycl5_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_inv(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_square(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_exp(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_expm1(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_log(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_log1p(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_sqrt(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_sin(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_cos(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_tan(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_asin(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_acos(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_atan(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_sinh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_cosh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_tanh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_asinh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_acosh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_atanh(info->k, ((PyCL5 *)self)->v));
}

static PyObject *
pycl5_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl5_from_value(Py_TYPE(self), clifford5_re(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl5_from_value(Py_TYPE(self), clifford5_im(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_real(PyObject *self, void *closure)
{
    return pycl5_re(self, NULL);
}

static PyObject *
pycl5_imag(PyObject *self, void *closure)
{
    return pycl5_im(self, NULL);
}

static PyObject *pycl5_grades(PyObject *self, void *closure);

static PyGetSetDef pycl5_getset[] = {
    {"grade", pycl5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycl5_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycl5_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycl5_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(clifford5_is_real(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(clifford5_to_re(((PyCL5 *)self)->v));
}

static PyObject *
pycl5_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cl5_from_value((PyTypeObject *)cls, clifford5_from_re(a));
}

static PyObject *
pycl5_lpow(PyObject *self, PyObject *other)
{
    clifford5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cl5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self),
        clifford5_lpow(info->k, ((PyCL5 *)self)->v, b));
}

static PyObject *
pycl5_rpow(PyObject *self, PyObject *other)
{
    clifford5 b;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    r = cl5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self),
        clifford5_rpow(info->k, ((PyCL5 *)self)->v, b));
}

static PyObject *
pycl5_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    clifford5 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cl5_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self),
        clifford5_mix(info->k, ((PyCL5 *)self)->v, b, t));
}

static PyObject *
pycl5_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pycl5_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 32);
}

static PyObject *
pycl5_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl5_from_value(Py_TYPE(self), clifford5_ipow(info->k, ((PyCL5 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycl5_grade(PyObject *self, Py_ssize_t k)
{
    clifford5 v = ((PyCL5 *)self)->v;
    switch (k) {
        case 0:
            return cl5_from_value(Py_TYPE(self), clifford5_grade0(v));
        case 1:
            return cl5_from_value(Py_TYPE(self), clifford5_grade1(v));
        case 2:
            return cl5_from_value(Py_TYPE(self), clifford5_grade2(v));
        case 3:
            return cl5_from_value(Py_TYPE(self), clifford5_grade3(v));
        case 4:
            return cl5_from_value(Py_TYPE(self), clifford5_grade4(v));
        case 5:
            return cl5_from_value(Py_TYPE(self), clifford5_grade5(v));
        default:
            return cl5_from_value(Py_TYPE(self), clifford5_zero());
    }
}

static PyObject *
pycl5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycl5_grade, 6);
}

static PyMethodDef pycl5_methods[] = {
    {"to_float", (PyCFunction)pycl5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycl5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycl5_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycl5_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycl5_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycl5_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"lpow", pycl5_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycl5_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycl5_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycl5_inv, METH_NOARGS, "inv."},
    {"square", pycl5_square, METH_NOARGS, "square."},
    {"exp", pycl5_exp, METH_NOARGS, "exp."},
    {"expm1", pycl5_expm1, METH_NOARGS, "expm1."},
    {"log", pycl5_log, METH_NOARGS, "log."},
    {"log1p", pycl5_log1p, METH_NOARGS, "log1p."},
    {"sqrt", pycl5_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pycl5_sin, METH_NOARGS, "sin."},
    {"cos", pycl5_cos, METH_NOARGS, "cos."},
    {"tan", pycl5_tan, METH_NOARGS, "tan."},
    {"asin", pycl5_asin, METH_NOARGS, "asin."},
    {"acos", pycl5_acos, METH_NOARGS, "acos."},
    {"atan", pycl5_atan, METH_NOARGS, "atan."},
    {"sinh", pycl5_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycl5_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycl5_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycl5_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycl5_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycl5_atanh, METH_NOARGS, "atanh."},
    {"re", pycl5_re, METH_NOARGS, "re."},
    {"im", pycl5_im, METH_NOARGS, "im."},
    {"arcsin", pycl5_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycl5_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycl5_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycl5_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycl5_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycl5_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycl5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycl5_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycl5_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycl5_gens, METH_NOARGS | METH_CLASS, "Algebra generators in parameter order; excludes the identity."},
    {"coeffs", pycl5_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycl5_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycl5_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pycl5_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycl5_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycl5_minpoly, METH_NOARGS, "Coefficients of the generic minimal polynomial."},
    {"center_minpoly", pycl5_center_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial over the centre, as an array of Type.center values in ascending order."},
    {"eivals", pycl5_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycl5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycl5_length(PyObject *self)
{
    return 32;
}

static PyObject *
pycl5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCL5 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycl5_subscript(PyObject *self, PyObject *key)
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
    return pycl5_item(self, i);
}

static int
pycl5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycl5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycl5_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    clifford5 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "clifford5 does not support modular power");
        return NULL;
    }
    if ((r = cl5_from_object(info, a, &x)) <= 0 ||
            (r = cl5_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cl5_from_value((PyTypeObject *)info->type, clifford5_rpow(info->k, x, y));
}

static PyType_Slot pycl5_slots[] = {
    {Py_tp_dealloc, pycl5_dealloc},
    {Py_tp_traverse, pycl5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycl5_repr},
    {Py_tp_str, pycl5_repr},
    {Py_tp_hash, pycl5_hash},
    {Py_tp_doc, (void *)"Clifford number with float64 components"},
    {Py_tp_richcompare, pycl5_richcompare},
    {Py_tp_methods, pycl5_methods},
    {Py_tp_getset, pycl5_getset},
    {Py_sq_length, pycl5_length},
    {Py_sq_item, pycl5_item},
    {Py_mp_subscript, pycl5_subscript},
    {Py_tp_new, pycl5_new},
    {Py_nb_add, pycl5_add},
    {Py_nb_subtract, pycl5_sub},
    {Py_nb_multiply, pycl5_mul},
    {Py_nb_true_divide, pycl5_div},
    {Py_nb_power, (void *)pycl5_pow},
    {Py_nb_positive, pycl5_pos},
    {Py_nb_negative, pycl5_neg},
    {Py_nb_bool, pycl5_bool},
    {0, NULL}
};

static PyType_Spec pycl5_spec = {
    .name = "algpy.algebras.clifford5",
    .basicsize = sizeof(PyCL5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycl5_slots,
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
npycl5_getitem(void *data, void *arr)
{
    clifford5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return cl5_from_value(descr->typeobj, v);
}

static int
npycl5_setitem(PyObject *item, void *data, void *arr)
{
    clifford5 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cl5_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to clifford5");
        double d[32];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "clifford5 needs exactly 32 components");
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
        v = clifford5_make(d);
    }
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
npycl5_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(clifford5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycl5_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(clifford5));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycl5_nonzero(void *data, void *arr)
{
    clifford5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycl5_arrfuncs;

typedef struct { char c; clifford5 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CL5_CAST_FROM(name, type, convert)                                          \
    static void                                                            \
    npycast_##name##_to_cl5(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        clifford5 *dst = (clifford5 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = clifford5_from_re(convert(src[i])); \
        }                                                                  \
    }

#define CL5_REAL_TYPES(X) \
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

CL5_REAL_TYPES(CL5_CAST_FROM)

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
                npycast_##name##_to_cl5) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CL5_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec5 *k;
    clifford5 (*unary)(clifford5);
    clifford5 (*unary_k)(vec5, clifford5);
    clifford5 (*binary)(vec5, clifford5, clifford5);
    int (*pred)(clifford5, clifford5);
} Op;

static clifford5 op_add(vec5 k, clifford5 a, clifford5 b)
{ return clifford5_add(a, b); }
static clifford5 op_sub(vec5 k, clifford5 a, clifford5 b)
{ return clifford5_sub(a, b); }
static int op_equal(clifford5 a, clifford5 b)
{ return clifford5_equal(a, b); }
static int op_not_equal(clifford5 a, clifford5 b)
{ return !clifford5_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        clifford5 a, r;
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
        clifford5 a, r;
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
        clifford5 a, b, r;
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
        clifford5 a, b;
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
        BIN("multiply", clifford5_mul) ||
        BIN("divide", clifford5_div) ||
        BIN("power", clifford5_rpow) ||
        UN("positive", clifford5_pos) ||
        UN("negative", clifford5_neg) ||
        UN("conjugate", clifford5_conjugate) ||
        UNK("reciprocal", clifford5_inv) ||
        UNK("square", clifford5_square) ||
        UNK("exp", clifford5_exp) ||
        UNK("expm1", clifford5_expm1) ||
        UNK("log", clifford5_log) ||
        UNK("log1p", clifford5_log1p) ||
        UNK("sqrt", clifford5_sqrt) ||
        UNK("sin", clifford5_sin) ||
        UNK("cos", clifford5_cos) ||
        UNK("tan", clifford5_tan) ||
        UNK("arcsin", clifford5_asin) ||
        UNK("arccos", clifford5_acos) ||
        UNK("arctan", clifford5_atan) ||
        UNK("sinh", clifford5_sinh) ||
        UNK("cosh", clifford5_cosh) ||
        UNK("tanh", clifford5_tanh) ||
        UNK("arcsinh", clifford5_asinh) ||
        UNK("arccosh", clifford5_acosh) ||
        UNK("arctanh", clifford5_atanh) ||
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
    if (info->center != NULL) {
        /* Keep the centre type's K = [w*w] in step with these parameters. */
        PyObject *center = center_type(info);
        if (center == NULL) {
            return -1;
        }
        Py_DECREF(center);
    }
    return 0;
}

static PyObject *
meta_get_center(PyObject *self, void *closure)
{
    Info *info = info_of_type((PyTypeObject *)self);
    return info == NULL ? NULL : center_type(info);
}

static PyGetSetDef meta_getset[] = {
    {"center", meta_get_center, NULL,
     "clifford1 type of the centre span{1, w}, w the pseudoscalar, with K = [w*w].", NULL},
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 5 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.clifford5.type",
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

    type = PyType_FromSpecWithBases(&pycl5_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.clifford5.info", NULL);
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
    proto->elsize = sizeof(clifford5);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycl5_arrfuncs;
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
cl5_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4;
    if (!PyArg_ParseTuple(args, "ddddd", &k0, &k1, &k2, &k3, &k4)) return NULL;
    return create_type(k0, k1, k2, k3, k4);
}

static PyObject *
cl5_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3, k4, c[32];
    PyObject *type, *result;
    if (!PyArg_ParseTuple(args, "ddddddddddddddddddddddddddddddddddddd", &k0, &k1, &k2, &k3, &k4,
        c + 0, c + 1, c + 2, c + 3, c + 4, c + 5, c + 6, c + 7, c + 8, c + 9, c + 10, c + 11, c + 12, c + 13, c + 14, c + 15, c + 16, c + 17, c + 18, c + 19, c + 20, c + 21, c + 22, c + 23, c + 24, c + 25, c + 26, c + 27, c + 28, c + 29, c + 30, c + 31)) return NULL;
    type = create_type(k0, k1, k2, k3, k4);
    if (type == NULL) return NULL;
    result = cl5_from_value((PyTypeObject *)type, clifford5_make(c));
    Py_DECREF(type);
    return result;
}

static PyMethodDef module_methods[] = {
    {"clifford5", cl5_factory, METH_VARARGS,
     "clifford5(k0, k1, k2, k3, k4): scalar type / dtype with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2, e3*e3 = k3, e4*e4 = k4."},
    {"_rebuild", cl5_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_clifford5",
    "Clifford number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__clifford5(void)
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
    PyArray_InitArrFuncs(&npycl5_arrfuncs);
    npycl5_arrfuncs.getitem = npycl5_getitem;
    npycl5_arrfuncs.setitem = npycl5_setitem;
    npycl5_arrfuncs.copyswapn = npycl5_copyswapn;
    npycl5_arrfuncs.copyswap = npycl5_copyswap;
    npycl5_arrfuncs.nonzero = npycl5_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
