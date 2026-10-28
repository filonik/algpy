/*
 * Multicomplex numbers sum_S x_S * e_S over subsets S of {0, 1, 2, 3}
 * (component index = bitmask of S; e0*e0 = k0, e1*e1 = k1, e2*e2 = k2,
 * e3*e3 = k3) exposed to Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `multicomplex4(k0, k1, k2, k3)` creates
 * a new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `multicomplex4.h`.
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

#include "multicomplex4.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec4 k;
} Info;

typedef struct {
    PyObject_HEAD
    multicomplex4 v;
} PyMC4;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.multicomplex4.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
mc4_from_value(PyTypeObject *type, multicomplex4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyMC4 *p = (PyMC4 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a multicomplex4 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
mc4_from_object(Info *info, PyObject *o, multicomplex4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyMC4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = multicomplex4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pymc4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12",
                             "x012", "x3", "x03", "x13", "x013", "x23",
                             "x023", "x123", "x0123", NULL};
    multicomplex4 v = {{0}};
    double *c = v.data;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
                                     c, c + 1, c + 2, c + 3, c + 4, c + 5,
                                     c + 6, c + 7, c + 8, c + 9, c + 10,
                                     c + 11, c + 12, c + 13, c + 14,
                                     c + 15)) {
        return NULL;
    }
    return mc4_from_value(type, v);
}

static PyObject *
pymc4_repr(PyObject *self)
{
    multicomplex4 v = ((PyMC4 *)self)->v;
    char *s[16] = {NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("multicomplex4(%s, %s, %s, %s, %s, %s, %s, %s, "
                                  "%s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6],
                                  s[7], s[8], s[9], s[10], s[11], s[12],
                                  s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pymc4_hash(PyObject *self)
{
    multicomplex4 v = ((PyMC4 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    if (multicomplex4_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddd)", v.data[0], v.data[1],
                            v.data[2], v.data[3], v.data[4], v.data[5],
                            v.data[6], v.data[7], v.data[8], v.data[9],
                            v.data[10], v.data[11], v.data[12], v.data[13],
                            v.data[14], v.data[15]);
    }
    if (obj == NULL) {
        return -1;
    }
    h = PyObject_Hash(obj);
    Py_DECREF(obj);
    return h;
}

/* Resolve the type info of whichever operand is a multicomplex4. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyObject *o = PyType_GetSlot(Py_TYPE(a), Py_tp_new) == (void *)pymc4_new
                  ? a : b;
    return info_of_type(Py_TYPE(o));
}

static PyObject *
pymc4_richcompare(PyObject *a, PyObject *b, int op)
{
    multicomplex4 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = mc4_from_object(info, a, &x)) <= 0 ||
            (r = mc4_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(multicomplex4_equal(x, y) == (op == Py_EQ));
}

#define MC4_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pymc4_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        multicomplex4 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = mc4_from_object(info, a, &x)) <= 0 ||                     \
                (r = mc4_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return mc4_from_value((PyTypeObject *)info->type, expr);           \
    }

MC4_BINARY(add, multicomplex4_add(x, y))
MC4_BINARY(sub, multicomplex4_sub(x, y))
MC4_BINARY(mul, multicomplex4_mul(info->k, x, y))
MC4_BINARY(div, multicomplex4_div(info->k, x, y))

static PyObject *
pymc4_pos(PyObject *self)
{
    return mc4_from_value(Py_TYPE(self), multicomplex4_pos(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_neg(PyObject *self)
{
    return mc4_from_value(Py_TYPE(self), multicomplex4_neg(((PyMC4 *)self)->v));
}

static int
pymc4_bool(PyObject *self)
{
    return vec16_any_nonzero(((PyMC4 *)self)->v);
}

static PyObject *
pymc4_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc4_from_value(Py_TYPE(self),
                          multicomplex4_conjugate(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc4_from_value((PyTypeObject *)cls, multicomplex4_zero());
}

static PyObject *
pymc4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc4_from_value((PyTypeObject *)cls, multicomplex4_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pymc4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    multicomplex4_lmat(info->k, ((PyMC4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    multicomplex4_rmat(info->k, ((PyMC4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {17};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex4_minpoly(info->k, ((PyMC4 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {16};
    complex ev[16];
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
    multicomplex4_eivals(info->k, ((PyMC4 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 16; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pymc4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    multicomplex4 v = ((PyMC4 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(dddddddddddddddddddd))", rebuild_func,
                         info->k.data[0], info->k.data[1], info->k.data[2],
                         info->k.data[3], v.data[0], v.data[1], v.data[2],
                         v.data[3], v.data[4], v.data[5], v.data[6],
                         v.data[7], v.data[8], v.data[9], v.data[10],
                         v.data[11], v.data[12], v.data[13], v.data[14],
                         v.data[15]);
}

static PyObject *
pymc4_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_MULTICOMPLEX4_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex4_basis((multicomplex4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc4_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_MULTICOMPLEX4_N, ALGPY_MULTICOMPLEX4_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex4_coeffs(
        info->k, (multicomplex4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = multicomplex4_norm(info->k, ((PyMC4 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pymc4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_inv(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_square(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_exp(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_expm1(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_log(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_log1p(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_sqrt(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_sin(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_cos(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_tan(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_asin(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_acos(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_atan(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_sinh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_cosh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_tanh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_asinh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_acosh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self), multicomplex4_atanh(info->k, ((PyMC4 *)self)->v));
}

static PyObject *
pymc4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc4_from_value(Py_TYPE(self), multicomplex4_re(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc4_from_value(Py_TYPE(self), multicomplex4_im(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_real(PyObject *self, void *closure)
{
    return pymc4_re(self, NULL);
}

static PyObject *
pymc4_imag(PyObject *self, void *closure)
{
    return pymc4_im(self, NULL);
}

static PyGetSetDef pymc4_getset[] = {
    {"real", pymc4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pymc4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pymc4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(multicomplex4_is_real(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(multicomplex4_to_re(((PyMC4 *)self)->v));
}

static PyObject *
pymc4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return mc4_from_value((PyTypeObject *)cls, multicomplex4_from_re(a));
}

static PyObject *
pymc4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    multicomplex4 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = mc4_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return mc4_from_value(Py_TYPE(self),
        multicomplex4_mix(info->k, ((PyMC4 *)self)->v, b, t));
}

static PyObject *
pymc4_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 16);
}

static PyObject *
pymc4_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 16);
}

static PyMethodDef pymc4_methods[] = {
    {"to_float", (PyCFunction)pymc4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pymc4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pymc4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pymc4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pymc4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pymc4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"inv", pymc4_inv, METH_NOARGS, "inv. See the C header for TODO status."},
    {"square", pymc4_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"exp", pymc4_exp, METH_NOARGS, "exp. See the C header for TODO status."},
    {"expm1", pymc4_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log", pymc4_log, METH_NOARGS, "log. See the C header for TODO status."},
    {"log1p", pymc4_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sqrt", pymc4_sqrt, METH_NOARGS, "sqrt. See the C header for TODO status."},
    {"sin", pymc4_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pymc4_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pymc4_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pymc4_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pymc4_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pymc4_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pymc4_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pymc4_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pymc4_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pymc4_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pymc4_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pymc4_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pymc4_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pymc4_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pymc4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pymc4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pymc4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pymc4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pymc4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pymc4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pymc4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pymc4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pymc4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pymc4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pymc4_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pymc4_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pymc4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pymc4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pymc4_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial."},
    {"eivals", pymc4_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pymc4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pymc4_length(PyObject *self)
{
    return 16;
}

static PyObject *
pymc4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyMC4 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pymc4_subscript(PyObject *self, PyObject *key)
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
    return pymc4_item(self, i);
}

static int
pymc4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pymc4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pymc4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    multicomplex4 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "multicomplex4 does not support modular power");
        return NULL;
    }
    if ((r = mc4_from_object(info, a, &x)) <= 0 ||
            (r = mc4_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return mc4_from_value((PyTypeObject *)info->type, multicomplex4_pow(info->k, x, y));
}

static PyType_Slot pymc4_slots[] = {
    {Py_tp_dealloc, pymc4_dealloc},
    {Py_tp_traverse, pymc4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pymc4_repr},
    {Py_tp_str, pymc4_repr},
    {Py_tp_hash, pymc4_hash},
    {Py_tp_doc, (void *)"Multicomplex number with float64 components"},
    {Py_tp_richcompare, pymc4_richcompare},
    {Py_tp_methods, pymc4_methods},
    {Py_tp_getset, pymc4_getset},
    {Py_sq_length, pymc4_length},
    {Py_sq_item, pymc4_item},
    {Py_mp_subscript, pymc4_subscript},
    {Py_tp_new, pymc4_new},
    {Py_nb_add, pymc4_add},
    {Py_nb_subtract, pymc4_sub},
    {Py_nb_multiply, pymc4_mul},
    {Py_nb_true_divide, pymc4_div},
    {Py_nb_power, (void *)pymc4_pow},
    {Py_nb_positive, pymc4_pos},
    {Py_nb_negative, pymc4_neg},
    {Py_nb_bool, pymc4_bool},
    {0, NULL}
};

static PyType_Spec pymc4_spec = {
    .name = "algpy.algebras.multicomplex4",
    .basicsize = sizeof(PyMC4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pymc4_slots,
};

/* NumPy dtype support */

static PyObject *
npymc4_getitem(void *data, void *arr)
{
    multicomplex4 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    return mc4_from_value(descr->typeobj, v);
}

static int
npymc4_setitem(PyObject *item, void *data, void *arr)
{
    multicomplex4 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = mc4_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to multicomplex4");
        double d[16];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "multicomplex4 needs exactly 16 components");
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
        memcpy(v.data, d, sizeof(d));
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

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

static void
npymc4_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(multicomplex4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npymc4_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(multicomplex4));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npymc4_nonzero(void *data, void *arr)
{
    multicomplex4 v;
    memcpy(&v, data, sizeof(v));
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npymc4_arrfuncs;

typedef struct { char c; multicomplex4 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define MC4_CAST_FROM(name, type)                                          \
    static void                                                            \
    npycast_##name##_to_mc4(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        multicomplex4 *dst = (multicomplex4 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = multicomplex4_from_re((double)src[i]);                              \
        }                                                                  \
    }

#define MC4_REAL_TYPES(X)                                                  \
    X(BOOL, npy_bool)                                                      \
    X(BYTE, npy_byte)                                                      \
    X(UBYTE, npy_ubyte)                                                    \
    X(SHORT, npy_short)                                                    \
    X(USHORT, npy_ushort)                                                  \
    X(INT, npy_int)                                                        \
    X(UINT, npy_uint)                                                      \
    X(LONG, npy_long)                                                      \
    X(ULONG, npy_ulong)                                                    \
    X(LONGLONG, npy_longlong)                                              \
    X(ULONGLONG, npy_ulonglong)                                            \
    X(FLOAT, npy_float)                                                    \
    X(DOUBLE, npy_double)

MC4_REAL_TYPES(MC4_CAST_FROM)

static int
register_casts(int typenum)
{
#define REGISTER_CAST(name, type)                                          \
    {                                                                      \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);           \
        if (from == NULL) {                                                \
            return -1;                                                     \
        }                                                                  \
        if (PyArray_RegisterCastFunc(from, typenum,                        \
                npycast_##name##_to_mc4) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    MC4_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec4 *k;
    multicomplex4 (*unary)(multicomplex4);
    multicomplex4 (*unary_k)(vec4, multicomplex4);
    multicomplex4 (*binary)(vec4, multicomplex4, multicomplex4);
    int (*pred)(multicomplex4, multicomplex4);
} Op;

static multicomplex4 op_add(vec4 k, multicomplex4 a, multicomplex4 b)
{ return multicomplex4_add(a, b); }
static multicomplex4 op_sub(vec4 k, multicomplex4 a, multicomplex4 b)
{ return multicomplex4_sub(a, b); }
static int op_equal(multicomplex4 a, multicomplex4 b)
{ return multicomplex4_equal(a, b); }
static int op_not_equal(multicomplex4 a, multicomplex4 b)
{ return !multicomplex4_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        multicomplex4 a, r;
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
        multicomplex4 a, r;
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
        multicomplex4 a, b, r;
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
        multicomplex4 a, b;
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
        BIN("multiply", multicomplex4_mul) ||
        BIN("divide", multicomplex4_div) ||
        BIN("power", multicomplex4_pow) ||
        UN("positive", multicomplex4_pos) ||
        UN("negative", multicomplex4_neg) ||
        UN("conjugate", multicomplex4_conjugate) ||
        UNK("reciprocal", multicomplex4_inv) ||
        UNK("square", multicomplex4_square) ||
        UNK("exp", multicomplex4_exp) ||
        UNK("expm1", multicomplex4_expm1) ||
        UNK("log", multicomplex4_log) ||
        UNK("log1p", multicomplex4_log1p) ||
        UNK("sqrt", multicomplex4_sqrt) ||
        UNK("sin", multicomplex4_sin) ||
        UNK("cos", multicomplex4_cos) ||
        UNK("tan", multicomplex4_tan) ||
        UNK("arcsin", multicomplex4_asin) ||
        UNK("arccos", multicomplex4_acos) ||
        UNK("arctan", multicomplex4_atan) ||
        UNK("sinh", multicomplex4_sinh) ||
        UNK("cosh", multicomplex4_cosh) ||
        UNK("tanh", multicomplex4_tanh) ||
        UNK("arcsinh", multicomplex4_asinh) ||
        UNK("arccosh", multicomplex4_acosh) ||
        UNK("arctanh", multicomplex4_atanh) ||
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
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 3 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.multicomplex4.type",
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

    type = PyType_FromSpecWithBases(&pymc4_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.multicomplex4.info", NULL);
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
    proto->elsize = sizeof(multicomplex4);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npymc4_arrfuncs;
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
mc4_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3;
    if (!PyArg_ParseTuple(args, "dddd", &k0, &k1, &k2, &k3)) {
        return NULL;
    }
    return create_type(k0, k1, k2, k3);
}

static PyObject *
mc4_rebuild(PyObject *module, PyObject *args)
{
    double k0, k1, k2, k3;
    multicomplex4 v;
    double *c = v.data;
    PyObject *type, *res;
    if (!PyArg_ParseTuple(args, "dddddddddddddddddddd", &k0, &k1, &k2, &k3,
                          c, c + 1, c + 2, c + 3, c + 4, c + 5, c + 6, c + 7,
                          c + 8, c + 9, c + 10, c + 11, c + 12, c + 13,
                          c + 14, c + 15)) {
        return NULL;
    }
    type = create_type(k0, k1, k2, k3);
    if (type == NULL) {
        return NULL;
    }
    res = mc4_from_value((PyTypeObject *)type, v);
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"multicomplex4", mc4_factory, METH_VARARGS,
     "multicomplex4(k0, k1, k2, k3): scalar type / dtype of numbers with e0*e0 = k0, ..., e3*e3 = k3."},
    {"_rebuild", mc4_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_multicomplex4",
    "Multicomplex number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__multicomplex4(void)
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
    PyArray_InitArrFuncs(&npymc4_arrfuncs);
    npymc4_arrfuncs.getitem = npymc4_getitem;
    npymc4_arrfuncs.setitem = npymc4_setitem;
    npymc4_arrfuncs.copyswapn = npymc4_copyswapn;
    npymc4_arrfuncs.copyswap = npymc4_copyswap;
    npymc4_arrfuncs.nonzero = npymc4_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
