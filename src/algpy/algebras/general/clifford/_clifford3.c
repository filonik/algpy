/*
 * Clifford numbers x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 +
 * x12*e1*e2 + x012*e0*e1*e2 (e0*e0 = k0, e1*e1 = k1, e2*e2 = k2) exposed to
 * Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `clifford3(k0, k1, k2)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `clifford3.h`.
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

#include "clifford3.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec3 k;
} Info;

typedef struct {
    PyObject_HEAD
    clifford3 v;
} PyCL3;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.clifford3.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cl3_from_value(PyTypeObject *type, clifford3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCL3 *p = (PyCL3 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a clifford3 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cl3_from_object(Info *info, PyObject *o, clifford3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCL3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = clifford3_make(d, 0, 0, 0, 0, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycl3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12",
                             "x012", NULL};
    double c[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
                                     c, c + 1, c + 2, c + 3, c + 4, c + 5,
                                     c + 6, c + 7)) {
        return NULL;
    }
    return cl3_from_value(type, clifford3_make(c[0], c[1], c[2], c[3],
                                                   c[4], c[5], c[6], c[7]));
}

static PyObject *
pycl3_repr(PyObject *self)
{
    clifford3 v = ((PyCL3 *)self)->v;
    char *s[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("clifford3(%s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6],
                                  s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycl3_hash(PyObject *self)
{
    clifford3 v = ((PyCL3 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    if (clifford3_is_real(v)) {
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

/* Resolve the type info of whichever operand is a clifford3. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyObject *o = PyType_GetSlot(Py_TYPE(a), Py_tp_new) == (void *)pycl3_new
                  ? a : b;
    return info_of_type(Py_TYPE(o));
}

static PyObject *
pycl3_richcompare(PyObject *a, PyObject *b, int op)
{
    clifford3 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cl3_from_object(info, a, &x)) <= 0 ||
            (r = cl3_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(clifford3_equal(x, y) == (op == Py_EQ));
}

#define CL3_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycl3_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        clifford3 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cl3_from_object(info, a, &x)) <= 0 ||                     \
                (r = cl3_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cl3_from_value((PyTypeObject *)info->type, expr);           \
    }

CL3_BINARY(add, clifford3_add(x, y))
CL3_BINARY(sub, clifford3_sub(x, y))
CL3_BINARY(mul, clifford3_mul(info->k, x, y))
CL3_BINARY(div, clifford3_div(info->k, x, y))

static PyObject *
pycl3_pos(PyObject *self)
{
    return cl3_from_value(Py_TYPE(self), clifford3_pos(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_neg(PyObject *self)
{
    return cl3_from_value(Py_TYPE(self), clifford3_neg(((PyCL3 *)self)->v));
}

static int
pycl3_bool(PyObject *self)
{
    return vec8_any_nonzero(((PyCL3 *)self)->v);
}

static PyObject *
pycl3_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl3_from_value(Py_TYPE(self),
                          clifford3_conjugate(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl3_from_value((PyTypeObject *)cls, clifford3_zero());
}

static PyObject *
pycl3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl3_from_value((PyTypeObject *)cls, clifford3_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycl3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford3_lmat(info->k, ((PyCL3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford3_rmat(info->k, ((PyCL3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {5};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    clifford3_minpoly(info->k, ((PyCL3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {4};
    complex ev[4];
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
    clifford3_eivals(info->k, ((PyCL3 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pycl3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    clifford3 v = ((PyCL3 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(ddddddddddd))", rebuild_func, info->k.data[0],
                         info->k.data[1], info->k.data[2], v.data[0],
                         v.data[1], v.data[2], v.data[3], v.data[4],
                         v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pycl3_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CLIFFORD3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford3_basis((clifford3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl3_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CLIFFORD3_N, ALGPY_CLIFFORD3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford3_coeffs(
        info->k, (clifford3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = clifford3_norm(info->k, ((PyCL3 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycl3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_inv(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_square(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_exp(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_expm1(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_log(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_log1p(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_sqrt(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_sin(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_cos(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_tan(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_asin(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_acos(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_atan(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_sinh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_cosh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_tanh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_asinh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_acosh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self), clifford3_atanh(info->k, ((PyCL3 *)self)->v));
}

static PyObject *
pycl3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl3_from_value(Py_TYPE(self), clifford3_re(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl3_from_value(Py_TYPE(self), clifford3_im(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_real(PyObject *self, void *closure)
{
    return pycl3_re(self, NULL);
}

static PyObject *
pycl3_imag(PyObject *self, void *closure)
{
    return pycl3_im(self, NULL);
}

static PyGetSetDef pycl3_getset[] = {
    {"real", pycl3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycl3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycl3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(clifford3_is_real(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(clifford3_to_re(((PyCL3 *)self)->v));
}

static PyObject *
pycl3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cl3_from_value((PyTypeObject *)cls, clifford3_from_re(a));
}

static PyObject *
pycl3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    clifford3 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cl3_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cl3_from_value(Py_TYPE(self),
        clifford3_mix(info->k, ((PyCL3 *)self)->v, b, t));
}

static PyObject *
pycl3_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 8);
}

static PyObject *
pycl3_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 8);
}

static PyMethodDef pycl3_methods[] = {
    {"to_float", (PyCFunction)pycl3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycl3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycl3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycl3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycl3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycl3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"inv", pycl3_inv, METH_NOARGS, "inv. See the C header for TODO status."},
    {"square", pycl3_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"exp", pycl3_exp, METH_NOARGS, "exp. See the C header for TODO status."},
    {"expm1", pycl3_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log", pycl3_log, METH_NOARGS, "log. See the C header for TODO status."},
    {"log1p", pycl3_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sqrt", pycl3_sqrt, METH_NOARGS, "sqrt. See the C header for TODO status."},
    {"sin", pycl3_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pycl3_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pycl3_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pycl3_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pycl3_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pycl3_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pycl3_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pycl3_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pycl3_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pycl3_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pycl3_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pycl3_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pycl3_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pycl3_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pycl3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycl3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycl3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycl3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycl3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycl3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycl3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycl3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycl3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycl3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycl3_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycl3_conjugate, METH_NOARGS, "Conjugate."},
    {"lmat", pycl3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycl3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycl3_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial."},
    {"eivals", pycl3_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycl3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycl3_length(PyObject *self)
{
    return 8;
}

static PyObject *
pycl3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCL3 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycl3_subscript(PyObject *self, PyObject *key)
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
    return pycl3_item(self, i);
}

static int
pycl3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycl3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycl3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    clifford3 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "clifford3 does not support modular power");
        return NULL;
    }
    if ((r = cl3_from_object(info, a, &x)) <= 0 ||
            (r = cl3_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cl3_from_value((PyTypeObject *)info->type, clifford3_pow(info->k, x, y));
}

static PyType_Slot pycl3_slots[] = {
    {Py_tp_dealloc, pycl3_dealloc},
    {Py_tp_traverse, pycl3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycl3_repr},
    {Py_tp_str, pycl3_repr},
    {Py_tp_hash, pycl3_hash},
    {Py_tp_doc, (void *)"Clifford number with float64 components"},
    {Py_tp_richcompare, pycl3_richcompare},
    {Py_tp_methods, pycl3_methods},
    {Py_tp_getset, pycl3_getset},
    {Py_sq_length, pycl3_length},
    {Py_sq_item, pycl3_item},
    {Py_mp_subscript, pycl3_subscript},
    {Py_tp_new, pycl3_new},
    {Py_nb_add, pycl3_add},
    {Py_nb_subtract, pycl3_sub},
    {Py_nb_multiply, pycl3_mul},
    {Py_nb_true_divide, pycl3_div},
    {Py_nb_power, (void *)pycl3_pow},
    {Py_nb_positive, pycl3_pos},
    {Py_nb_negative, pycl3_neg},
    {Py_nb_bool, pycl3_bool},
    {0, NULL}
};

static PyType_Spec pycl3_spec = {
    .name = "algpy.algebras.clifford3",
    .basicsize = sizeof(PyCL3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycl3_slots,
};

/* NumPy dtype support */

static PyObject *
npycl3_getitem(void *data, void *arr)
{
    clifford3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    return cl3_from_value(descr->typeobj, v);
}

static int
npycl3_setitem(PyObject *item, void *data, void *arr)
{
    clifford3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cl3_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to clifford3");
        double d[8];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "clifford3 needs exactly 8 components");
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
        v = clifford3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

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

static void
npycl3_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(clifford3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycl3_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(clifford3));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycl3_nonzero(void *data, void *arr)
{
    clifford3 v;
    memcpy(&v, data, sizeof(v));
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycl3_arrfuncs;

typedef struct { char c; clifford3 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CL3_CAST_FROM(name, type)                                          \
    static void                                                            \
    npycast_##name##_to_cl3(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        clifford3 *dst = (clifford3 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = clifford3_make((double)src[i], 0, 0, 0, 0, 0, 0, 0);                \
        }                                                                  \
    }

#define CL3_REAL_TYPES(X)                                                  \
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

CL3_REAL_TYPES(CL3_CAST_FROM)

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
                npycast_##name##_to_cl3) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CL3_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec3 *k;
    clifford3 (*unary)(clifford3);
    clifford3 (*unary_k)(vec3, clifford3);
    clifford3 (*binary)(vec3, clifford3, clifford3);
    int (*pred)(clifford3, clifford3);
} Op;

static clifford3 op_add(vec3 k, clifford3 a, clifford3 b)
{ return clifford3_add(a, b); }
static clifford3 op_sub(vec3 k, clifford3 a, clifford3 b)
{ return clifford3_sub(a, b); }
static int op_equal(clifford3 a, clifford3 b)
{ return clifford3_equal(a, b); }
static int op_not_equal(clifford3 a, clifford3 b)
{ return !clifford3_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        clifford3 a, r;
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
        clifford3 a, r;
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
        clifford3 a, b, r;
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
        clifford3 a, b;
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
        BIN("multiply", clifford3_mul) ||
        BIN("divide", clifford3_div) ||
        BIN("power", clifford3_pow) ||
        UN("positive", clifford3_pos) ||
        UN("negative", clifford3_neg) ||
        UN("conjugate", clifford3_conjugate) ||
        UNK("reciprocal", clifford3_inv) ||
        UNK("square", clifford3_square) ||
        UNK("exp", clifford3_exp) ||
        UNK("expm1", clifford3_expm1) ||
        UNK("log", clifford3_log) ||
        UNK("log1p", clifford3_log1p) ||
        UNK("sqrt", clifford3_sqrt) ||
        UNK("sin", clifford3_sin) ||
        UNK("cos", clifford3_cos) ||
        UNK("tan", clifford3_tan) ||
        UNK("arcsin", clifford3_asin) ||
        UNK("arccos", clifford3_acos) ||
        UNK("arctan", clifford3_atan) ||
        UNK("sinh", clifford3_sinh) ||
        UNK("cosh", clifford3_cosh) ||
        UNK("tanh", clifford3_tanh) ||
        UNK("arcsinh", clifford3_asinh) ||
        UNK("arccosh", clifford3_acosh) ||
        UNK("arctanh", clifford3_atanh) ||
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
    .tp_name = "algpy.clifford3.type",
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

    type = PyType_FromSpecWithBases(&pycl3_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.clifford3.info", NULL);
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
    proto->elsize = sizeof(clifford3);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycl3_arrfuncs;
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
cl3_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2;
    if (!PyArg_ParseTuple(args, "ddd", &k0, &k1, &k2)) {
        return NULL;
    }
    return create_type(k0, k1, k2);
}

static PyObject *
cl3_rebuild(PyObject *module, PyObject *args)
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
    res = cl3_from_value((PyTypeObject *)type, clifford3_make(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]));
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"clifford3", cl3_factory, METH_VARARGS,
     "clifford3(k0, k1, k2): scalar type / dtype of numbers with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2."},
    {"_rebuild", cl3_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_clifford3",
    "Clifford number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__clifford3(void)
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
    PyArray_InitArrFuncs(&npycl3_arrfuncs);
    npycl3_arrfuncs.getitem = npycl3_getitem;
    npycl3_arrfuncs.setitem = npycl3_setitem;
    npycl3_arrfuncs.copyswapn = npycl3_copyswapn;
    npycl3_arrfuncs.copyswap = npycl3_copyswap;
    npycl3_arrfuncs.nonzero = npycl3_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
