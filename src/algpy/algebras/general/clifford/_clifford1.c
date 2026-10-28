/*
 * Multicomplex numbers x + y*e (e*e = k) exposed to Python as NumPy
 * user-defined dtypes.
 *
 * The parameter k is not part of the data: `clifford1(k)` creates a new scalar type and NumPy dtype
 * on every call (not cached; keep and reuse it). The parameters live in the
 * mutable class attribute `K` (a list). The arithmetic
 * itself lives in `clifford1.h`.
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

#include "clifford1.h"

typedef struct {
    PyObject *type;
    int typenum;
    vec1 k;
} Info;

typedef struct {
    PyObject_HEAD
    clifford1 v;
} PyCL1;

static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.clifford1.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
cl1_from_value(PyTypeObject *type, clifford1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCL1 *p = (PyCL1 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a clifford1 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
cl1_from_object(Info *info, PyObject *o, clifford1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyCL1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = clifford1_make(d, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycl1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "y", NULL};
    double x = 0, y = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist, &x, &y)) {
        return NULL;
    }
    return cl1_from_value(type, clifford1_make(x, y));
}

static PyObject *
pycl1_repr(PyObject *self)
{
    clifford1 v = ((PyCL1 *)self)->v;
    char *s[2] = {NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 2; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("clifford1" "(%s, %s)", s[0], s[1]);
done:
    for (i = 0; i < 2; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycl1_hash(PyObject *self)
{
    clifford1 v = ((PyCL1 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    if (v.data[1] == 0) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dd)", v.data[0], v.data[1]);
    }
    if (obj == NULL) {
        return -1;
    }
    h = PyObject_Hash(obj);
    Py_DECREF(obj);
    return h;
}

/* Resolve the type info of whichever operand is a clifford1. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyObject *o = PyType_GetSlot(Py_TYPE(a), Py_tp_new) == (void *)pycl1_new
                  ? a : b;
    return info_of_type(Py_TYPE(o));
}

static PyObject *
pycl1_richcompare(PyObject *a, PyObject *b, int op)
{
    clifford1 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = cl1_from_object(info, a, &x)) <= 0 ||
            (r = cl1_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(clifford1_equal(x, y) == (op == Py_EQ));
}

#define CL1_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pycl1_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        clifford1 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = cl1_from_object(info, a, &x)) <= 0 ||                     \
                (r = cl1_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return cl1_from_value((PyTypeObject *)info->type, expr);           \
    }

CL1_BINARY(add, clifford1_add(x, y))
CL1_BINARY(sub, clifford1_sub(x, y))
CL1_BINARY(mul, clifford1_mul(info->k, x, y))
CL1_BINARY(div, clifford1_div(info->k, x, y))

static PyObject *
pycl1_pos(PyObject *self)
{
    return cl1_from_value(Py_TYPE(self), clifford1_pos(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_neg(PyObject *self)
{
    return cl1_from_value(Py_TYPE(self), clifford1_neg(((PyCL1 *)self)->v));
}

static int
pycl1_bool(PyObject *self)
{
    return vec2_any_nonzero(((PyCL1 *)self)->v);
}

static PyObject *
pycl1_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl1_from_value(Py_TYPE(self),
                          clifford1_conjugate(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return PyFloat_FromDouble(
        clifford1_norm2(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl1_from_value((PyTypeObject *)cls, clifford1_zero());
}

static PyObject *
pycl1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return cl1_from_value((PyTypeObject *)cls, clifford1_make(1, 0));
}

static PyObject *
pycl1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {2, 2};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    clifford1_lmat(info->k, ((PyCL1 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[2] = {2, 2};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    clifford1_rmat(info->k, ((PyCL1 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford1_minpoly(info->k, ((PyCL1 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycl1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
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
    clifford1_eivals(info->k, ((PyCL1 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *rebuild_func = NULL;

static PyObject *
pycl1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    clifford1 v = ((PyCL1 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(ddd))", rebuild_func, info->k.data[0],
                         v.data[0], v.data[1]);
}

static PyObject *
pycl1_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_CLIFFORD1_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford1_basis((clifford1 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl1_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_CLIFFORD1_N, ALGPY_CLIFFORD1_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    clifford1_coeffs(
        info->k, (clifford1 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pycl1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    complex n = clifford1_norm(info->k, ((PyCL1 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycl1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_inv(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_square(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_exp(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_expm1(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_log(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_log1p(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_sqrt(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_sqrt(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_sin(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_cos(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_tan(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_asin(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_acos(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_atan(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_sinh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_cosh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_tanh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_asinh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_acosh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self), clifford1_atanh(info->k, ((PyCL1 *)self)->v));
}

static PyObject *
pycl1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl1_from_value(Py_TYPE(self), clifford1_re(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return cl1_from_value(Py_TYPE(self), clifford1_im(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_real(PyObject *self, void *closure)
{
    return pycl1_re(self, NULL);
}

static PyObject *
pycl1_imag(PyObject *self, void *closure)
{
    return pycl1_im(self, NULL);
}

static PyGetSetDef pycl1_getset[] = {
    {"real", pycl1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycl1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycl1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(clifford1_is_real(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(clifford1_to_re(((PyCL1 *)self)->v));
}

static PyObject *
pycl1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return cl1_from_value((PyTypeObject *)cls, clifford1_from_re(a));
}

static PyObject *
pycl1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    clifford1 b;
    double t;
    int r;
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cl1_from_object(info, other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return cl1_from_value(Py_TYPE(self),
        clifford1_mix(info->k, ((PyCL1 *)self)->v, b, t));
}

static PyObject *
pycl1_to_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_to_float(args, kwargs, info->typenum, 2);
}

static PyObject *
pycl1_from_float(PyObject *cls, PyObject *args, PyObject *kwargs)
{
    Info *info = info_of_type((PyTypeObject *)cls);
    if (info == NULL) {
        return NULL;
    }
    return algpy_from_float(args, kwargs, info->typenum, 2);
}

static PyMethodDef pycl1_methods[] = {
    {"to_float", (PyCFunction)pycl1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycl1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycl1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycl1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycl1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycl1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"inv", pycl1_inv, METH_NOARGS, "inv. See the C header for TODO status."},
    {"square", pycl1_square, METH_NOARGS, "square. See the C header for TODO status."},
    {"exp", pycl1_exp, METH_NOARGS, "exp. See the C header for TODO status."},
    {"expm1", pycl1_expm1, METH_NOARGS, "expm1. See the C header for TODO status."},
    {"log", pycl1_log, METH_NOARGS, "log. See the C header for TODO status."},
    {"log1p", pycl1_log1p, METH_NOARGS, "log1p. See the C header for TODO status."},
    {"sqrt", pycl1_sqrt, METH_NOARGS, "sqrt. See the C header for TODO status."},
    {"sin", pycl1_sin, METH_NOARGS, "sin. See the C header for TODO status."},
    {"cos", pycl1_cos, METH_NOARGS, "cos. See the C header for TODO status."},
    {"tan", pycl1_tan, METH_NOARGS, "tan. See the C header for TODO status."},
    {"asin", pycl1_asin, METH_NOARGS, "asin. See the C header for TODO status."},
    {"acos", pycl1_acos, METH_NOARGS, "acos. See the C header for TODO status."},
    {"atan", pycl1_atan, METH_NOARGS, "atan. See the C header for TODO status."},
    {"sinh", pycl1_sinh, METH_NOARGS, "sinh. See the C header for TODO status."},
    {"cosh", pycl1_cosh, METH_NOARGS, "cosh. See the C header for TODO status."},
    {"tanh", pycl1_tanh, METH_NOARGS, "tanh. See the C header for TODO status."},
    {"asinh", pycl1_asinh, METH_NOARGS, "asinh. See the C header for TODO status."},
    {"acosh", pycl1_acosh, METH_NOARGS, "acosh. See the C header for TODO status."},
    {"atanh", pycl1_atanh, METH_NOARGS, "atanh. See the C header for TODO status."},
    {"re", pycl1_re, METH_NOARGS, "re. See the C header for TODO status."},
    {"im", pycl1_im, METH_NOARGS, "im. See the C header for TODO status."},
    {"arcsin", pycl1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycl1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycl1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycl1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycl1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycl1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycl1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycl1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycl1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycl1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycl1_norm, METH_NOARGS, "Norm as a Python complex."},
    {"conjugate", pycl1_conjugate, METH_NOARGS, "Conjugate."},
    {"norm2", pycl1_norm2, METH_NOARGS, "Norm x*x - k*y*y."},
    {"lmat", pycl1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycl1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycl1_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial."},
    {"eivals", pycl1_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pycl1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycl1_length(PyObject *self)
{
    return 2;
}

static PyObject *
pycl1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCL1 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pycl1_subscript(PyObject *self, PyObject *key)
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
        i += 2;
    }
    return pycl1_item(self, i);
}

static int
pycl1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycl1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyObject *
pycl1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    clifford1 x, y;
    int r;
    Info *info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "clifford1 does not support modular power");
        return NULL;
    }
    if ((r = cl1_from_object(info, a, &x)) <= 0 ||
            (r = cl1_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return cl1_from_value((PyTypeObject *)info->type, clifford1_pow(info->k, x, y));
}

static PyType_Slot pycl1_slots[] = {
    {Py_tp_dealloc, pycl1_dealloc},
    {Py_tp_traverse, pycl1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycl1_repr},
    {Py_tp_str, pycl1_repr},
    {Py_tp_hash, pycl1_hash},
    {Py_tp_doc, (void *)"One-dimensional algebra number with float64 components"},
    {Py_tp_richcompare, pycl1_richcompare},
    {Py_tp_methods, pycl1_methods},
    {Py_tp_getset, pycl1_getset},
    {Py_sq_length, pycl1_length},
    {Py_sq_item, pycl1_item},
    {Py_mp_subscript, pycl1_subscript},
    {Py_tp_new, pycl1_new},
    {Py_nb_add, pycl1_add},
    {Py_nb_subtract, pycl1_sub},
    {Py_nb_multiply, pycl1_mul},
    {Py_nb_true_divide, pycl1_div},
    {Py_nb_power, (void *)pycl1_pow},
    {Py_nb_positive, pycl1_pos},
    {Py_nb_negative, pycl1_neg},
    {Py_nb_bool, pycl1_bool},
    {0, NULL}
};

static PyType_Spec pycl1_spec = {
    .name = "algpy.algebras.general.clifford.clifford1",
    .basicsize = sizeof(PyCL1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycl1_slots,
};

/* NumPy dtype support */

static PyObject *
npycl1_getitem(void *data, void *arr)
{
    clifford1 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    return cl1_from_value(descr->typeobj, v);
}

static int
npycl1_setitem(PyObject *item, void *data, void *arr)
{
    clifford1 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = cl1_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to clifford1");
        double d[2];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "clifford1 needs exactly 2 components");
            return -1;
        }
        for (i = 0; i < 2; i++) {
            d[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (d[i] == -1.0 && PyErr_Occurred()) {
                Py_DECREF(seq);
                return -1;
            }
        }
        Py_DECREF(seq);
        v = clifford1_make(d[0], d[1]);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static void
npycl1_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(clifford1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npycl1_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(clifford1));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npycl1_nonzero(void *data, void *arr)
{
    clifford1 v;
    memcpy(&v, data, sizeof(v));
    return vec2_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npycl1_arrfuncs;

typedef struct { char c; clifford1 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define CL1_CAST_FROM(name, type)                                          \
    static void                                                            \
    npycast_##name##_to_cl1(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        clifford1 *dst = (clifford1 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = clifford1_make((double)src[i], 0);                \
        }                                                                  \
    }

#define CL1_REAL_TYPES(X)                                                  \
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

CL1_REAL_TYPES(CL1_CAST_FROM)

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
                npycast_##name##_to_cl1) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    CL1_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec1 *k;
    clifford1 (*unary)(clifford1);
    clifford1 (*unary_k)(vec1, clifford1);
    clifford1 (*binary)(vec1, clifford1, clifford1);
    int (*pred)(clifford1, clifford1);
} Op;

static clifford1 op_add(vec1 k, clifford1 a, clifford1 b)
{ return clifford1_add(a, b); }
static clifford1 op_sub(vec1 k, clifford1 a, clifford1 b)
{ return clifford1_sub(a, b); }
static int op_equal(clifford1 a, clifford1 b)
{ return clifford1_equal(a, b); }
static int op_not_equal(clifford1 a, clifford1 b)
{ return !clifford1_equal(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        clifford1 a, r;
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
        clifford1 a, r;
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
        clifford1 a, b, r;
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
        clifford1 a, b;
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
        BIN("multiply", clifford1_mul) ||
        BIN("divide", clifford1_div) ||
        BIN("power", clifford1_pow) ||
        UN("positive", clifford1_pos) ||
        UN("negative", clifford1_neg) ||
        UN("conjugate", clifford1_conjugate) ||
        UNK("reciprocal", clifford1_inv) ||
        UNK("square", clifford1_square) ||
        UNK("exp", clifford1_exp) ||
        UNK("expm1", clifford1_expm1) ||
        UNK("log", clifford1_log) ||
        UNK("log1p", clifford1_log1p) ||
        UNK("sqrt", clifford1_sqrt) ||
        UNK("sin", clifford1_sin) ||
        UNK("cos", clifford1_cos) ||
        UNK("tan", clifford1_tan) ||
        UNK("arcsin", clifford1_asin) ||
        UNK("arccos", clifford1_acos) ||
        UNK("arctan", clifford1_atan) ||
        UNK("sinh", clifford1_sinh) ||
        UNK("cosh", clifford1_cosh) ||
        UNK("tanh", clifford1_tanh) ||
        UNK("arcsinh", clifford1_asinh) ||
        UNK("arccosh", clifford1_acosh) ||
        UNK("arctanh", clifford1_atanh) ||
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
    list = PyList_New(1);
    for (i = 0; list != NULL && i < 1; i++) {
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
    vec1 k;
    int i;
    if (info == NULL) {
        return -1;
    }
    if (value == NULL) {
        PyErr_SetString(PyExc_AttributeError, "cannot delete K");
        return -1;
    }
    seq = PySequence_Fast(value, "K must be a sequence of 1 numbers");
    if (seq == NULL) {
        return -1;
    }
    if (PySequence_Fast_GET_SIZE(seq) != 1) {
        PyErr_SetString(PyExc_ValueError, "K must have exactly 1 element(s)");
        Py_DECREF(seq);
        return -1;
    }
    for (i = 0; i < 1; i++) {
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
    {"K", meta_get_K, meta_set_K, "Mutable parameters (list of 1 floats).",
     NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.clifford1.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = meta_getset,
};

/* Type factory */

static PyObject *
create_type(double k)
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
    info->k.data[0] = k;

    type = PyType_FromSpecWithBases(&pycl1_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.clifford1.info", NULL);
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
    proto->elsize = sizeof(clifford1);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npycl1_arrfuncs;
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

    if (algpy_set_dimension(type, 2) < 0) {
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
cl1_factory(PyObject *module, PyObject *args)
{
    double k;
    if (!PyArg_ParseTuple(args, "d", &k)) {
        return NULL;
    }
    return create_type(k);
}

static PyObject *
cl1_rebuild(PyObject *module, PyObject *args)
{
    double k, x, y;
    PyObject *type, *res;
    if (!PyArg_ParseTuple(args, "ddd", &k, &x, &y)) {
        return NULL;
    }
    type = create_type(k);
    if (type == NULL) {
        return NULL;
    }
    res = cl1_from_value((PyTypeObject *)type, clifford1_make(x, y));
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"clifford1", cl1_factory, METH_VARARGS,
     "Factory for one-dimensional algebra scalar types / NumPy dtypes."},
    {"_rebuild", cl1_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_clifford1",
    "One-dimensional algebra number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__clifford1(void)
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
    PyArray_InitArrFuncs(&npycl1_arrfuncs);
    npycl1_arrfuncs.getitem = npycl1_getitem;
    npycl1_arrfuncs.setitem = npycl1_setitem;
    npycl1_arrfuncs.copyswapn = npycl1_copyswapn;
    npycl1_arrfuncs.copyswap = npycl1_copyswap;
    npycl1_arrfuncs.nonzero = npycl1_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
