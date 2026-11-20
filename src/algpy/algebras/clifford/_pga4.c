/*
 * pga4 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in pga4.h.
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

#include "pga4.h"

typedef struct {
    PyObject_HEAD
    pga4 v;
} PyPga4;

static PyObject *PyPga4_Type = NULL;
static int npy_pga4 = -1;

static PyObject *
PyPga4_FromValue(PyTypeObject *type, pga4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyPga4 *p = (PyPga4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_PGA4(v) \
    PyPga4_FromValue((PyTypeObject *)PyPga4_Type, (v))

static int
pga4_from_object(PyObject *o, pga4 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyPga4_Type)) {
        *out = ((PyPga4 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = pga4_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pypga4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyPga4_FromValue(type, pga4_make(d));
}

static PyObject *
pypga4_repr(PyObject *self)
{
    pga4 v = ((PyPga4 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("pga4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pypga4_hash(PyObject *self)
{
    pga4 v = ((PyPga4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_PGA4_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (pga4_is_real(v)) {
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

static PyObject *
pypga4_richcompare(PyObject *a, PyObject *b, int op)
{
    pga4 x, y;
    int r;
    if ((r = pga4_from_object(a, &x)) <= 0 ||
            (r = pga4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(pga4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!pga4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define PGA4_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pypga4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        pga4 x, y;                                                       \
        int r;                                                               \
        if ((r = pga4_from_object(a, &x)) <= 0 ||                        \
                (r = pga4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_PGA4(expr(x, y));                                     \
    }

PGA4_BINARY(add, pga4_add)
PGA4_BINARY(sub, pga4_sub)
PGA4_BINARY(mul, pga4_mul)
PGA4_BINARY(div, pga4_div)
PGA4_BINARY(join, pga4_join)
PGA4_BINARY(meet, pga4_meet)

static PyObject *
pypga4_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    pga4 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "pga4 does not support modular power");
        return NULL;
    }
    if ((r = pga4_from_object(a, &x)) <= 0 ||
            (r = pga4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_PGA4(pga4_rpow(x, y));
}

#define PGA4_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pypga4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_PGA4(expr(((PyPga4 *)self)->v));                  \
    }

PGA4_UNARY(pos, pga4_pos)
PGA4_UNARY(neg, pga4_neg)
PGA4_UNARY(invert, pga4_dual)

static PyObject *
pypga4_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue((PyTypeObject *)cls, pga4_zero());
}

static PyObject *
pypga4_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue((PyTypeObject *)cls, pga4_one());
}

static PyObject *
pypga4_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue((PyTypeObject *)cls, pga4_bottom());
}

static PyObject *
pypga4_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue((PyTypeObject *)cls, pga4_top());
}

static PyObject *
pypga4_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_ldual(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_rdual(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_dual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_dual(((PyPga4 *)self)->v));
}

static int
pypga4_bool(PyObject *self)
{
    pga4 v = ((PyPga4 *)self)->v;
    return vec32_any_nonzero(v);
}

#define PGA4_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pypga4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_PGA4(op(((PyPga4 *)self)->v));                   \
    }

PGA4_METHOD(conjugate, pga4_conjugate)
PGA4_METHOD(exp, pga4_exp)
PGA4_METHOD(log, pga4_log)
PGA4_METHOD(sqrt, pga4_sqrt)

static PyObject *
pypga4_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga4_norm(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_PGA4(pga4_inv(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga4_lmat(((PyPga4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga4_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga4_rmat(((PyPga4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

/* The type of the centre span{1, w}, w the pseudoscalar with w*w = 0:
 * dual1, imported on first use so this extension loads on its own. */
static PyObject *pga4_center_type = NULL;

static PyObject *
pypga4_center(void)
{
    if (pga4_center_type == NULL) {
        pga4_center_type = algpy_import_attribute("algpy.algebras.multicomplex._dual1", "dual1");
        if (pga4_center_type == NULL) {
            return NULL;
        }
    }
    return Py_NewRef(pga4_center_type);
}

static PyObject *
pypga4_center_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 p[5];
    PyObject *center = pypga4_center(), *result;
    if (center == NULL) {
        return NULL;
    }
    pga4_center_minpoly(((PyPga4 *)self)->v, p);
    result = algpy_array_of_type(center, p, 5, sizeof(dual1));
    Py_DECREF(center);
    return result;
}

/* Metatype providing the class attribute pga4.center. */
static PyObject *
pypga4_meta_center(PyObject *Py_UNUSED(self), void *Py_UNUSED(closure))
{
    return pypga4_center();
}

static PyGetSetDef pypga4_meta_getset[] = {
    {"center", pypga4_meta_center, NULL, "Type of the centre span{1, w}: dual1.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PyPga4_MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.pga4.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = pypga4_meta_getset,
};

static PyObject *
pypga4_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {9};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga4_minpoly(((PyPga4 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga4_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[8];
    npy_intp dims[1] = {8};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    pga4_eivals(((PyPga4 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pypga4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga4);
    pga4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga4 *)PyArray_DATA(arr);
    pga4_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga4_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga4);
    pga4 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga4 *)PyArray_DATA(arr);
    pga4_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga4_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_pga4);
    pga4 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (pga4 *)PyArray_DATA(arr);
    pga4_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pypga4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    pga4 v = ((PyPga4 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pypga4_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_square(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_expm1(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_log1p(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_sin(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_cos(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_tan(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_asin(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_acos(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_atan(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_sinh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_cosh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_tanh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_asinh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_acosh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_atanh(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_re(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga4_FromValue(Py_TYPE(self), pga4_im(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_real(PyObject *self, void *closure)
{
    return pypga4_re(self, NULL);
}

static PyObject *
pypga4_imag(PyObject *self, void *closure)
{
    return pypga4_im(self, NULL);
}

static PyObject *pypga4_grades(PyObject *self, void *closure);

static PyGetSetDef pypga4_getset[] = {
    {"grade", pypga4_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pypga4_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pypga4_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pypga4_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(pga4_is_real(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga4_to_re(((PyPga4 *)self)->v));
}

static PyObject *
pypga4_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyPga4_FromValue((PyTypeObject *)cls, pga4_from_re(a));
}

static PyObject *
pypga4_lpow(PyObject *self, PyObject *other)
{
    pga4 b;
    int r;
    r = pga4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga4_FromValue(Py_TYPE(self),
        pga4_lpow(((PyPga4 *)self)->v, b));
}

static PyObject *
pypga4_rpow(PyObject *self, PyObject *other)
{
    pga4 b;
    int r;
    r = pga4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga4_FromValue(Py_TYPE(self),
        pga4_rpow(((PyPga4 *)self)->v, b));
}

static PyObject *
pypga4_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    pga4 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = pga4_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga4_FromValue(Py_TYPE(self),
        pga4_mix(((PyPga4 *)self)->v, b, t));
}

static PyObject *
pypga4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_pga4, 32);
}

static PyObject *
pypga4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_pga4, 32);
}

static PyObject *
pypga4_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_PGA4(pga4_ipow(((PyPga4 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pypga4_grade(PyObject *self, Py_ssize_t k)
{
    pga4 v = ((PyPga4 *)self)->v;
    switch (k) {
        case 0:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade0(v));
        case 1:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade1(v));
        case 2:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade2(v));
        case 3:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade3(v));
        case 4:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade4(v));
        case 5:
            return PyPga4_FromValue(Py_TYPE(self), pga4_grade5(v));
        default:
            return PyPga4_FromValue(Py_TYPE(self), pga4_zero());
    }
}

static PyObject *
pypga4_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pypga4_grade, 6);
}

static PyMethodDef pypga4_methods[] = {
    {"to_float", (PyCFunction)pypga4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pypga4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pypga4_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pypga4_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pypga4_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pypga4_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pypga4_square, METH_NOARGS, "square."},
    {"expm1", pypga4_expm1, METH_NOARGS, "expm1."},
    {"log1p", pypga4_log1p, METH_NOARGS, "log1p."},
    {"sin", pypga4_sin, METH_NOARGS, "sin."},
    {"cos", pypga4_cos, METH_NOARGS, "cos."},
    {"tan", pypga4_tan, METH_NOARGS, "tan."},
    {"asin", pypga4_asin, METH_NOARGS, "asin."},
    {"acos", pypga4_acos, METH_NOARGS, "acos."},
    {"atan", pypga4_atan, METH_NOARGS, "atan."},
    {"sinh", pypga4_sinh, METH_NOARGS, "sinh."},
    {"cosh", pypga4_cosh, METH_NOARGS, "cosh."},
    {"tanh", pypga4_tanh, METH_NOARGS, "tanh."},
    {"asinh", pypga4_asinh, METH_NOARGS, "asinh."},
    {"acosh", pypga4_acosh, METH_NOARGS, "acosh."},
    {"atanh", pypga4_atanh, METH_NOARGS, "atanh."},
    {"re", pypga4_re, METH_NOARGS, "re."},
    {"im", pypga4_im, METH_NOARGS, "im."},
    {"arcsin", pypga4_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pypga4_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pypga4_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pypga4_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pypga4_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pypga4_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pypga4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pypga4_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pypga4_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pypga4_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pypga4_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pypga4_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"dual", pypga4_dual, METH_NOARGS, "Poincare dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"basis", pypga4_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pypga4_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pypga4_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pypga4_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pypga4_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pypga4_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pypga4_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pypga4_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pypga4_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pypga4_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pypga4_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pypga4_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"center_minpoly", pypga4_center_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial over the centre, as an array of dual1 values in ascending order."},
    {"eivals", pypga4_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pypga4_exp, METH_NOARGS, "exp."},
    {"log", pypga4_log, METH_NOARGS, "log."},
    {"sqrt", pypga4_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pypga4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pypga4_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pypga4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyPga4 *)self)->v.data[i]);
}

static PyObject *
pypga4_subscript(PyObject *self, PyObject *key)
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
    return pypga4_item(self, i);
}

static int
pypga4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pypga4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pypga4_slots[] = {
    {Py_tp_dealloc, pypga4_dealloc},
    {Py_tp_traverse, pypga4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pypga4_repr},
    {Py_tp_str, pypga4_repr},
    {Py_tp_hash, pypga4_hash},
    {Py_tp_doc, (void *)"Projective geometric algebra in 4D with 32 float64 coefficients"},
    {Py_tp_richcompare, pypga4_richcompare},
    {Py_tp_methods, pypga4_methods},
    {Py_tp_getset, pypga4_getset},
    {Py_sq_length, pypga4_length},
    {Py_sq_item, pypga4_item},
    {Py_mp_subscript, pypga4_subscript},
    {Py_tp_new, pypga4_new},
    {Py_nb_add, pypga4_add},
    {Py_nb_subtract, pypga4_sub},
    {Py_nb_multiply, pypga4_mul},
    {Py_nb_true_divide, pypga4_div},
    {Py_nb_power, (void *)pypga4_pow},
    {Py_nb_positive, pypga4_pos},
    {Py_nb_negative, pypga4_neg},
    {Py_nb_or, pypga4_join},
    {Py_nb_and, pypga4_meet},
    {Py_nb_invert, pypga4_invert},
    {Py_nb_bool, pypga4_bool},
    {0, NULL}
};

static PyType_Spec pypga4_spec = {
    .name = "algpy.algebras.clifford.pga4",
    .basicsize = sizeof(PyPga4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pypga4_slots,
};

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
npypga4_getitem(void *data, void *arr)
{
    pga4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_PGA4(v);
}

static int
npypga4_setitem(PyObject *item, void *data, void *arr)
{
    pga4 v;
    int r = pga4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to pga4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "pga4 needs exactly 32 components");
            return -1;
        }
        for (int i = 0; i < 32; i++) {
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
pypga4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(pga4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pypga4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pypga4_copyswapn(dst, sizeof(pga4), src, sizeof(pga4), 1,
                         swap, arr);
}

static npy_bool
pypga4_nonzero(void *data, void *arr)
{
    pga4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pypga4_arrfuncs;
typedef struct { char c; pga4 v; } pypga4_align_test;

static PyArray_DescrProto pypga4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(pga4), offsetof(pypga4_align_test, v), 0, 0, 0,
    &pypga4_arrfuncs,
};

#define PGA4_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pypga4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        pga4 *dst = (pga4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = pga4_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define PGA4_REAL_TYPES(X) \
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

PGA4_REAL_TYPES(PGA4_CAST_FROM)

static int
pypga4_register_casts(void)
{
#define REGISTER_PGA4_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_pga4,                    \
                                     pypga4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_pga4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    PGA4_REAL_TYPES(REGISTER_PGA4_CAST)
#undef REGISTER_PGA4_CAST
    return 0;
}

typedef pga4 (*pga4_unary_op)(pga4);
typedef pga4 (*pga4_binary_op)(pga4, pga4);
typedef int (*pga4_predicate)(pga4, pga4);

static void
loop_pga4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    pga4_unary_op op = (pga4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_pga4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    pga4_binary_op op = (pga4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga4 a, b, r;
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
loop_pga4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    pga4_predicate op = (pga4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pypga4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_pga4;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_pga4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pypga4_register_ufuncs(PyObject *numpy)
{
    return
        pypga4_register_loop(numpy, "add", loop_pga4_binary, 3, 0, pga4_add) ||
        pypga4_register_loop(numpy, "subtract", loop_pga4_binary, 3, 0, pga4_sub) ||
        pypga4_register_loop(numpy, "multiply", loop_pga4_binary, 3, 0, pga4_mul) ||
        pypga4_register_loop(numpy, "divide", loop_pga4_binary, 3, 0, pga4_div) ||
        pypga4_register_loop(numpy, "power", loop_pga4_binary, 3, 0, pga4_rpow) ||
        pypga4_register_loop(numpy, "positive", loop_pga4_unary, 2, 0, pga4_pos) ||
        pypga4_register_loop(numpy, "negative", loop_pga4_unary, 2, 0, pga4_neg) ||
        pypga4_register_loop(numpy, "conjugate", loop_pga4_unary, 2, 0, pga4_conjugate) ||
        pypga4_register_loop(numpy, "reciprocal", loop_pga4_unary, 2, 0, pga4_inv) ||
        pypga4_register_loop(numpy, "square", loop_pga4_unary, 2, 0, pga4_square) ||
        pypga4_register_loop(numpy, "exp", loop_pga4_unary, 2, 0, pga4_exp) ||
        pypga4_register_loop(numpy, "expm1", loop_pga4_unary, 2, 0, pga4_expm1) ||
        pypga4_register_loop(numpy, "log", loop_pga4_unary, 2, 0, pga4_log) ||
        pypga4_register_loop(numpy, "log1p", loop_pga4_unary, 2, 0, pga4_log1p) ||
        pypga4_register_loop(numpy, "sqrt", loop_pga4_unary, 2, 0, pga4_sqrt) ||
        pypga4_register_loop(numpy, "sin", loop_pga4_unary, 2, 0, pga4_sin) ||
        pypga4_register_loop(numpy, "cos", loop_pga4_unary, 2, 0, pga4_cos) ||
        pypga4_register_loop(numpy, "tan", loop_pga4_unary, 2, 0, pga4_tan) ||
        pypga4_register_loop(numpy, "arcsin", loop_pga4_unary, 2, 0, pga4_asin) ||
        pypga4_register_loop(numpy, "arccos", loop_pga4_unary, 2, 0, pga4_acos) ||
        pypga4_register_loop(numpy, "arctan", loop_pga4_unary, 2, 0, pga4_atan) ||
        pypga4_register_loop(numpy, "sinh", loop_pga4_unary, 2, 0, pga4_sinh) ||
        pypga4_register_loop(numpy, "cosh", loop_pga4_unary, 2, 0, pga4_cosh) ||
        pypga4_register_loop(numpy, "tanh", loop_pga4_unary, 2, 0, pga4_tanh) ||
        pypga4_register_loop(numpy, "arcsinh", loop_pga4_unary, 2, 0, pga4_asinh) ||
        pypga4_register_loop(numpy, "arccosh", loop_pga4_unary, 2, 0, pga4_acosh) ||
        pypga4_register_loop(numpy, "arctanh", loop_pga4_unary, 2, 0, pga4_atanh) ||
        pypga4_register_loop(numpy, "equal", loop_pga4_predicate, 3, 1, pga4_equal) ||
        pypga4_register_loop(numpy, "bitwise_or", loop_pga4_binary, 3, 0, pga4_join) ||
        pypga4_register_loop(numpy, "bitwise_and", loop_pga4_binary, 3, 0, pga4_meet) ||
        pypga4_register_loop(numpy, "invert", loop_pga4_unary, 2, 0, pga4_dual) ||
        pypga4_register_loop(numpy, "not_equal", loop_pga4_predicate, 3, 1, pga4_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pypga4_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pypga4_register_loop(ufuncs, "ldual", loop_pga4_unary, 2, 0, pga4_ldual) ||
        pypga4_register_loop(ufuncs, "rdual", loop_pga4_unary, 2, 0, pga4_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_pga4", "pga4 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__pga4(void)
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
    PyPga4_MetaType.tp_base = &PyType_Type;  /* not a constant initializer on MSVC */
    if (PyType_Ready(&PyPga4_MetaType) < 0) {
        goto fail;
    }
    PyPga4_Type = PyType_FromModuleAndSpec(
        m, &pypga4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyPga4_Type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(PyPga4_Type, &PyPga4_MetaType);
    PyArray_InitArrFuncs(&pypga4_arrfuncs);
    pypga4_arrfuncs.getitem = npypga4_getitem;
    pypga4_arrfuncs.setitem = npypga4_setitem;
    pypga4_arrfuncs.copyswapn = pypga4_copyswapn;
    pypga4_arrfuncs.copyswap = pypga4_copyswap;
    pypga4_arrfuncs.nonzero = pypga4_nonzero;
    pypga4_descr_proto.typeobj = (PyTypeObject *)PyPga4_Type;
    Py_SET_TYPE((PyObject *)&pypga4_descr_proto, &PyArrayDescr_Type);
    npy_pga4 = PyArray_RegisterDataType(&pypga4_descr_proto);
    if (npy_pga4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_pga4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyPga4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyPga4_Type, 32) < 0) {
        goto fail;
    }
    if (pypga4_register_casts() < 0 ||
            pypga4_register_ufuncs(numpy) < 0 ||
            pypga4_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "pga4", PyPga4_Type) < 0) {
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
