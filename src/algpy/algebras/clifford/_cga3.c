/*
 * cga3 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in cga3.h.
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

#include "cga3.h"

typedef struct {
    PyObject_HEAD
    cga3 v;
} PyCga3;

static PyObject *PyCga3_Type = NULL;
static int npy_cga3 = -1;

static PyObject *
PyCga3_FromValue(PyTypeObject *type, cga3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCga3 *p = (PyCga3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_CGA3(v) \
    PyCga3_FromValue((PyTypeObject *)PyCga3_Type, (v))

static int
cga3_from_object(PyObject *o, cga3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyCga3_Type)) {
        *out = ((PyCga3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cga3_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycga3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyCga3_FromValue(type, cga3_make(d));
}

static PyObject *
pycga3_repr(PyObject *self)
{
    cga3 v = ((PyCga3 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cga3(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycga3_hash(PyObject *self)
{
    cga3 v = ((PyCga3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CGA3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cga3_is_real(v)) {
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
pycga3_richcompare(PyObject *a, PyObject *b, int op)
{
    cga3 x, y;
    int r;
    if ((r = cga3_from_object(a, &x)) <= 0 ||
            (r = cga3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(cga3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!cga3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define CGA3_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycga3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        cga3 x, y;                                                       \
        int r;                                                               \
        if ((r = cga3_from_object(a, &x)) <= 0 ||                        \
                (r = cga3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_CGA3(expr(x, y));                                     \
    }

CGA3_BINARY(add, cga3_add)
CGA3_BINARY(sub, cga3_sub)
CGA3_BINARY(mul, cga3_mul)
CGA3_BINARY(div, cga3_div)
CGA3_BINARY(join, cga3_join)
CGA3_BINARY(meet, cga3_meet)

static PyObject *
pycga3_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cga3 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cga3 does not support modular power");
        return NULL;
    }
    if ((r = cga3_from_object(a, &x)) <= 0 ||
            (r = cga3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_CGA3(cga3_rpow(x, y));
}

#define CGA3_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycga3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_CGA3(expr(((PyCga3 *)self)->v));                  \
    }

CGA3_UNARY(pos, cga3_pos)
CGA3_UNARY(neg, cga3_neg)
CGA3_UNARY(invert, cga3_dual)

static PyObject *
pycga3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue((PyTypeObject *)cls, cga3_zero());
}

static PyObject *
pycga3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue((PyTypeObject *)cls, cga3_one());
}

static PyObject *
pycga3_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue((PyTypeObject *)cls, cga3_bottom());
}

static PyObject *
pycga3_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue((PyTypeObject *)cls, cga3_top());
}

static PyObject *
pycga3_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_ldual(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_rdual(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_dual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_dual(((PyCga3 *)self)->v));
}

static int
pycga3_bool(PyObject *self)
{
    cga3 v = ((PyCga3 *)self)->v;
    return vec32_any_nonzero(v);
}

#define CGA3_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycga3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_CGA3(op(((PyCga3 *)self)->v));                   \
    }

CGA3_METHOD(conjugate, cga3_conjugate)
CGA3_METHOD(exp, cga3_exp)
CGA3_METHOD(log, cga3_log)
CGA3_METHOD(sqrt, cga3_sqrt)

static PyObject *
pycga3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cga3_norm(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_CGA3(cga3_inv(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga3_lmat(((PyCga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga3_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga3_rmat(((PyCga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

/* The type of the centre span{1, w}, w the pseudoscalar with w*w = -1:
 * complex1, imported on first use so this extension loads on its own. */
static PyObject *cga3_center_type = NULL;

static PyObject *
pycga3_center(void)
{
    if (cga3_center_type == NULL) {
        cga3_center_type = algpy_import_attribute("algpy.algebras.multicomplex._complex1", "complex1");
        if (cga3_center_type == NULL) {
            return NULL;
        }
    }
    return Py_NewRef(cga3_center_type);
}

static PyObject *
pycga3_center_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex1 p[5];
    PyObject *center = pycga3_center(), *result;
    if (center == NULL) {
        return NULL;
    }
    cga3_center_minpoly(((PyCga3 *)self)->v, p);
    result = algpy_array_of_type(center, p, 5, sizeof(complex1));
    Py_DECREF(center);
    return result;
}

/* Metatype providing the class attribute cga3.center. */
static PyObject *
pycga3_meta_center(PyObject *Py_UNUSED(self), void *Py_UNUSED(closure))
{
    return pycga3_center();
}

static PyGetSetDef pycga3_meta_getset[] = {
    {"center", pycga3_meta_center, NULL, "Type of the centre span{1, w}: complex1.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PyCga3_MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.cga3.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = pycga3_meta_getset,
};

static PyObject *
pycga3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {9};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga3_minpoly(((PyCga3 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[8];
    npy_intp dims[1] = {8};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    cga3_eivals(((PyCga3 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycga3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga3);
    cga3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga3 *)PyArray_DATA(arr);
    cga3_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga3_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga3);
    cga3 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga3 *)PyArray_DATA(arr);
    cga3_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga3_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {32, 32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_cga3);
    cga3 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (cga3 *)PyArray_DATA(arr);
    cga3_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycga3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    cga3 v = ((PyCga3 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pycga3_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_square(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_expm1(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_log1p(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_sin(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_cos(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_tan(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_asin(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_acos(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_atan(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_sinh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_cosh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_tanh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_asinh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_acosh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_atanh(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_re(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga3_FromValue(Py_TYPE(self), cga3_im(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_real(PyObject *self, void *closure)
{
    return pycga3_re(self, NULL);
}

static PyObject *
pycga3_imag(PyObject *self, void *closure)
{
    return pycga3_im(self, NULL);
}

static PyObject *pycga3_grades(PyObject *self, void *closure);

static PyGetSetDef pycga3_getset[] = {
    {"grade", pycga3_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycga3_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycga3_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycga3_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cga3_is_real(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cga3_to_re(((PyCga3 *)self)->v));
}

static PyObject *
pycga3_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyCga3_FromValue((PyTypeObject *)cls, cga3_from_re(a));
}

static PyObject *
pycga3_lpow(PyObject *self, PyObject *other)
{
    cga3 b;
    int r;
    r = cga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga3_FromValue(Py_TYPE(self),
        cga3_lpow(((PyCga3 *)self)->v, b));
}

static PyObject *
pycga3_rpow(PyObject *self, PyObject *other)
{
    cga3 b;
    int r;
    r = cga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga3_FromValue(Py_TYPE(self),
        cga3_rpow(((PyCga3 *)self)->v, b));
}

static PyObject *
pycga3_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cga3 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cga3_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga3_FromValue(Py_TYPE(self),
        cga3_mix(((PyCga3 *)self)->v, b, t));
}

static PyObject *
pycga3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_cga3, 32);
}

static PyObject *
pycga3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_cga3, 32);
}

static PyObject *
pycga3_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_CGA3(cga3_ipow(((PyCga3 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycga3_grade(PyObject *self, Py_ssize_t k)
{
    cga3 v = ((PyCga3 *)self)->v;
    switch (k) {
        case 0:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade0(v));
        case 1:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade1(v));
        case 2:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade2(v));
        case 3:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade3(v));
        case 4:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade4(v));
        case 5:
            return PyCga3_FromValue(Py_TYPE(self), cga3_grade5(v));
        default:
            return PyCga3_FromValue(Py_TYPE(self), cga3_zero());
    }
}

static PyObject *
pycga3_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycga3_grade, 6);
}

static PyMethodDef pycga3_methods[] = {
    {"to_float", (PyCFunction)pycga3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycga3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycga3_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycga3_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycga3_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycga3_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycga3_square, METH_NOARGS, "square."},
    {"expm1", pycga3_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycga3_log1p, METH_NOARGS, "log1p."},
    {"sin", pycga3_sin, METH_NOARGS, "sin."},
    {"cos", pycga3_cos, METH_NOARGS, "cos."},
    {"tan", pycga3_tan, METH_NOARGS, "tan."},
    {"asin", pycga3_asin, METH_NOARGS, "asin."},
    {"acos", pycga3_acos, METH_NOARGS, "acos."},
    {"atan", pycga3_atan, METH_NOARGS, "atan."},
    {"sinh", pycga3_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycga3_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycga3_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycga3_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycga3_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycga3_atanh, METH_NOARGS, "atanh."},
    {"re", pycga3_re, METH_NOARGS, "re."},
    {"im", pycga3_im, METH_NOARGS, "im."},
    {"arcsin", pycga3_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycga3_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycga3_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycga3_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycga3_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycga3_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycga3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycga3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pycga3_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pycga3_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pycga3_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pycga3_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"dual", pycga3_dual, METH_NOARGS, "Poincare dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"basis", pycga3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycga3_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pycga3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycga3_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pycga3_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycga3_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycga3_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycga3_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycga3_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pycga3_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycga3_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycga3_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"center_minpoly", pycga3_center_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial over the centre, as an array of complex1 values in ascending order."},
    {"eivals", pycga3_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycga3_exp, METH_NOARGS, "exp."},
    {"log", pycga3_log, METH_NOARGS, "log."},
    {"sqrt", pycga3_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycga3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycga3_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pycga3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCga3 *)self)->v.data[i]);
}

static PyObject *
pycga3_subscript(PyObject *self, PyObject *key)
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
    return pycga3_item(self, i);
}

static int
pycga3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycga3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycga3_slots[] = {
    {Py_tp_dealloc, pycga3_dealloc},
    {Py_tp_traverse, pycga3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycga3_repr},
    {Py_tp_str, pycga3_repr},
    {Py_tp_hash, pycga3_hash},
    {Py_tp_doc, (void *)"Conformal geometric algebra in 3D with 32 float64 coefficients"},
    {Py_tp_richcompare, pycga3_richcompare},
    {Py_tp_methods, pycga3_methods},
    {Py_tp_getset, pycga3_getset},
    {Py_sq_length, pycga3_length},
    {Py_sq_item, pycga3_item},
    {Py_mp_subscript, pycga3_subscript},
    {Py_tp_new, pycga3_new},
    {Py_nb_add, pycga3_add},
    {Py_nb_subtract, pycga3_sub},
    {Py_nb_multiply, pycga3_mul},
    {Py_nb_true_divide, pycga3_div},
    {Py_nb_power, (void *)pycga3_pow},
    {Py_nb_positive, pycga3_pos},
    {Py_nb_negative, pycga3_neg},
    {Py_nb_or, pycga3_join},
    {Py_nb_and, pycga3_meet},
    {Py_nb_invert, pycga3_invert},
    {Py_nb_bool, pycga3_bool},
    {0, NULL}
};

static PyType_Spec pycga3_spec = {
    .name = "algpy.algebras.clifford.cga3",
    .basicsize = sizeof(PyCga3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycga3_slots,
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
npycga3_getitem(void *data, void *arr)
{
    cga3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_CGA3(v);
}

static int
npycga3_setitem(PyObject *item, void *data, void *arr)
{
    cga3 v;
    int r = cga3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cga3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cga3 needs exactly 32 components");
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
pycga3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cga3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycga3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycga3_copyswapn(dst, sizeof(cga3), src, sizeof(cga3), 1,
                         swap, arr);
}

static npy_bool
pycga3_nonzero(void *data, void *arr)
{
    cga3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycga3_arrfuncs;
typedef struct { char c; cga3 v; } pycga3_align_test;

static PyArray_DescrProto pycga3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(cga3), offsetof(pycga3_align_test, v), 0, 0, 0,
    &pycga3_arrfuncs,
};

#define CGA3_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycga3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        cga3 *dst = (cga3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = cga3_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define CGA3_REAL_TYPES(X) \
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

CGA3_REAL_TYPES(CGA3_CAST_FROM)

static int
pycga3_register_casts(void)
{
#define REGISTER_CGA3_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_cga3,                    \
                                     pycga3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_cga3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    CGA3_REAL_TYPES(REGISTER_CGA3_CAST)
#undef REGISTER_CGA3_CAST
    return 0;
}

typedef cga3 (*cga3_unary_op)(cga3);
typedef cga3 (*cga3_binary_op)(cga3, cga3);
typedef int (*cga3_predicate)(cga3, cga3);

static void
loop_cga3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    cga3_unary_op op = (cga3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_cga3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    cga3_binary_op op = (cga3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga3 a, b, r;
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
loop_cga3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    cga3_predicate op = (cga3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycga3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_cga3;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_cga3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycga3_register_ufuncs(PyObject *numpy)
{
    return
        pycga3_register_loop(numpy, "add", loop_cga3_binary, 3, 0, cga3_add) ||
        pycga3_register_loop(numpy, "subtract", loop_cga3_binary, 3, 0, cga3_sub) ||
        pycga3_register_loop(numpy, "multiply", loop_cga3_binary, 3, 0, cga3_mul) ||
        pycga3_register_loop(numpy, "divide", loop_cga3_binary, 3, 0, cga3_div) ||
        pycga3_register_loop(numpy, "power", loop_cga3_binary, 3, 0, cga3_rpow) ||
        pycga3_register_loop(numpy, "positive", loop_cga3_unary, 2, 0, cga3_pos) ||
        pycga3_register_loop(numpy, "negative", loop_cga3_unary, 2, 0, cga3_neg) ||
        pycga3_register_loop(numpy, "conjugate", loop_cga3_unary, 2, 0, cga3_conjugate) ||
        pycga3_register_loop(numpy, "reciprocal", loop_cga3_unary, 2, 0, cga3_inv) ||
        pycga3_register_loop(numpy, "square", loop_cga3_unary, 2, 0, cga3_square) ||
        pycga3_register_loop(numpy, "exp", loop_cga3_unary, 2, 0, cga3_exp) ||
        pycga3_register_loop(numpy, "expm1", loop_cga3_unary, 2, 0, cga3_expm1) ||
        pycga3_register_loop(numpy, "log", loop_cga3_unary, 2, 0, cga3_log) ||
        pycga3_register_loop(numpy, "log1p", loop_cga3_unary, 2, 0, cga3_log1p) ||
        pycga3_register_loop(numpy, "sqrt", loop_cga3_unary, 2, 0, cga3_sqrt) ||
        pycga3_register_loop(numpy, "sin", loop_cga3_unary, 2, 0, cga3_sin) ||
        pycga3_register_loop(numpy, "cos", loop_cga3_unary, 2, 0, cga3_cos) ||
        pycga3_register_loop(numpy, "tan", loop_cga3_unary, 2, 0, cga3_tan) ||
        pycga3_register_loop(numpy, "arcsin", loop_cga3_unary, 2, 0, cga3_asin) ||
        pycga3_register_loop(numpy, "arccos", loop_cga3_unary, 2, 0, cga3_acos) ||
        pycga3_register_loop(numpy, "arctan", loop_cga3_unary, 2, 0, cga3_atan) ||
        pycga3_register_loop(numpy, "sinh", loop_cga3_unary, 2, 0, cga3_sinh) ||
        pycga3_register_loop(numpy, "cosh", loop_cga3_unary, 2, 0, cga3_cosh) ||
        pycga3_register_loop(numpy, "tanh", loop_cga3_unary, 2, 0, cga3_tanh) ||
        pycga3_register_loop(numpy, "arcsinh", loop_cga3_unary, 2, 0, cga3_asinh) ||
        pycga3_register_loop(numpy, "arccosh", loop_cga3_unary, 2, 0, cga3_acosh) ||
        pycga3_register_loop(numpy, "arctanh", loop_cga3_unary, 2, 0, cga3_atanh) ||
        pycga3_register_loop(numpy, "equal", loop_cga3_predicate, 3, 1, cga3_equal) ||
        pycga3_register_loop(numpy, "bitwise_or", loop_cga3_binary, 3, 0, cga3_join) ||
        pycga3_register_loop(numpy, "bitwise_and", loop_cga3_binary, 3, 0, cga3_meet) ||
        pycga3_register_loop(numpy, "invert", loop_cga3_unary, 2, 0, cga3_dual) ||
        pycga3_register_loop(numpy, "not_equal", loop_cga3_predicate, 3, 1, cga3_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pycga3_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pycga3_register_loop(ufuncs, "ldual", loop_cga3_unary, 2, 0, cga3_ldual) ||
        pycga3_register_loop(ufuncs, "rdual", loop_cga3_unary, 2, 0, cga3_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cga3", "cga3 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cga3(void)
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
    PyCga3_MetaType.tp_base = &PyType_Type;  /* not a constant initializer on MSVC */
    if (PyType_Ready(&PyCga3_MetaType) < 0) {
        goto fail;
    }
    PyCga3_Type = PyType_FromModuleAndSpec(
        m, &pycga3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyCga3_Type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(PyCga3_Type, &PyCga3_MetaType);
    PyArray_InitArrFuncs(&pycga3_arrfuncs);
    pycga3_arrfuncs.getitem = npycga3_getitem;
    pycga3_arrfuncs.setitem = npycga3_setitem;
    pycga3_arrfuncs.copyswapn = pycga3_copyswapn;
    pycga3_arrfuncs.copyswap = pycga3_copyswap;
    pycga3_arrfuncs.nonzero = pycga3_nonzero;
    pycga3_descr_proto.typeobj = (PyTypeObject *)PyCga3_Type;
    Py_SET_TYPE((PyObject *)&pycga3_descr_proto, &PyArrayDescr_Type);
    npy_cga3 = PyArray_RegisterDataType(&pycga3_descr_proto);
    if (npy_cga3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_cga3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyCga3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyCga3_Type, 32) < 0) {
        goto fail;
    }
    if (pycga3_register_casts() < 0 ||
            pycga3_register_ufuncs(numpy) < 0 ||
            pycga3_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "cga3", PyCga3_Type) < 0) {
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
