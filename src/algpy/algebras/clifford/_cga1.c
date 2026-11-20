/*
 * cga1 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in cga1.h.
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

#include "cga1.h"

typedef struct {
    PyObject_HEAD
    cga1 v;
} PyCga1;

static PyObject *PyCga1_Type = NULL;
static int npy_cga1 = -1;

static PyObject *
PyCga1_FromValue(PyTypeObject *type, cga1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyCga1 *p = (PyCga1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_CGA1(v) \
    PyCga1_FromValue((PyTypeObject *)PyCga1_Type, (v))

static int
cga1_from_object(PyObject *o, cga1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyCga1_Type)) {
        *out = ((PyCga1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = cga1_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pycga1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyCga1_FromValue(type, cga1_make(d));
}

static PyObject *
pycga1_repr(PyObject *self)
{
    cga1 v = ((PyCga1 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("cga1(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycga1_hash(PyObject *self)
{
    cga1 v = ((PyCga1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_CGA1_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (cga1_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pycga1_richcompare(PyObject *a, PyObject *b, int op)
{
    cga1 x, y;
    int r;
    if ((r = cga1_from_object(a, &x)) <= 0 ||
            (r = cga1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(cga1_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!cga1_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define CGA1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycga1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        cga1 x, y;                                                       \
        int r;                                                               \
        if ((r = cga1_from_object(a, &x)) <= 0 ||                        \
                (r = cga1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_CGA1(expr(x, y));                                     \
    }

CGA1_BINARY(add, cga1_add)
CGA1_BINARY(sub, cga1_sub)
CGA1_BINARY(mul, cga1_mul)
CGA1_BINARY(div, cga1_div)
CGA1_BINARY(join, cga1_join)
CGA1_BINARY(meet, cga1_meet)

static PyObject *
pycga1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    cga1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "cga1 does not support modular power");
        return NULL;
    }
    if ((r = cga1_from_object(a, &x)) <= 0 ||
            (r = cga1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_CGA1(cga1_rpow(x, y));
}

#define CGA1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycga1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_CGA1(expr(((PyCga1 *)self)->v));                  \
    }

CGA1_UNARY(pos, cga1_pos)
CGA1_UNARY(neg, cga1_neg)
CGA1_UNARY(invert, cga1_dual)

static PyObject *
pycga1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue((PyTypeObject *)cls, cga1_zero());
}

static PyObject *
pycga1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue((PyTypeObject *)cls, cga1_one());
}

static PyObject *
pycga1_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue((PyTypeObject *)cls, cga1_bottom());
}

static PyObject *
pycga1_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue((PyTypeObject *)cls, cga1_top());
}

static PyObject *
pycga1_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_ldual(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_rdual(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_dual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_dual(((PyCga1 *)self)->v));
}

static int
pycga1_bool(PyObject *self)
{
    cga1 v = ((PyCga1 *)self)->v;
    return vec8_any_nonzero(v);
}

#define CGA1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycga1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_CGA1(op(((PyCga1 *)self)->v));                   \
    }

CGA1_METHOD(conjugate, cga1_conjugate)
CGA1_METHOD(exp, cga1_exp)
CGA1_METHOD(log, cga1_log)
CGA1_METHOD(sqrt, cga1_sqrt)

static PyObject *
pycga1_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = cga1_norm(((PyCga1 *)self)->v);
    return PyComplex_FromDoubles(n.data[0], n.data[1]);
}

static PyObject *
pycga1_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_CGA1(cga1_inv(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga1_lmat(((PyCga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga1_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga1_rmat(((PyCga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

/* The type of the centre span{1, w}, w the pseudoscalar with w*w = 1:
 * splitcomplex1, imported on first use so this extension loads on its own. */
static PyObject *cga1_center_type = NULL;

static PyObject *
pycga1_center(void)
{
    if (cga1_center_type == NULL) {
        cga1_center_type = algpy_import_attribute("algpy.algebras.multicomplex._splitcomplex1", "splitcomplex1");
        if (cga1_center_type == NULL) {
            return NULL;
        }
    }
    return Py_NewRef(cga1_center_type);
}

static PyObject *
pycga1_center_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex1 p[3];
    PyObject *center = pycga1_center(), *result;
    if (center == NULL) {
        return NULL;
    }
    cga1_center_minpoly(((PyCga1 *)self)->v, p);
    result = algpy_array_of_type(center, p, 3, sizeof(splitcomplex1));
    Py_DECREF(center);
    return result;
}

/* Metatype providing the class attribute cga1.center. */
static PyObject *
pycga1_meta_center(PyObject *Py_UNUSED(self), void *Py_UNUSED(closure))
{
    return pycga1_center();
}

static PyGetSetDef pycga1_meta_getset[] = {
    {"center", pycga1_meta_center, NULL, "Type of the centre span{1, w}: splitcomplex1.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PyCga1_MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.cga1.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = pycga1_meta_getset,
};

static PyObject *
pycga1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    cga1_minpoly(((PyCga1 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycga1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    cga1_eivals(((PyCga1 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycga1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga1);
    cga1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga1 *)PyArray_DATA(arr);
    cga1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga1_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_cga1);
    cga1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (cga1 *)PyArray_DATA(arr);
    cga1_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pycga1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_cga1);
    cga1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (cga1 *)PyArray_DATA(arr);
    cga1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycga1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    cga1 v = ((PyCga1 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pycga1_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_square(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_expm1(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_log1p(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_sin(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_cos(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_tan(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_asin(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_acos(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_atan(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_sinh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_cosh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_tanh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_asinh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_acosh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_atanh(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_re(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyCga1_FromValue(Py_TYPE(self), cga1_im(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_real(PyObject *self, void *closure)
{
    return pycga1_re(self, NULL);
}

static PyObject *
pycga1_imag(PyObject *self, void *closure)
{
    return pycga1_im(self, NULL);
}

static PyObject *pycga1_grades(PyObject *self, void *closure);

static PyGetSetDef pycga1_getset[] = {
    {"grade", pycga1_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pycga1_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pycga1_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pycga1_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(cga1_is_real(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(cga1_to_re(((PyCga1 *)self)->v));
}

static PyObject *
pycga1_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyCga1_FromValue((PyTypeObject *)cls, cga1_from_re(a));
}

static PyObject *
pycga1_lpow(PyObject *self, PyObject *other)
{
    cga1 b;
    int r;
    r = cga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga1_FromValue(Py_TYPE(self),
        cga1_lpow(((PyCga1 *)self)->v, b));
}

static PyObject *
pycga1_rpow(PyObject *self, PyObject *other)
{
    cga1 b;
    int r;
    r = cga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga1_FromValue(Py_TYPE(self),
        cga1_rpow(((PyCga1 *)self)->v, b));
}

static PyObject *
pycga1_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    cga1 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = cga1_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyCga1_FromValue(Py_TYPE(self),
        cga1_mix(((PyCga1 *)self)->v, b, t));
}

static PyObject *
pycga1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_cga1, 8);
}

static PyObject *
pycga1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_cga1, 8);
}

static PyObject *
pycga1_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_CGA1(cga1_ipow(((PyCga1 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pycga1_grade(PyObject *self, Py_ssize_t k)
{
    cga1 v = ((PyCga1 *)self)->v;
    switch (k) {
        case 0:
            return PyCga1_FromValue(Py_TYPE(self), cga1_grade0(v));
        case 1:
            return PyCga1_FromValue(Py_TYPE(self), cga1_grade1(v));
        case 2:
            return PyCga1_FromValue(Py_TYPE(self), cga1_grade2(v));
        case 3:
            return PyCga1_FromValue(Py_TYPE(self), cga1_grade3(v));
        default:
            return PyCga1_FromValue(Py_TYPE(self), cga1_zero());
    }
}

static PyObject *
pycga1_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pycga1_grade, 4);
}

static PyMethodDef pycga1_methods[] = {
    {"to_float", (PyCFunction)pycga1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pycga1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pycga1_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pycga1_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pycga1_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pycga1_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pycga1_square, METH_NOARGS, "square."},
    {"expm1", pycga1_expm1, METH_NOARGS, "expm1."},
    {"log1p", pycga1_log1p, METH_NOARGS, "log1p."},
    {"sin", pycga1_sin, METH_NOARGS, "sin."},
    {"cos", pycga1_cos, METH_NOARGS, "cos."},
    {"tan", pycga1_tan, METH_NOARGS, "tan."},
    {"asin", pycga1_asin, METH_NOARGS, "asin."},
    {"acos", pycga1_acos, METH_NOARGS, "acos."},
    {"atan", pycga1_atan, METH_NOARGS, "atan."},
    {"sinh", pycga1_sinh, METH_NOARGS, "sinh."},
    {"cosh", pycga1_cosh, METH_NOARGS, "cosh."},
    {"tanh", pycga1_tanh, METH_NOARGS, "tanh."},
    {"asinh", pycga1_asinh, METH_NOARGS, "asinh."},
    {"acosh", pycga1_acosh, METH_NOARGS, "acosh."},
    {"atanh", pycga1_atanh, METH_NOARGS, "atanh."},
    {"re", pycga1_re, METH_NOARGS, "re."},
    {"im", pycga1_im, METH_NOARGS, "im."},
    {"arcsin", pycga1_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pycga1_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pycga1_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pycga1_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pycga1_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pycga1_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pycga1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycga1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pycga1_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pycga1_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pycga1_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pycga1_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"dual", pycga1_dual, METH_NOARGS, "Poincare dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"basis", pycga1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pycga1_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pycga1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pycga1_norm, METH_NOARGS, "Norm as a Python complex."},
    {"lpow", pycga1_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pycga1_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pycga1_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pycga1_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycga1_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pycga1_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pycga1_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pycga1_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"center_minpoly", pycga1_center_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial over the centre, as an array of splitcomplex1 values in ascending order."},
    {"eivals", pycga1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycga1_exp, METH_NOARGS, "exp."},
    {"log", pycga1_log, METH_NOARGS, "log."},
    {"sqrt", pycga1_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycga1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycga1_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pycga1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyCga1 *)self)->v.data[i]);
}

static PyObject *
pycga1_subscript(PyObject *self, PyObject *key)
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
    return pycga1_item(self, i);
}

static int
pycga1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycga1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycga1_slots[] = {
    {Py_tp_dealloc, pycga1_dealloc},
    {Py_tp_traverse, pycga1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycga1_repr},
    {Py_tp_str, pycga1_repr},
    {Py_tp_hash, pycga1_hash},
    {Py_tp_doc, (void *)"Conformal geometric algebra in 1D with 8 float64 coefficients"},
    {Py_tp_richcompare, pycga1_richcompare},
    {Py_tp_methods, pycga1_methods},
    {Py_tp_getset, pycga1_getset},
    {Py_sq_length, pycga1_length},
    {Py_sq_item, pycga1_item},
    {Py_mp_subscript, pycga1_subscript},
    {Py_tp_new, pycga1_new},
    {Py_nb_add, pycga1_add},
    {Py_nb_subtract, pycga1_sub},
    {Py_nb_multiply, pycga1_mul},
    {Py_nb_true_divide, pycga1_div},
    {Py_nb_power, (void *)pycga1_pow},
    {Py_nb_positive, pycga1_pos},
    {Py_nb_negative, pycga1_neg},
    {Py_nb_or, pycga1_join},
    {Py_nb_and, pycga1_meet},
    {Py_nb_invert, pycga1_invert},
    {Py_nb_bool, pycga1_bool},
    {0, NULL}
};

static PyType_Spec pycga1_spec = {
    .name = "algpy.algebras.clifford.cga1",
    .basicsize = sizeof(PyCga1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycga1_slots,
};

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
npycga1_getitem(void *data, void *arr)
{
    cga1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_CGA1(v);
}

static int
npycga1_setitem(PyObject *item, void *data, void *arr)
{
    cga1 v;
    int r = cga1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to cga1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "cga1 needs exactly 8 components");
            return -1;
        }
        for (int i = 0; i < 8; i++) {
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
pycga1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(cga1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pycga1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycga1_copyswapn(dst, sizeof(cga1), src, sizeof(cga1), 1,
                         swap, arr);
}

static npy_bool
pycga1_nonzero(void *data, void *arr)
{
    cga1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycga1_arrfuncs;
typedef struct { char c; cga1 v; } pycga1_align_test;

static PyArray_DescrProto pycga1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(cga1), offsetof(pycga1_align_test, v), 0, 0, 0,
    &pycga1_arrfuncs,
};

#define CGA1_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pycga1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        cga1 *dst = (cga1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = cga1_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define CGA1_REAL_TYPES(X) \
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

CGA1_REAL_TYPES(CGA1_CAST_FROM)

static int
pycga1_register_casts(void)
{
#define REGISTER_CGA1_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_cga1,                    \
                                     pycga1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_cga1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    CGA1_REAL_TYPES(REGISTER_CGA1_CAST)
#undef REGISTER_CGA1_CAST
    return 0;
}

typedef cga1 (*cga1_unary_op)(cga1);
typedef cga1 (*cga1_binary_op)(cga1, cga1);
typedef int (*cga1_predicate)(cga1, cga1);

static void
loop_cga1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    cga1_unary_op op = (cga1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_cga1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    cga1_binary_op op = (cga1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga1 a, b, r;
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
loop_cga1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    cga1_predicate op = (cga1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        cga1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycga1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_cga1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_cga1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycga1_register_ufuncs(PyObject *numpy)
{
    return
        pycga1_register_loop(numpy, "add", loop_cga1_binary, 3, 0, cga1_add) ||
        pycga1_register_loop(numpy, "subtract", loop_cga1_binary, 3, 0, cga1_sub) ||
        pycga1_register_loop(numpy, "multiply", loop_cga1_binary, 3, 0, cga1_mul) ||
        pycga1_register_loop(numpy, "divide", loop_cga1_binary, 3, 0, cga1_div) ||
        pycga1_register_loop(numpy, "power", loop_cga1_binary, 3, 0, cga1_rpow) ||
        pycga1_register_loop(numpy, "positive", loop_cga1_unary, 2, 0, cga1_pos) ||
        pycga1_register_loop(numpy, "negative", loop_cga1_unary, 2, 0, cga1_neg) ||
        pycga1_register_loop(numpy, "conjugate", loop_cga1_unary, 2, 0, cga1_conjugate) ||
        pycga1_register_loop(numpy, "reciprocal", loop_cga1_unary, 2, 0, cga1_inv) ||
        pycga1_register_loop(numpy, "square", loop_cga1_unary, 2, 0, cga1_square) ||
        pycga1_register_loop(numpy, "exp", loop_cga1_unary, 2, 0, cga1_exp) ||
        pycga1_register_loop(numpy, "expm1", loop_cga1_unary, 2, 0, cga1_expm1) ||
        pycga1_register_loop(numpy, "log", loop_cga1_unary, 2, 0, cga1_log) ||
        pycga1_register_loop(numpy, "log1p", loop_cga1_unary, 2, 0, cga1_log1p) ||
        pycga1_register_loop(numpy, "sqrt", loop_cga1_unary, 2, 0, cga1_sqrt) ||
        pycga1_register_loop(numpy, "sin", loop_cga1_unary, 2, 0, cga1_sin) ||
        pycga1_register_loop(numpy, "cos", loop_cga1_unary, 2, 0, cga1_cos) ||
        pycga1_register_loop(numpy, "tan", loop_cga1_unary, 2, 0, cga1_tan) ||
        pycga1_register_loop(numpy, "arcsin", loop_cga1_unary, 2, 0, cga1_asin) ||
        pycga1_register_loop(numpy, "arccos", loop_cga1_unary, 2, 0, cga1_acos) ||
        pycga1_register_loop(numpy, "arctan", loop_cga1_unary, 2, 0, cga1_atan) ||
        pycga1_register_loop(numpy, "sinh", loop_cga1_unary, 2, 0, cga1_sinh) ||
        pycga1_register_loop(numpy, "cosh", loop_cga1_unary, 2, 0, cga1_cosh) ||
        pycga1_register_loop(numpy, "tanh", loop_cga1_unary, 2, 0, cga1_tanh) ||
        pycga1_register_loop(numpy, "arcsinh", loop_cga1_unary, 2, 0, cga1_asinh) ||
        pycga1_register_loop(numpy, "arccosh", loop_cga1_unary, 2, 0, cga1_acosh) ||
        pycga1_register_loop(numpy, "arctanh", loop_cga1_unary, 2, 0, cga1_atanh) ||
        pycga1_register_loop(numpy, "equal", loop_cga1_predicate, 3, 1, cga1_equal) ||
        pycga1_register_loop(numpy, "bitwise_or", loop_cga1_binary, 3, 0, cga1_join) ||
        pycga1_register_loop(numpy, "bitwise_and", loop_cga1_binary, 3, 0, cga1_meet) ||
        pycga1_register_loop(numpy, "invert", loop_cga1_unary, 2, 0, cga1_dual) ||
        pycga1_register_loop(numpy, "not_equal", loop_cga1_predicate, 3, 1, cga1_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pycga1_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pycga1_register_loop(ufuncs, "ldual", loop_cga1_unary, 2, 0, cga1_ldual) ||
        pycga1_register_loop(ufuncs, "rdual", loop_cga1_unary, 2, 0, cga1_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_cga1", "cga1 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__cga1(void)
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
    PyCga1_MetaType.tp_base = &PyType_Type;  /* not a constant initializer on MSVC */
    if (PyType_Ready(&PyCga1_MetaType) < 0) {
        goto fail;
    }
    PyCga1_Type = PyType_FromModuleAndSpec(
        m, &pycga1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyCga1_Type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(PyCga1_Type, &PyCga1_MetaType);
    PyArray_InitArrFuncs(&pycga1_arrfuncs);
    pycga1_arrfuncs.getitem = npycga1_getitem;
    pycga1_arrfuncs.setitem = npycga1_setitem;
    pycga1_arrfuncs.copyswapn = pycga1_copyswapn;
    pycga1_arrfuncs.copyswap = pycga1_copyswap;
    pycga1_arrfuncs.nonzero = pycga1_nonzero;
    pycga1_descr_proto.typeobj = (PyTypeObject *)PyCga1_Type;
    Py_SET_TYPE((PyObject *)&pycga1_descr_proto, &PyArrayDescr_Type);
    npy_cga1 = PyArray_RegisterDataType(&pycga1_descr_proto);
    if (npy_cga1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_cga1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyCga1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyCga1_Type, 8) < 0) {
        goto fail;
    }
    if (pycga1_register_casts() < 0 ||
            pycga1_register_ufuncs(numpy) < 0 ||
            pycga1_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "cga1", PyCga1_Type) < 0) {
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
