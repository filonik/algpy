/*
 * pga2 numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in pga2.h.
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

#include "pga2.h"

typedef struct {
    PyObject_HEAD
    pga2 v;
} PyPga2;

static PyObject *PyPga2_Type = NULL;
static int npy_pga2 = -1;

static PyObject *
PyPga2_FromValue(PyTypeObject *type, pga2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyPga2 *p = (PyPga2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_PGA2(v) \
    PyPga2_FromValue((PyTypeObject *)PyPga2_Type, (v))

static int
pga2_from_object(PyObject *o, pga2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyPga2_Type)) {
        *out = ((PyPga2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = pga2_from_re(d);
        return 1;
    }
    return 0;
}

static PyObject *
pypga2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyPga2_FromValue(type, pga2_make(d));
}

static PyObject *
pypga2_repr(PyObject *self)
{
    pga2 v = ((PyPga2 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("pga2(%s, %s, %s, %s, %s, %s, %s, %s)",
        s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pypga2_hash(PyObject *self)
{
    pga2 v = ((PyPga2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_PGA2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (pga2_is_real(v)) {
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
pypga2_richcompare(PyObject *a, PyObject *b, int op)
{
    pga2 x, y;
    int r;
    if ((r = pga2_from_object(a, &x)) <= 0 ||
            (r = pga2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(pga2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!pga2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define PGA2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pypga2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        pga2 x, y;                                                       \
        int r;                                                               \
        if ((r = pga2_from_object(a, &x)) <= 0 ||                        \
                (r = pga2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_PGA2(expr(x, y));                                     \
    }

PGA2_BINARY(add, pga2_add)
PGA2_BINARY(sub, pga2_sub)
PGA2_BINARY(mul, pga2_mul)
PGA2_BINARY(div, pga2_div)
PGA2_BINARY(join, pga2_join)
PGA2_BINARY(meet, pga2_meet)

static PyObject *
pypga2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    pga2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "pga2 does not support modular power");
        return NULL;
    }
    if ((r = pga2_from_object(a, &x)) <= 0 ||
            (r = pga2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_PGA2(pga2_rpow(x, y));
}

#define PGA2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pypga2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_PGA2(expr(((PyPga2 *)self)->v));                  \
    }

PGA2_UNARY(pos, pga2_pos)
PGA2_UNARY(neg, pga2_neg)
PGA2_UNARY(invert, pga2_dual)

static PyObject *
pypga2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue((PyTypeObject *)cls, pga2_zero());
}

static PyObject *
pypga2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue((PyTypeObject *)cls, pga2_one());
}

static PyObject *
pypga2_bottom(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue((PyTypeObject *)cls, pga2_bottom());
}

static PyObject *
pypga2_top(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue((PyTypeObject *)cls, pga2_top());
}

static PyObject *
pypga2_ldual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_ldual(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_rdual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_rdual(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_dual(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_dual(((PyPga2 *)self)->v));
}

static int
pypga2_bool(PyObject *self)
{
    pga2 v = ((PyPga2 *)self)->v;
    return vec8_any_nonzero(v);
}

#define PGA2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pypga2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_PGA2(op(((PyPga2 *)self)->v));                   \
    }

PGA2_METHOD(conjugate, pga2_conjugate)
PGA2_METHOD(exp, pga2_exp)
PGA2_METHOD(log, pga2_log)
PGA2_METHOD(sqrt, pga2_sqrt)

static PyObject *
pypga2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga2_norm(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_PGA2(pga2_inv(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_lmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga2_lmat(((PyPga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga2_rmat(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyObject *arr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga2_rmat(((PyPga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

/* The type of the centre span{1, w}, w the pseudoscalar with w*w = 0:
 * dual1, imported on first use so this extension loads on its own. */
static PyObject *pga2_center_type = NULL;

static PyObject *
pypga2_center(void)
{
    if (pga2_center_type == NULL) {
        pga2_center_type = algpy_import_attribute("algpy.algebras.multicomplex._dual1", "dual1");
        if (pga2_center_type == NULL) {
            return NULL;
        }
    }
    return Py_NewRef(pga2_center_type);
}

static PyObject *
pypga2_center_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    dual1 p[3];
    PyObject *center = pypga2_center(), *result;
    if (center == NULL) {
        return NULL;
    }
    pga2_center_minpoly(((PyPga2 *)self)->v, p);
    result = algpy_array_of_type(center, p, 3, sizeof(dual1));
    Py_DECREF(center);
    return result;
}

/* Metatype providing the class attribute pga2.center. */
static PyObject *
pypga2_meta_center(PyObject *Py_UNUSED(self), void *Py_UNUSED(closure))
{
    return pypga2_center();
}

static PyGetSetDef pypga2_meta_getset[] = {
    {"center", pypga2_meta_center, NULL, "Type of the centre span{1, w}: dual1.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyTypeObject PyPga2_MetaType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "algpy.pga2.type",
    .tp_basicsize = sizeof(PyHeapTypeObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_getset = pypga2_meta_getset,
};

static PyObject *
pypga2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    pga2_minpoly(((PyPga2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pypga2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    pga2_eivals(((PyPga2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pypga2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga2);
    pga2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga2 *)PyArray_DATA(arr);
    pga2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga2_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_pga2);
    pga2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (pga2 *)PyArray_DATA(arr);
    pga2_gens(basis);
    return (PyObject *)arr;
}

static PyObject *
pypga2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {8, 8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_pga2);
    pga2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (pga2 *)PyArray_DATA(arr);
    pga2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pypga2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    pga2 v = ((PyPga2 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
        v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pypga2_square(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_square(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_expm1(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_expm1(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_log1p(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_log1p(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_sin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_sin(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_cos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_cos(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_tan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_tan(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_asin(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_asin(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_acos(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_acos(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_atan(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_atan(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_sinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_sinh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_cosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_cosh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_tanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_tanh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_asinh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_asinh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_acosh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_acosh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_atanh(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_atanh(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_re(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_im(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyPga2_FromValue(Py_TYPE(self), pga2_im(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_real(PyObject *self, void *closure)
{
    return pypga2_re(self, NULL);
}

static PyObject *
pypga2_imag(PyObject *self, void *closure)
{
    return pypga2_im(self, NULL);
}

static PyObject *pypga2_grades(PyObject *self, void *closure);

static PyGetSetDef pypga2_getset[] = {
    {"grade", pypga2_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {"real", pypga2_real, NULL, "Real part in the algebra.", NULL},
    {"imag", pypga2_imag, NULL, "Non-real part in the algebra.", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pypga2_is_real(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyBool_FromLong(pga2_is_real(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_to_re(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(pga2_to_re(((PyPga2 *)self)->v));
}

static PyObject *
pypga2_from_re(PyObject *cls, PyObject *value)
{
    double a = PyFloat_AsDouble(value);
    if (a == -1.0 && PyErr_Occurred()) {
        return NULL;
    }
    return PyPga2_FromValue((PyTypeObject *)cls, pga2_from_re(a));
}

static PyObject *
pypga2_lpow(PyObject *self, PyObject *other)
{
    pga2 b;
    int r;
    r = pga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "lpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga2_FromValue(Py_TYPE(self),
        pga2_lpow(((PyPga2 *)self)->v, b));
}

static PyObject *
pypga2_rpow(PyObject *self, PyObject *other)
{
    pga2 b;
    int r;
    r = pga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "rpow expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga2_FromValue(Py_TYPE(self),
        pga2_rpow(((PyPga2 *)self)->v, b));
}

static PyObject *
pypga2_mix(PyObject *self, PyObject *args)
{
    PyObject *other;
    pga2 b;
    double t;
    int r;
    if (!PyArg_ParseTuple(args, "Od", &other, &t)) {
        return NULL;
    }
    r = pga2_from_object(other, &b);
    if (r < 0) {
        return NULL;
    }
    if (r == 0) {
        PyErr_SetString(PyExc_TypeError, "mix expects a value of the same algebra or a real number");
        return NULL;
    }
    return PyPga2_FromValue(Py_TYPE(self),
        pga2_mix(((PyPga2 *)self)->v, b, t));
}

static PyObject *
pypga2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_pga2, 8);
}

static PyObject *
pypga2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_pga2, 8);
}

static PyObject *
pypga2_ipow(PyObject *self, PyObject *arg)
{
    long long exponent = PyLong_AsLongLong(arg);
    if (exponent == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return NEW_PGA2(pga2_ipow(((PyPga2 *)self)->v, (int64_t)exponent));
}

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pypga2_grade(PyObject *self, Py_ssize_t k)
{
    pga2 v = ((PyPga2 *)self)->v;
    switch (k) {
        case 0:
            return PyPga2_FromValue(Py_TYPE(self), pga2_grade0(v));
        case 1:
            return PyPga2_FromValue(Py_TYPE(self), pga2_grade1(v));
        case 2:
            return PyPga2_FromValue(Py_TYPE(self), pga2_grade2(v));
        case 3:
            return PyPga2_FromValue(Py_TYPE(self), pga2_grade3(v));
        default:
            return PyPga2_FromValue(Py_TYPE(self), pga2_zero());
    }
}

static PyObject *
pypga2_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pypga2_grade, 4);
}

static PyMethodDef pypga2_methods[] = {
    {"to_float", (PyCFunction)pypga2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert algebra values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pypga2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the algebra dtype."},
    {"is_real", pypga2_is_real, METH_NOARGS, "Whether all non-real components are zero."},
    {"to_re", pypga2_to_re, METH_NOARGS, "Real coefficient."},
    {"from_re", pypga2_from_re, METH_O | METH_CLASS, "Embed a real number."},
    {"mix", pypga2_mix, METH_VARARGS, "mix(other, t): self * (inv(self) * other)**t."},
    {"square", pypga2_square, METH_NOARGS, "square."},
    {"expm1", pypga2_expm1, METH_NOARGS, "expm1."},
    {"log1p", pypga2_log1p, METH_NOARGS, "log1p."},
    {"sin", pypga2_sin, METH_NOARGS, "sin."},
    {"cos", pypga2_cos, METH_NOARGS, "cos."},
    {"tan", pypga2_tan, METH_NOARGS, "tan."},
    {"asin", pypga2_asin, METH_NOARGS, "asin."},
    {"acos", pypga2_acos, METH_NOARGS, "acos."},
    {"atan", pypga2_atan, METH_NOARGS, "atan."},
    {"sinh", pypga2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pypga2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pypga2_tanh, METH_NOARGS, "tanh."},
    {"asinh", pypga2_asinh, METH_NOARGS, "asinh."},
    {"acosh", pypga2_acosh, METH_NOARGS, "acosh."},
    {"atanh", pypga2_atanh, METH_NOARGS, "atanh."},
    {"re", pypga2_re, METH_NOARGS, "re."},
    {"im", pypga2_im, METH_NOARGS, "im."},
    {"arcsin", pypga2_asin, METH_NOARGS, "Alias for asin."},
    {"arccos", pypga2_acos, METH_NOARGS, "Alias for acos."},
    {"arctan", pypga2_atan, METH_NOARGS, "Alias for atan."},
    {"arcsinh", pypga2_asinh, METH_NOARGS, "Alias for asinh."},
    {"arccosh", pypga2_acosh, METH_NOARGS, "Alias for acosh."},
    {"arctanh", pypga2_atanh, METH_NOARGS, "Alias for atanh."},

    {"zero", pypga2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pypga2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"bottom", pypga2_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pypga2_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"ldual", pypga2_ldual, METH_NOARGS, "Left Poincare dual: a | a.ldual() == top for basis blades."},
    {"rdual", pypga2_rdual, METH_NOARGS, "Right Poincare dual: a.rdual() | a == top for basis blades."},
    {"dual", pypga2_dual, METH_NOARGS, "Poincare dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"basis", pypga2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"gens", pypga2_gens, METH_NOARGS | METH_CLASS, "Algebra generators in basis order; excludes the identity."},
    {"coeffs", pypga2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pypga2_norm, METH_NOARGS, "Norm as a Python float."},
    {"lpow", pypga2_lpow, METH_O, "lpow(exponent): exp(exponent * log(self))."},
    {"rpow", pypga2_rpow, METH_O, "rpow(exponent): exp(log(self) * exponent); default power."},
    {"ipow", pypga2_ipow, METH_O, "ipow(exponent): integer power by squaring; exponent must fit signed 64 bits."},
    {"inv", pypga2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pypga2_conjugate, METH_NOARGS, "Clifford conjugation."},
    {"lmat", pypga2_lmat, METH_NOARGS, "Left multiplication matrix acting on coefficient columns."},
    {"rmat", pypga2_rmat, METH_NOARGS, "Right multiplication matrix acting on coefficient rows."},
    {"minpoly", pypga2_minpoly, METH_NOARGS, "Generic minimal polynomial coefficients."},
    {"center_minpoly", pypga2_center_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial over the centre, as an array of dual1 values in ascending order."},
    {"eivals", pypga2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pypga2_exp, METH_NOARGS, "exp."},
    {"log", pypga2_log, METH_NOARGS, "log."},
    {"sqrt", pypga2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pypga2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pypga2_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pypga2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyPga2 *)self)->v.data[i]);
}

static PyObject *
pypga2_subscript(PyObject *self, PyObject *key)
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
    return pypga2_item(self, i);
}

static int
pypga2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pypga2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pypga2_slots[] = {
    {Py_tp_dealloc, pypga2_dealloc},
    {Py_tp_traverse, pypga2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pypga2_repr},
    {Py_tp_str, pypga2_repr},
    {Py_tp_hash, pypga2_hash},
    {Py_tp_doc, (void *)"Projective geometric algebra in 2D with 8 float64 coefficients"},
    {Py_tp_richcompare, pypga2_richcompare},
    {Py_tp_methods, pypga2_methods},
    {Py_tp_getset, pypga2_getset},
    {Py_sq_length, pypga2_length},
    {Py_sq_item, pypga2_item},
    {Py_mp_subscript, pypga2_subscript},
    {Py_tp_new, pypga2_new},
    {Py_nb_add, pypga2_add},
    {Py_nb_subtract, pypga2_sub},
    {Py_nb_multiply, pypga2_mul},
    {Py_nb_true_divide, pypga2_div},
    {Py_nb_power, (void *)pypga2_pow},
    {Py_nb_positive, pypga2_pos},
    {Py_nb_negative, pypga2_neg},
    {Py_nb_or, pypga2_join},
    {Py_nb_and, pypga2_meet},
    {Py_nb_invert, pypga2_invert},
    {Py_nb_bool, pypga2_bool},
    {0, NULL}
};

static PyType_Spec pypga2_spec = {
    .name = "algpy.algebras.clifford.pga2",
    .basicsize = sizeof(PyPga2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pypga2_slots,
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
npypga2_getitem(void *data, void *arr)
{
    pga2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_PGA2(v);
}

static int
npypga2_setitem(PyObject *item, void *data, void *arr)
{
    pga2 v;
    int r = pga2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to pga2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "pga2 needs exactly 8 components");
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
pypga2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(pga2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pypga2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pypga2_copyswapn(dst, sizeof(pga2), src, sizeof(pga2), 1,
                         swap, arr);
}

static npy_bool
pypga2_nonzero(void *data, void *arr)
{
    pga2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pypga2_arrfuncs;
typedef struct { char c; pga2 v; } pypga2_align_test;

static PyArray_DescrProto pypga2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(pga2), offsetof(pypga2_align_test, v), 0, 0, 0,
    &pypga2_arrfuncs,
};

#define PGA2_CAST_FROM(name, type, convert)                                       \
    static void                                                              \
    pypga2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        pga2 *dst = (pga2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = pga2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define PGA2_REAL_TYPES(X) \
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

PGA2_REAL_TYPES(PGA2_CAST_FROM)

static int
pypga2_register_casts(void)
{
#define REGISTER_PGA2_CAST(name, type, convert)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_pga2,                    \
                                     pypga2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_pga2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    PGA2_REAL_TYPES(REGISTER_PGA2_CAST)
#undef REGISTER_PGA2_CAST
    return 0;
}

typedef pga2 (*pga2_unary_op)(pga2);
typedef pga2 (*pga2_binary_op)(pga2, pga2);
typedef int (*pga2_predicate)(pga2, pga2);

static void
loop_pga2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    pga2_unary_op op = (pga2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_pga2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    pga2_binary_op op = (pga2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga2 a, b, r;
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
loop_pga2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    pga2_predicate op = (pga2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        pga2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pypga2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_pga2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_pga2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pypga2_register_ufuncs(PyObject *numpy)
{
    return
        pypga2_register_loop(numpy, "add", loop_pga2_binary, 3, 0, pga2_add) ||
        pypga2_register_loop(numpy, "subtract", loop_pga2_binary, 3, 0, pga2_sub) ||
        pypga2_register_loop(numpy, "multiply", loop_pga2_binary, 3, 0, pga2_mul) ||
        pypga2_register_loop(numpy, "divide", loop_pga2_binary, 3, 0, pga2_div) ||
        pypga2_register_loop(numpy, "power", loop_pga2_binary, 3, 0, pga2_rpow) ||
        pypga2_register_loop(numpy, "positive", loop_pga2_unary, 2, 0, pga2_pos) ||
        pypga2_register_loop(numpy, "negative", loop_pga2_unary, 2, 0, pga2_neg) ||
        pypga2_register_loop(numpy, "conjugate", loop_pga2_unary, 2, 0, pga2_conjugate) ||
        pypga2_register_loop(numpy, "reciprocal", loop_pga2_unary, 2, 0, pga2_inv) ||
        pypga2_register_loop(numpy, "square", loop_pga2_unary, 2, 0, pga2_square) ||
        pypga2_register_loop(numpy, "exp", loop_pga2_unary, 2, 0, pga2_exp) ||
        pypga2_register_loop(numpy, "expm1", loop_pga2_unary, 2, 0, pga2_expm1) ||
        pypga2_register_loop(numpy, "log", loop_pga2_unary, 2, 0, pga2_log) ||
        pypga2_register_loop(numpy, "log1p", loop_pga2_unary, 2, 0, pga2_log1p) ||
        pypga2_register_loop(numpy, "sqrt", loop_pga2_unary, 2, 0, pga2_sqrt) ||
        pypga2_register_loop(numpy, "sin", loop_pga2_unary, 2, 0, pga2_sin) ||
        pypga2_register_loop(numpy, "cos", loop_pga2_unary, 2, 0, pga2_cos) ||
        pypga2_register_loop(numpy, "tan", loop_pga2_unary, 2, 0, pga2_tan) ||
        pypga2_register_loop(numpy, "arcsin", loop_pga2_unary, 2, 0, pga2_asin) ||
        pypga2_register_loop(numpy, "arccos", loop_pga2_unary, 2, 0, pga2_acos) ||
        pypga2_register_loop(numpy, "arctan", loop_pga2_unary, 2, 0, pga2_atan) ||
        pypga2_register_loop(numpy, "sinh", loop_pga2_unary, 2, 0, pga2_sinh) ||
        pypga2_register_loop(numpy, "cosh", loop_pga2_unary, 2, 0, pga2_cosh) ||
        pypga2_register_loop(numpy, "tanh", loop_pga2_unary, 2, 0, pga2_tanh) ||
        pypga2_register_loop(numpy, "arcsinh", loop_pga2_unary, 2, 0, pga2_asinh) ||
        pypga2_register_loop(numpy, "arccosh", loop_pga2_unary, 2, 0, pga2_acosh) ||
        pypga2_register_loop(numpy, "arctanh", loop_pga2_unary, 2, 0, pga2_atanh) ||
        pypga2_register_loop(numpy, "equal", loop_pga2_predicate, 3, 1, pga2_equal) ||
        pypga2_register_loop(numpy, "bitwise_or", loop_pga2_binary, 3, 0, pga2_join) ||
        pypga2_register_loop(numpy, "bitwise_and", loop_pga2_binary, 3, 0, pga2_meet) ||
        pypga2_register_loop(numpy, "invert", loop_pga2_unary, 2, 0, pga2_dual) ||
        pypga2_register_loop(numpy, "not_equal", loop_pga2_predicate, 3, 1, pga2_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pypga2_register_lattice_ufuncs(PyObject *ufuncs)
{
    return
        pypga2_register_loop(ufuncs, "ldual", loop_pga2_unary, 2, 0, pga2_ldual) ||
        pypga2_register_loop(ufuncs, "rdual", loop_pga2_unary, 2, 0, pga2_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_pga2", "pga2 number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__pga2(void)
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
    PyPga2_MetaType.tp_base = &PyType_Type;  /* not a constant initializer on MSVC */
    if (PyType_Ready(&PyPga2_MetaType) < 0) {
        goto fail;
    }
    PyPga2_Type = PyType_FromModuleAndSpec(
        m, &pypga2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyPga2_Type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(PyPga2_Type, &PyPga2_MetaType);
    PyArray_InitArrFuncs(&pypga2_arrfuncs);
    pypga2_arrfuncs.getitem = npypga2_getitem;
    pypga2_arrfuncs.setitem = npypga2_setitem;
    pypga2_arrfuncs.copyswapn = pypga2_copyswapn;
    pypga2_arrfuncs.copyswap = pypga2_copyswap;
    pypga2_arrfuncs.nonzero = pypga2_nonzero;
    pypga2_descr_proto.typeobj = (PyTypeObject *)PyPga2_Type;
    Py_SET_TYPE((PyObject *)&pypga2_descr_proto, &PyArrayDescr_Type);
    npy_pga2 = PyArray_RegisterDataType(&pypga2_descr_proto);
    if (npy_pga2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_pga2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyPga2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyPga2_Type, 8) < 0) {
        goto fail;
    }
    if (pypga2_register_casts() < 0 ||
            pypga2_register_ufuncs(numpy) < 0 ||
            pypga2_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "pga2", PyPga2_Type) < 0) {
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
