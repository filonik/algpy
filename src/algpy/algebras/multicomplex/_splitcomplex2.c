/*
 * Bi-split-complex numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in splitcomplex2.h.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "splitcomplex2.h"
#include "../norm.h"

typedef struct {
    PyObject_HEAD
    splitcomplex2 v;
} PySplitcomplex2;

static PyObject *PySplitcomplex2_Type = NULL;
static int npy_splitcomplex2 = -1;

static PyObject *
PySplitcomplex2_FromValue(PyTypeObject *type, splitcomplex2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySplitcomplex2 *p = (PySplitcomplex2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SPLITCOMPLEX2(v) \
    PySplitcomplex2_FromValue((PyTypeObject *)PySplitcomplex2_Type, (v))

static int
splitcomplex2_from_object(PyObject *o, splitcomplex2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySplitcomplex2_Type)) {
        *out = ((PySplitcomplex2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = splitcomplex2_make(d, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pysplitcomplex2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"w", "x", "y", "z", NULL};
    double w = 0, x = 0, y = 0, z = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                     &w, &x, &y, &z)) {
        return NULL;
    }
    return PySplitcomplex2_FromValue(type, splitcomplex2_make(w, x, y, z));
}

static PyObject *
pysplitcomplex2_repr(PyObject *self)
{
    splitcomplex2 v = ((PySplitcomplex2 *)self)->v;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("splitcomplex2(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pysplitcomplex2_hash(PyObject *self)
{
    splitcomplex2 v = ((PySplitcomplex2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (v.data[1] == 0 && v.data[2] == 0 && v.data[3] == 0) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddd)", v.data[0], v.data[1], v.data[2], v.data[3]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pysplitcomplex2_richcompare(PyObject *a, PyObject *b, int op)
{
    splitcomplex2 x, y;
    int r;
    if ((r = splitcomplex2_from_object(a, &x)) <= 0 ||
            (r = splitcomplex2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(vec4_all_eq(x, y));
        case Py_NE:
            return PyBool_FromLong(!vec4_all_eq(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SPLITCOMPLEX2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysplitcomplex2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        splitcomplex2 x, y;                                                       \
        int r;                                                               \
        if ((r = splitcomplex2_from_object(a, &x)) <= 0 ||                        \
                (r = splitcomplex2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SPLITCOMPLEX2(expr(x, y));                                     \
    }

SPLITCOMPLEX2_BINARY(add, vec4_add)
SPLITCOMPLEX2_BINARY(sub, vec4_sub)
SPLITCOMPLEX2_BINARY(mul, splitcomplex2_mul)
SPLITCOMPLEX2_BINARY(div, splitcomplex2_div)

static PyObject *
pysplitcomplex2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    splitcomplex2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "splitcomplex2 does not support modular power");
        return NULL;
    }
    if ((r = splitcomplex2_from_object(a, &x)) <= 0 ||
            (r = splitcomplex2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SPLITCOMPLEX2(splitcomplex2_pow(x, y));
}

#define SPLITCOMPLEX2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysplitcomplex2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SPLITCOMPLEX2(expr(((PySplitcomplex2 *)self)->v));                  \
    }

SPLITCOMPLEX2_UNARY(pos, vec4_pos)
SPLITCOMPLEX2_UNARY(neg, vec4_neg)

static PyObject *
pysplitcomplex2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex2_FromValue((PyTypeObject *)cls, vec4_zero());
}

static PyObject *
pysplitcomplex2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex2_FromValue((PyTypeObject *)cls, splitcomplex2_one());
}

static int
pysplitcomplex2_bool(PyObject *self)
{
    splitcomplex2 v = ((PySplitcomplex2 *)self)->v;
    return vec4_any_nonzero(v);
}

#define SPLITCOMPLEX2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysplitcomplex2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SPLITCOMPLEX2(op(((PySplitcomplex2 *)self)->v));                   \
    }

SPLITCOMPLEX2_METHOD(conjugate, splitcomplex2_conjugate)
SPLITCOMPLEX2_METHOD(exp, splitcomplex2_exp)
SPLITCOMPLEX2_METHOD(log, splitcomplex2_log)
SPLITCOMPLEX2_METHOD(sqrt, splitcomplex2_sqrt)
SPLITCOMPLEX2_METHOD(sin, splitcomplex2_sin)
SPLITCOMPLEX2_METHOD(cos, splitcomplex2_cos)
SPLITCOMPLEX2_METHOD(tan, splitcomplex2_tan)
SPLITCOMPLEX2_METHOD(sinh, splitcomplex2_sinh)
SPLITCOMPLEX2_METHOD(cosh, splitcomplex2_cosh)
SPLITCOMPLEX2_METHOD(tanh, splitcomplex2_tanh)

static PyObject *
pysplitcomplex2_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex1 n = splitcomplex2_norm2(((PySplitcomplex2 *)self)->v);
    return Py_BuildValue("(dd)", n.data[0], n.data[1]);
}

static PyObject *
pysplitcomplex2_norm4(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(splitcomplex2_norm4(((PySplitcomplex2 *)self)->v));
}

static PyObject *
pysplitcomplex2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return algpy_norm_to_py(splitcomplex2_norm(((PySplitcomplex2 *)self)->v));
}

static PyObject *
pysplitcomplex2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_SPLITCOMPLEX2(splitcomplex2_inv(((PySplitcomplex2 *)self)->v));
}

static PyObject *
pysplitcomplex2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex2_minpoly(((PySplitcomplex2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pysplitcomplex2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    splitcomplex2_eivals(((PySplitcomplex2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pysplitcomplex2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex2);
    splitcomplex2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex2 *)PyArray_DATA(arr);
    splitcomplex2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_splitcomplex2);
    splitcomplex2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (splitcomplex2 *)PyArray_DATA(arr);
    splitcomplex2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex2 v = ((PySplitcomplex2 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyMethodDef pysplitcomplex2_methods[] = {
    {"zero", pysplitcomplex2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysplitcomplex2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysplitcomplex2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pysplitcomplex2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pysplitcomplex2_norm2, METH_NOARGS, "Squared norm (a dual pair)."},
    {"norm4", pysplitcomplex2_norm4, METH_NOARGS, "Fourth power of the norm."},
    {"norm", pysplitcomplex2_norm, METH_NOARGS, "Norm."},
    {"inv", pysplitcomplex2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pysplitcomplex2_conjugate, METH_NOARGS, "Conjugate of the e1 part."},
    {"minpoly", pysplitcomplex2_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pysplitcomplex2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysplitcomplex2_exp, METH_NOARGS, "exp."},
    {"log", pysplitcomplex2_log, METH_NOARGS, "log."},
    {"sqrt", pysplitcomplex2_sqrt, METH_NOARGS, "sqrt."},
    {"sin", pysplitcomplex2_sin, METH_NOARGS, "sin."},
    {"cos", pysplitcomplex2_cos, METH_NOARGS, "cos."},
    {"tan", pysplitcomplex2_tan, METH_NOARGS, "tan."},
    {"sinh", pysplitcomplex2_sinh, METH_NOARGS, "sinh."},
    {"cosh", pysplitcomplex2_cosh, METH_NOARGS, "cosh."},
    {"tanh", pysplitcomplex2_tanh, METH_NOARGS, "tanh."},
    {"__reduce__", pysplitcomplex2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysplitcomplex2_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pysplitcomplex2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySplitcomplex2 *)self)->v.data[i]);
}

static PyObject *
pysplitcomplex2_subscript(PyObject *self, PyObject *key)
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
        i += 4;
    }
    return pysplitcomplex2_item(self, i);
}

static int
pysplitcomplex2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysplitcomplex2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysplitcomplex2_slots[] = {
    {Py_tp_dealloc, pysplitcomplex2_dealloc},
    {Py_tp_traverse, pysplitcomplex2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysplitcomplex2_repr},
    {Py_tp_str, pysplitcomplex2_repr},
    {Py_tp_hash, pysplitcomplex2_hash},
    {Py_tp_doc, (void *)"Bi-split-complex number with float64 components"},
    {Py_tp_richcompare, pysplitcomplex2_richcompare},
    {Py_tp_methods, pysplitcomplex2_methods},
    {Py_sq_length, pysplitcomplex2_length},
    {Py_sq_item, pysplitcomplex2_item},
    {Py_mp_subscript, pysplitcomplex2_subscript},
    {Py_tp_new, pysplitcomplex2_new},
    {Py_nb_add, pysplitcomplex2_add},
    {Py_nb_subtract, pysplitcomplex2_sub},
    {Py_nb_multiply, pysplitcomplex2_mul},
    {Py_nb_true_divide, pysplitcomplex2_div},
    {Py_nb_power, (void *)pysplitcomplex2_pow},
    {Py_nb_positive, pysplitcomplex2_pos},
    {Py_nb_negative, pysplitcomplex2_neg},
    {Py_nb_bool, pysplitcomplex2_bool},
    {0, NULL}
};

static PyType_Spec pysplitcomplex2_spec = {
    .name = "algpy.algebras.multicomplex.splitcomplex2",
    .basicsize = sizeof(PySplitcomplex2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysplitcomplex2_slots,
};

static PyObject *
npysplitcomplex2_getitem(void *data, void *arr)
{
    splitcomplex2 v;
    memcpy(&v, data, sizeof(v));
    return NEW_SPLITCOMPLEX2(v);
}

static int
npysplitcomplex2_setitem(PyObject *item, void *data, void *arr)
{
    splitcomplex2 v;
    int r = splitcomplex2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to splitcomplex2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "splitcomplex2 needs exactly 4 components");
            return -1;
        }
        for (int i = 0; i < 4; i++) {
            v.data[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        }
        Py_DECREF(seq);
        if (PyErr_Occurred()) {
            return -1;
        }
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pysplitcomplex2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(splitcomplex2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 4; j++) {
                char *d = p + 8 * j;
                int b;
                for (b = 0; b < 4; b++) {
                    char t = d[b];
                    d[b] = d[7 - b];
                    d[7 - b] = t;
                }
            }
        }
    }
}

static void
pysplitcomplex2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysplitcomplex2_copyswapn(dst, sizeof(splitcomplex2), src, sizeof(splitcomplex2), 1,
                         swap, arr);
}

static npy_bool
pysplitcomplex2_nonzero(void *data, void *arr)
{
    splitcomplex2 v;
    memcpy(&v, data, sizeof(v));
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysplitcomplex2_arrfuncs;
typedef struct { char c; splitcomplex2 v; } pysplitcomplex2_align_test;

static PyArray_DescrProto pysplitcomplex2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(splitcomplex2), offsetof(pysplitcomplex2_align_test, v), 0, 0, 0,
    &pysplitcomplex2_arrfuncs,
};

#define SPLITCOMPLEX2_CAST_FROM(name, type)                                       \
    static void                                                              \
    pysplitcomplex2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        splitcomplex2 *dst = (splitcomplex2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = splitcomplex2_make((double)src[i], 0, 0, 0);                      \
        }                                                                    \
    }

#define SPLITCOMPLEX2_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

SPLITCOMPLEX2_REAL_TYPES(SPLITCOMPLEX2_CAST_FROM)

static int
pysplitcomplex2_register_casts(void)
{
#define REGISTER_SPLITCOMPLEX2_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_splitcomplex2,                    \
                                     pysplitcomplex2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_splitcomplex2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SPLITCOMPLEX2_REAL_TYPES(REGISTER_SPLITCOMPLEX2_CAST)
#undef REGISTER_SPLITCOMPLEX2_CAST
    return 0;
}

typedef splitcomplex2 (*splitcomplex2_unary_op)(splitcomplex2);
typedef splitcomplex2 (*splitcomplex2_binary_op)(splitcomplex2, splitcomplex2);
typedef int (*splitcomplex2_predicate)(splitcomplex2, splitcomplex2);

static void
loop_splitcomplex2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    splitcomplex2_unary_op op = (splitcomplex2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_splitcomplex2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    splitcomplex2_binary_op op = (splitcomplex2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex2 a, b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = op(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
splitcomplex2_equal(splitcomplex2 a, splitcomplex2 b)
{
    return vec4_all_eq(a, b);
}

static int
splitcomplex2_not_equal(splitcomplex2 a, splitcomplex2 b)
{
    return !vec4_all_eq(a, b);
}

static void
loop_splitcomplex2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    splitcomplex2_predicate op = (splitcomplex2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysplitcomplex2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_splitcomplex2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_splitcomplex2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysplitcomplex2_register_ufuncs(PyObject *numpy)
{
    return pysplitcomplex2_register_loop(numpy, "add", loop_splitcomplex2_binary, 3, 0,
                                    vec4_add) ||
        pysplitcomplex2_register_loop(numpy, "subtract", loop_splitcomplex2_binary, 3, 0,
                                 vec4_sub) ||
        pysplitcomplex2_register_loop(numpy, "multiply", loop_splitcomplex2_binary, 3, 0,
                                 splitcomplex2_mul) ||
        pysplitcomplex2_register_loop(numpy, "divide", loop_splitcomplex2_binary, 3, 0,
                                 splitcomplex2_div) ||
        pysplitcomplex2_register_loop(numpy, "positive", loop_splitcomplex2_unary, 2, 0,
                                 vec4_pos) ||
        pysplitcomplex2_register_loop(numpy, "negative", loop_splitcomplex2_unary, 2, 0,
                                 vec4_neg) ||
        pysplitcomplex2_register_loop(numpy, "power", loop_splitcomplex2_binary, 3, 0,
                                 splitcomplex2_pow) ||
        pysplitcomplex2_register_loop(numpy, "exp", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_exp) ||
        pysplitcomplex2_register_loop(numpy, "log", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_log) ||
        pysplitcomplex2_register_loop(numpy, "sqrt", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_sqrt) ||
        pysplitcomplex2_register_loop(numpy, "sin", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_sin) ||
        pysplitcomplex2_register_loop(numpy, "cos", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_cos) ||
        pysplitcomplex2_register_loop(numpy, "tan", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_tan) ||
        pysplitcomplex2_register_loop(numpy, "sinh", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_sinh) ||
        pysplitcomplex2_register_loop(numpy, "cosh", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_cosh) ||
        pysplitcomplex2_register_loop(numpy, "tanh", loop_splitcomplex2_unary, 2, 0,
                                 splitcomplex2_tanh) ||
        pysplitcomplex2_register_loop(numpy, "equal", loop_splitcomplex2_predicate, 3, 1,
                                 splitcomplex2_equal) ||
        pysplitcomplex2_register_loop(numpy, "not_equal", loop_splitcomplex2_predicate, 3, 1,
                                 splitcomplex2_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_splitcomplex2", "Bi-split-complex number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__splitcomplex2(void)
{
    PyObject *m, *numpy;
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
    PySplitcomplex2_Type = PyType_FromModuleAndSpec(
        m, &pysplitcomplex2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySplitcomplex2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysplitcomplex2_arrfuncs);
    pysplitcomplex2_arrfuncs.getitem = npysplitcomplex2_getitem;
    pysplitcomplex2_arrfuncs.setitem = npysplitcomplex2_setitem;
    pysplitcomplex2_arrfuncs.copyswapn = pysplitcomplex2_copyswapn;
    pysplitcomplex2_arrfuncs.copyswap = pysplitcomplex2_copyswap;
    pysplitcomplex2_arrfuncs.nonzero = pysplitcomplex2_nonzero;
    pysplitcomplex2_descr_proto.typeobj = (PyTypeObject *)PySplitcomplex2_Type;
    Py_SET_TYPE((PyObject *)&pysplitcomplex2_descr_proto, &PyArrayDescr_Type);
    npy_splitcomplex2 = PyArray_RegisterDataType(&pysplitcomplex2_descr_proto);
    if (npy_splitcomplex2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_splitcomplex2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySplitcomplex2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (pysplitcomplex2_register_casts() < 0 ||
            pysplitcomplex2_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "splitcomplex2", PySplitcomplex2_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
