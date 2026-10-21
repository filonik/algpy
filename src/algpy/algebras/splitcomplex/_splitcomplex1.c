/*
 * Split-complex numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in splitcomplex1.h.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "splitcomplex1.h"

typedef struct {
    PyObject_HEAD
    splitcomplex1 v;
} PySplitcomplex1;

static PyObject *PySplitcomplex1_Type = NULL;
static int npy_splitcomplex1 = -1;

static PyObject *
PySplitcomplex1_FromValue(PyTypeObject *type, splitcomplex1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PySplitcomplex1 *p = (PySplitcomplex1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_SPLITCOMPLEX1(v) \
    PySplitcomplex1_FromValue((PyTypeObject *)PySplitcomplex1_Type, (v))

static int
splitcomplex1_from_object(PyObject *o, splitcomplex1 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PySplitcomplex1_Type)) {
        *out = ((PySplitcomplex1 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = splitcomplex1_make(d, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pysplitcomplex1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"real", "imag", NULL};
    double re = 0, im = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist, &re, &im)) {
        return NULL;
    }
    return PySplitcomplex1_FromValue(type, splitcomplex1_make(re, im));
}

static PyObject *
pysplitcomplex1_repr(PyObject *self)
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    char *re = PyOS_double_to_string(v.data[0], 'g', 6, 0, NULL);
    char *im = PyOS_double_to_string(v.data[1], 'g', 6, 0, NULL);
    PyObject *result;
    if (re == NULL || im == NULL) {
        PyMem_Free(re);
        PyMem_Free(im);
        return NULL;
    }
    result = PyUnicode_FromFormat("splitcomplex1(%s, %s)", re, im);
    PyMem_Free(re);
    PyMem_Free(im);
    return result;
}

static Py_hash_t
pysplitcomplex1_hash(PyObject *self)
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    if (v.data[1] == 0) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dd)", v.data[0], v.data[1]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pysplitcomplex1_richcompare(PyObject *a, PyObject *b, int op)
{
    splitcomplex1 x, y;
    int r;
    if ((r = splitcomplex1_from_object(a, &x)) <= 0 ||
            (r = splitcomplex1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(x.data[0] == y.data[0] &&
                                   x.data[1] == y.data[1]);
        case Py_NE:
            return PyBool_FromLong(x.data[0] != y.data[0] ||
                                   x.data[1] != y.data[1]);
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define SPLITCOMPLEX1_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pysplitcomplex1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        splitcomplex1 x, y;                                                       \
        int r;                                                               \
        if ((r = splitcomplex1_from_object(a, &x)) <= 0 ||                        \
                (r = splitcomplex1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_SPLITCOMPLEX1(expr(x, y));                                     \
    }

SPLITCOMPLEX1_BINARY(add, vec2_add)
SPLITCOMPLEX1_BINARY(sub, vec2_sub)
SPLITCOMPLEX1_BINARY(mul, splitcomplex1_mul)
SPLITCOMPLEX1_BINARY(div, splitcomplex1_div)

static PyObject *
pysplitcomplex1_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    splitcomplex1 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "splitcomplex1 does not support modular power");
        return NULL;
    }
    if ((r = splitcomplex1_from_object(a, &x)) <= 0 ||
            (r = splitcomplex1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_SPLITCOMPLEX1(splitcomplex1_pow(x, y));
}

#define SPLITCOMPLEX1_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pysplitcomplex1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_SPLITCOMPLEX1(expr(((PySplitcomplex1 *)self)->v));                  \
    }

SPLITCOMPLEX1_UNARY(pos, vec2_pos)
SPLITCOMPLEX1_UNARY(neg, vec2_neg)

static PyObject *
pysplitcomplex1_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex1_FromValue((PyTypeObject *)cls, vec2_zero());
}

static PyObject *
pysplitcomplex1_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PySplitcomplex1_FromValue((PyTypeObject *)cls, splitcomplex1_one());
}

static int
pysplitcomplex1_bool(PyObject *self)
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    return v.data[0] != 0 || v.data[1] != 0;
}

#define SPLITCOMPLEX1_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pysplitcomplex1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_SPLITCOMPLEX1(op(((PySplitcomplex1 *)self)->v));                   \
    }

SPLITCOMPLEX1_METHOD(exp, splitcomplex1_exp)
SPLITCOMPLEX1_METHOD(log, splitcomplex1_log)
SPLITCOMPLEX1_METHOD(sin, splitcomplex1_sin)
SPLITCOMPLEX1_METHOD(cos, splitcomplex1_cos)
SPLITCOMPLEX1_METHOD(tan, splitcomplex1_tan)
SPLITCOMPLEX1_METHOD(sinh, splitcomplex1_sinh)
SPLITCOMPLEX1_METHOD(cosh, splitcomplex1_cosh)
SPLITCOMPLEX1_METHOD(tanh, splitcomplex1_tanh)
SPLITCOMPLEX1_METHOD(sqrt, splitcomplex1_sqrt)
SPLITCOMPLEX1_METHOD(conjugate, splitcomplex1_conjugate)

static PyObject *
pysplitcomplex1_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(splitcomplex1_norm2(((PySplitcomplex1 *)self)->v));
}

static PyObject *
pysplitcomplex1_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double *d;
    if (arr == NULL) {
        return NULL;
    }
    d = (double *)PyArray_DATA((PyArrayObject *)arr);
    d[0] = splitcomplex1_norm2(v);
    d[1] = -2 * v.data[0];
    d[2] = 1;
    return arr;
}

static PyObject *
pysplitcomplex1_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    npy_intp dims[1] = {2};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    npy_csetreal(&d[0], v.data[0] - v.data[1]);
    npy_csetimag(&d[0], 0);
    npy_csetreal(&d[1], v.data[0] + v.data[1]);
    npy_csetimag(&d[1], 0);
    return arr;
}

static PyObject *
pysplitcomplex1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_splitcomplex1);
    splitcomplex1 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (splitcomplex1 *)PyArray_DATA(arr);
    splitcomplex1_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex1_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {2, 2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_splitcomplex1);
    splitcomplex1 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (splitcomplex1 *)PyArray_DATA(arr);
    splitcomplex1_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pysplitcomplex1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    splitcomplex1 v = ((PySplitcomplex1 *)self)->v;
    return Py_BuildValue("(O(dd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1]);
}

static PyMethodDef pysplitcomplex1_methods[] = {
    {"zero", pysplitcomplex1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pysplitcomplex1_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pysplitcomplex1_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pysplitcomplex1_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pysplitcomplex1_norm2, METH_NOARGS, "Squared norm."},
    {"conjugate", pysplitcomplex1_conjugate, METH_NOARGS, "Split-complex conjugate."},
    {"minpoly", pysplitcomplex1_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pysplitcomplex1_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pysplitcomplex1_exp, METH_NOARGS, "Split-complex exponential."},
    {"log", pysplitcomplex1_log, METH_NOARGS, "Principal split-complex logarithm."},
    {"sin", pysplitcomplex1_sin, METH_NOARGS, "Split-complex sine."},
    {"cos", pysplitcomplex1_cos, METH_NOARGS, "Split-complex cosine."},
    {"tan", pysplitcomplex1_tan, METH_NOARGS, "Split-complex tangent."},
    {"sinh", pysplitcomplex1_sinh, METH_NOARGS, "Split-complex hyperbolic sine."},
    {"cosh", pysplitcomplex1_cosh, METH_NOARGS, "Split-complex hyperbolic cosine."},
    {"tanh", pysplitcomplex1_tanh, METH_NOARGS, "Split-complex hyperbolic tangent."},
    {"sqrt", pysplitcomplex1_sqrt, METH_NOARGS, "Principal split-complex square root."},
    {"__reduce__", pysplitcomplex1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pysplitcomplex1_length(PyObject *Py_UNUSED(self))
{
    return 2;
}

static PyObject *
pysplitcomplex1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PySplitcomplex1 *)self)->v.data[i]);
}

static PyObject *
pysplitcomplex1_subscript(PyObject *self, PyObject *key)
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
    return pysplitcomplex1_item(self, i);
}

static int
pysplitcomplex1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pysplitcomplex1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pysplitcomplex1_slots[] = {
    {Py_tp_dealloc, pysplitcomplex1_dealloc},
    {Py_tp_traverse, pysplitcomplex1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pysplitcomplex1_repr},
    {Py_tp_str, pysplitcomplex1_repr},
    {Py_tp_hash, pysplitcomplex1_hash},
    {Py_tp_doc, (void *)"Split-complex number with float64 components"},
    {Py_tp_richcompare, pysplitcomplex1_richcompare},
    {Py_tp_methods, pysplitcomplex1_methods},
    {Py_sq_length, pysplitcomplex1_length},
    {Py_sq_item, pysplitcomplex1_item},
    {Py_mp_subscript, pysplitcomplex1_subscript},
    {Py_tp_new, pysplitcomplex1_new},
    {Py_nb_add, pysplitcomplex1_add},
    {Py_nb_subtract, pysplitcomplex1_sub},
    {Py_nb_multiply, pysplitcomplex1_mul},
    {Py_nb_true_divide, pysplitcomplex1_div},
    {Py_nb_power, (void *)pysplitcomplex1_pow},
    {Py_nb_positive, pysplitcomplex1_pos},
    {Py_nb_negative, pysplitcomplex1_neg},
    {Py_nb_bool, pysplitcomplex1_bool},
    {0, NULL}
};

static PyType_Spec pysplitcomplex1_spec = {
    .name = "algpy.algebras.splitcomplex.splitcomplex1",
    .basicsize = sizeof(PySplitcomplex1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pysplitcomplex1_slots,
};

static PyObject *
npysplitcomplex1_getitem(void *data, void *arr)
{
    splitcomplex1 v;
    memcpy(&v, data, sizeof(v));
    return NEW_SPLITCOMPLEX1(v);
}

static int
npysplitcomplex1_setitem(PyObject *item, void *data, void *arr)
{
    splitcomplex1 v;
    int r = splitcomplex1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to splitcomplex1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "splitcomplex1 needs exactly 2 components");
            return -1;
        }
        v.data[0] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 0));
        v.data[1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, 1));
        Py_DECREF(seq);
        if (PyErr_Occurred()) {
            return -1;
        }
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pysplitcomplex1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(splitcomplex1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            char *p = (char *)dst + i * dstride;
            int j;
            for (j = 0; j < 2; j++) {
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
pysplitcomplex1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pysplitcomplex1_copyswapn(dst, sizeof(splitcomplex1), src, sizeof(splitcomplex1), 1,
                         swap, arr);
}

static npy_bool
pysplitcomplex1_nonzero(void *data, void *arr)
{
    splitcomplex1 v;
    memcpy(&v, data, sizeof(v));
    return (v.data[0] != 0 || v.data[1] != 0) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pysplitcomplex1_arrfuncs;
typedef struct { char c; splitcomplex1 v; } pysplitcomplex1_align_test;

static PyArray_DescrProto pysplitcomplex1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(splitcomplex1), offsetof(pysplitcomplex1_align_test, v), 0, 0, 0,
    &pysplitcomplex1_arrfuncs,
};

#define SPLITCOMPLEX1_CAST_FROM(name, type)                                       \
    static void                                                              \
    pysplitcomplex1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        splitcomplex1 *dst = (splitcomplex1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = splitcomplex1_make((double)src[i], 0);                      \
        }                                                                    \
    }

#define SPLITCOMPLEX1_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

SPLITCOMPLEX1_REAL_TYPES(SPLITCOMPLEX1_CAST_FROM)

static int
pysplitcomplex1_register_casts(void)
{
#define REGISTER_SPLITCOMPLEX1_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_splitcomplex1,                    \
                                     pysplitcomplex1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_splitcomplex1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    SPLITCOMPLEX1_REAL_TYPES(REGISTER_SPLITCOMPLEX1_CAST)
#undef REGISTER_SPLITCOMPLEX1_CAST
    return 0;
}

typedef splitcomplex1 (*splitcomplex1_unary_op)(splitcomplex1);
typedef splitcomplex1 (*splitcomplex1_binary_op)(splitcomplex1, splitcomplex1);
typedef int (*splitcomplex1_predicate)(splitcomplex1, splitcomplex1);

static void
loop_splitcomplex1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    splitcomplex1_unary_op op = (splitcomplex1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_splitcomplex1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    splitcomplex1_binary_op op = (splitcomplex1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex1 a, b, r;
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
splitcomplex1_equal(splitcomplex1 a, splitcomplex1 b)
{
    return a.data[0] == b.data[0] && a.data[1] == b.data[1];
}

static int
splitcomplex1_not_equal(splitcomplex1 a, splitcomplex1 b)
{
    return a.data[0] != b.data[0] || a.data[1] != b.data[1];
}

static void
loop_splitcomplex1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    splitcomplex1_predicate op = (splitcomplex1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        splitcomplex1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pysplitcomplex1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_splitcomplex1;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_splitcomplex1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pysplitcomplex1_register_ufuncs(PyObject *numpy)
{
    return pysplitcomplex1_register_loop(numpy, "add", loop_splitcomplex1_binary, 3, 0,
                                    vec2_add) ||
        pysplitcomplex1_register_loop(numpy, "subtract", loop_splitcomplex1_binary, 3, 0,
                                 vec2_sub) ||
        pysplitcomplex1_register_loop(numpy, "multiply", loop_splitcomplex1_binary, 3, 0,
                                 splitcomplex1_mul) ||
        pysplitcomplex1_register_loop(numpy, "divide", loop_splitcomplex1_binary, 3, 0,
                                 splitcomplex1_div) ||
        pysplitcomplex1_register_loop(numpy, "power", loop_splitcomplex1_binary, 3, 0,
                                 splitcomplex1_pow) ||
        pysplitcomplex1_register_loop(numpy, "positive", loop_splitcomplex1_unary, 2, 0,
                                 vec2_pos) ||
        pysplitcomplex1_register_loop(numpy, "negative", loop_splitcomplex1_unary, 2, 0,
                                 vec2_neg) ||
        pysplitcomplex1_register_loop(numpy, "exp", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_exp) ||
        pysplitcomplex1_register_loop(numpy, "log", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_log) ||
        pysplitcomplex1_register_loop(numpy, "sin", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_sin) ||
        pysplitcomplex1_register_loop(numpy, "cos", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_cos) ||
        pysplitcomplex1_register_loop(numpy, "tan", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_tan) ||
        pysplitcomplex1_register_loop(numpy, "sinh", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_sinh) ||
        pysplitcomplex1_register_loop(numpy, "cosh", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_cosh) ||
        pysplitcomplex1_register_loop(numpy, "tanh", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_tanh) ||
        pysplitcomplex1_register_loop(numpy, "sqrt", loop_splitcomplex1_unary, 2, 0,
                                 splitcomplex1_sqrt) ||
        pysplitcomplex1_register_loop(numpy, "equal", loop_splitcomplex1_predicate, 3, 1,
                                 splitcomplex1_equal) ||
        pysplitcomplex1_register_loop(numpy, "not_equal", loop_splitcomplex1_predicate, 3, 1,
                                 splitcomplex1_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_splitcomplex1", "Split-complex number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__splitcomplex1(void)
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
    PySplitcomplex1_Type = PyType_FromModuleAndSpec(
        m, &pysplitcomplex1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PySplitcomplex1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pysplitcomplex1_arrfuncs);
    pysplitcomplex1_arrfuncs.getitem = npysplitcomplex1_getitem;
    pysplitcomplex1_arrfuncs.setitem = npysplitcomplex1_setitem;
    pysplitcomplex1_arrfuncs.copyswapn = pysplitcomplex1_copyswapn;
    pysplitcomplex1_arrfuncs.copyswap = pysplitcomplex1_copyswap;
    pysplitcomplex1_arrfuncs.nonzero = pysplitcomplex1_nonzero;
    pysplitcomplex1_descr_proto.typeobj = (PyTypeObject *)PySplitcomplex1_Type;
    Py_SET_TYPE((PyObject *)&pysplitcomplex1_descr_proto, &PyArrayDescr_Type);
    npy_splitcomplex1 = PyArray_RegisterDataType(&pysplitcomplex1_descr_proto);
    if (npy_splitcomplex1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_splitcomplex1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PySplitcomplex1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (pysplitcomplex1_register_casts() < 0 ||
            pysplitcomplex1_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "splitcomplex1", PySplitcomplex1_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
