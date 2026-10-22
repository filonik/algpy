/*
 * Bi-complex numbers exposed to Python as a NumPy user-defined dtype.
 * The arithmetic itself lives in complex2.h.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "complex2.h"
#include "../norm.h"

typedef struct {
    PyObject_HEAD
    complex2 v;
} PyComplex2;

static PyObject *PyComplex2_Type = NULL;
static int npy_complex2 = -1;

static PyObject *
PyComplex2_FromValue(PyTypeObject *type, complex2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyComplex2 *p = (PyComplex2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_COMPLEX2(v) \
    PyComplex2_FromValue((PyTypeObject *)PyComplex2_Type, (v))

static int
complex2_from_object(PyObject *o, complex2 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyComplex2_Type)) {
        *out = ((PyComplex2 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = complex2_make(d, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pycomplex2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"w", "x", "y", "z", NULL};
    double w = 0, x = 0, y = 0, z = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                     &w, &x, &y, &z)) {
        return NULL;
    }
    return PyComplex2_FromValue(type, complex2_make(w, x, y, z));
}

static PyObject *
pycomplex2_repr(PyObject *self)
{
    complex2 v = ((PyComplex2 *)self)->v;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("complex2(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pycomplex2_hash(PyObject *self)
{
    complex2 v = ((PyComplex2 *)self)->v;
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
pycomplex2_richcompare(PyObject *a, PyObject *b, int op)
{
    complex2 x, y;
    int r;
    if ((r = complex2_from_object(a, &x)) <= 0 ||
            (r = complex2_from_object(b, &y)) <= 0) {
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

#define COMPLEX2_BINARY(name, expr)                                          \
    static PyObject *                                                        \
    pycomplex2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        complex2 x, y;                                                       \
        int r;                                                               \
        if ((r = complex2_from_object(a, &x)) <= 0 ||                        \
                (r = complex2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_COMPLEX2(expr(x, y));                                     \
    }

COMPLEX2_BINARY(add, vec4_add)
COMPLEX2_BINARY(sub, vec4_sub)
COMPLEX2_BINARY(mul, complex2_mul)
COMPLEX2_BINARY(div, complex2_div)

static PyObject *
pycomplex2_pow(PyObject *a, PyObject *b, PyObject *mod)
{
    complex2 x, y;
    int r;
    if (mod != Py_None) {
        PyErr_SetString(PyExc_TypeError, "complex2 does not support modular power");
        return NULL;
    }
    if ((r = complex2_from_object(a, &x)) <= 0 ||
            (r = complex2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_COMPLEX2(complex2_pow(x, y));
}

#define COMPLEX2_UNARY(name, expr)                                           \
    static PyObject *                                                        \
    pycomplex2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_COMPLEX2(expr(((PyComplex2 *)self)->v));                  \
    }

COMPLEX2_UNARY(pos, vec4_pos)
COMPLEX2_UNARY(neg, vec4_neg)

static PyObject *
pycomplex2_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex2_FromValue((PyTypeObject *)cls, vec4_zero());
}

static PyObject *
pycomplex2_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyComplex2_FromValue((PyTypeObject *)cls, complex2_one());
}

static int
pycomplex2_bool(PyObject *self)
{
    complex2 v = ((PyComplex2 *)self)->v;
    return vec4_any_nonzero(v);
}

#define COMPLEX2_METHOD(name, op)                                            \
    static PyObject *                                                        \
    pycomplex2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return NEW_COMPLEX2(op(((PyComplex2 *)self)->v));                   \
    }

COMPLEX2_METHOD(conjugate, complex2_conjugate)
COMPLEX2_METHOD(exp, complex2_exp)
COMPLEX2_METHOD(log, complex2_log)
COMPLEX2_METHOD(sqrt, complex2_sqrt)

static PyObject *
pycomplex2_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex n = complex2_norm2(((PyComplex2 *)self)->v);
    return Py_BuildValue("(dd)", n.data[0], n.data[1]);
}

static PyObject *
pycomplex2_norm4(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(complex2_norm4(((PyComplex2 *)self)->v));
}

static PyObject *
pycomplex2_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return algpy_norm_to_py(complex2_norm(((PyComplex2 *)self)->v));
}

static PyObject *
pycomplex2_inv(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_COMPLEX2(complex2_inv(((PyComplex2 *)self)->v));
}

static PyObject *
pycomplex2_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    complex2_minpoly(((PyComplex2 *)self)->v,
                     (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pycomplex2_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex e[4];
    npy_intp dims[1] = {4};
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    complex2_eivals(((PyComplex2 *)self)->v, e);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 4; i++) {
        npy_csetreal(&d[i], e[i].data[0]);
        npy_csetimag(&d[i], e[i].data[1]);
    }
    return arr;
}

static PyObject *
pycomplex2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_complex2);
    complex2 *basis;
    if (arr == NULL) {
        return NULL;
    }
    basis = (complex2 *)PyArray_DATA(arr);
    complex2_basis(basis);
    return (PyObject *)arr;
}

static PyObject *
pycomplex2_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {4, 4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_complex2);
    complex2 *coeffs;
    if (arr == NULL) {
        return NULL;
    }
    coeffs = (complex2 *)PyArray_DATA(arr);
    complex2_coeffs(coeffs);
    return (PyObject *)arr;
}

static PyObject *
pycomplex2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    complex2 v = ((PyComplex2 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyMethodDef pycomplex2_methods[] = {
    {"zero", pycomplex2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pycomplex2_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pycomplex2_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pycomplex2_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pycomplex2_norm2, METH_NOARGS, "Squared norm (a complex pair)."},
    {"norm4", pycomplex2_norm4, METH_NOARGS, "Fourth power of the norm."},
    {"norm", pycomplex2_norm, METH_NOARGS, "Norm."},
    {"inv", pycomplex2_inv, METH_NOARGS, "Multiplicative inverse."},
    {"conjugate", pycomplex2_conjugate, METH_NOARGS, "Conjugate of the i1 part."},
    {"minpoly", pycomplex2_minpoly, METH_NOARGS, "Minimal polynomial coefficients."},
    {"eivals", pycomplex2_eivals, METH_NOARGS, "Eigenvalues."},
    {"exp", pycomplex2_exp, METH_NOARGS, "exp."},
    {"log", pycomplex2_log, METH_NOARGS, "log."},
    {"sqrt", pycomplex2_sqrt, METH_NOARGS, "sqrt."},
    {"__reduce__", pycomplex2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pycomplex2_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pycomplex2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyComplex2 *)self)->v.data[i]);
}

static PyObject *
pycomplex2_subscript(PyObject *self, PyObject *key)
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
    return pycomplex2_item(self, i);
}

static int
pycomplex2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pycomplex2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pycomplex2_slots[] = {
    {Py_tp_dealloc, pycomplex2_dealloc},
    {Py_tp_traverse, pycomplex2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pycomplex2_repr},
    {Py_tp_str, pycomplex2_repr},
    {Py_tp_hash, pycomplex2_hash},
    {Py_tp_doc, (void *)"Bi-complex number with float64 components"},
    {Py_tp_richcompare, pycomplex2_richcompare},
    {Py_tp_methods, pycomplex2_methods},
    {Py_sq_length, pycomplex2_length},
    {Py_sq_item, pycomplex2_item},
    {Py_mp_subscript, pycomplex2_subscript},
    {Py_tp_new, pycomplex2_new},
    {Py_nb_add, pycomplex2_add},
    {Py_nb_subtract, pycomplex2_sub},
    {Py_nb_multiply, pycomplex2_mul},
    {Py_nb_true_divide, pycomplex2_div},
    {Py_nb_power, (void *)pycomplex2_pow},
    {Py_nb_positive, pycomplex2_pos},
    {Py_nb_negative, pycomplex2_neg},
    {Py_nb_bool, pycomplex2_bool},
    {0, NULL}
};

static PyType_Spec pycomplex2_spec = {
    .name = "algpy.algebras.multicomplex.complex2",
    .basicsize = sizeof(PyComplex2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pycomplex2_slots,
};

static PyObject *
npycomplex2_getitem(void *data, void *arr)
{
    complex2 v;
    memcpy(&v, data, sizeof(v));
    return NEW_COMPLEX2(v);
}

static int
npycomplex2_setitem(PyObject *item, void *data, void *arr)
{
    complex2 v;
    int r = complex2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to complex2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "complex2 needs exactly 4 components");
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
pycomplex2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(complex2));
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
pycomplex2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pycomplex2_copyswapn(dst, sizeof(complex2), src, sizeof(complex2), 1,
                         swap, arr);
}

static npy_bool
pycomplex2_nonzero(void *data, void *arr)
{
    complex2 v;
    memcpy(&v, data, sizeof(v));
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pycomplex2_arrfuncs;
typedef struct { char c; complex2 v; } pycomplex2_align_test;

static PyArray_DescrProto pycomplex2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(complex2), offsetof(pycomplex2_align_test, v), 0, 0, 0,
    &pycomplex2_arrfuncs,
};

#define COMPLEX2_CAST_FROM(name, type)                                       \
    static void                                                              \
    pycomplex2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        complex2 *dst = (complex2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = complex2_make((double)src[i], 0, 0, 0);                      \
        }                                                                    \
    }

#define COMPLEX2_REAL_TYPES(X)                                               \
    X(BOOL, npy_bool) X(BYTE, npy_byte) X(UBYTE, npy_ubyte)                 \
    X(SHORT, npy_short) X(USHORT, npy_ushort) X(INT, npy_int)              \
    X(UINT, npy_uint) X(LONG, npy_long) X(ULONG, npy_ulong)                 \
    X(LONGLONG, npy_longlong) X(ULONGLONG, npy_ulonglong)                   \
    X(FLOAT, npy_float) X(DOUBLE, npy_double)

COMPLEX2_REAL_TYPES(COMPLEX2_CAST_FROM)

static int
pycomplex2_register_casts(void)
{
#define REGISTER_COMPLEX2_CAST(name, type)                                   \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_complex2,                    \
                                     pycomplex2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_complex2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    COMPLEX2_REAL_TYPES(REGISTER_COMPLEX2_CAST)
#undef REGISTER_COMPLEX2_CAST
    return 0;
}

typedef complex2 (*complex2_unary_op)(complex2);
typedef complex2 (*complex2_binary_op)(complex2, complex2);
typedef int (*complex2_predicate)(complex2, complex2);

static void
loop_complex2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    complex2_unary_op op = (complex2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_complex2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    complex2_binary_op op = (complex2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex2 a, b, r;
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
complex2_equal(complex2 a, complex2 b)
{
    return vec4_all_eq(a, b);
}

static int
complex2_not_equal(complex2 a, complex2 b)
{
    return !vec4_all_eq(a, b);
}

static void
loop_complex2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    complex2_predicate op = (complex2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        complex2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
pycomplex2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs, int out_bool,
                         void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_complex2;
    }
    if (out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_complex2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pycomplex2_register_ufuncs(PyObject *numpy)
{
    return pycomplex2_register_loop(numpy, "add", loop_complex2_binary, 3, 0,
                                    vec4_add) ||
        pycomplex2_register_loop(numpy, "subtract", loop_complex2_binary, 3, 0,
                                 vec4_sub) ||
        pycomplex2_register_loop(numpy, "multiply", loop_complex2_binary, 3, 0,
                                 complex2_mul) ||
        pycomplex2_register_loop(numpy, "divide", loop_complex2_binary, 3, 0,
                                 complex2_div) ||
        pycomplex2_register_loop(numpy, "positive", loop_complex2_unary, 2, 0,
                                 vec4_pos) ||
        pycomplex2_register_loop(numpy, "negative", loop_complex2_unary, 2, 0,
                                 vec4_neg) ||
        pycomplex2_register_loop(numpy, "power", loop_complex2_binary, 3, 0,
                                 complex2_pow) ||
        pycomplex2_register_loop(numpy, "exp", loop_complex2_unary, 2, 0,
                                 complex2_exp) ||
        pycomplex2_register_loop(numpy, "log", loop_complex2_unary, 2, 0,
                                 complex2_log) ||
        pycomplex2_register_loop(numpy, "sqrt", loop_complex2_unary, 2, 0,
                                 complex2_sqrt) ||
        pycomplex2_register_loop(numpy, "equal", loop_complex2_predicate, 3, 1,
                                 complex2_equal) ||
        pycomplex2_register_loop(numpy, "not_equal", loop_complex2_predicate, 3, 1,
                                 complex2_not_equal)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_complex2", "Bi-complex number type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__complex2(void)
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
    PyComplex2_Type = PyType_FromModuleAndSpec(
        m, &pycomplex2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyComplex2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pycomplex2_arrfuncs);
    pycomplex2_arrfuncs.getitem = npycomplex2_getitem;
    pycomplex2_arrfuncs.setitem = npycomplex2_setitem;
    pycomplex2_arrfuncs.copyswapn = pycomplex2_copyswapn;
    pycomplex2_arrfuncs.copyswap = pycomplex2_copyswap;
    pycomplex2_arrfuncs.nonzero = pycomplex2_nonzero;
    pycomplex2_descr_proto.typeobj = (PyTypeObject *)PyComplex2_Type;
    Py_SET_TYPE((PyObject *)&pycomplex2_descr_proto, &PyArrayDescr_Type);
    npy_complex2 = PyArray_RegisterDataType(&pycomplex2_descr_proto);
    if (npy_complex2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_complex2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyComplex2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (pycomplex2_register_casts() < 0 ||
            pycomplex2_register_ufuncs(numpy) < 0 ||
            PyModule_AddObjectRef(m, "complex2", PyComplex2_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;
fail:
    Py_XDECREF(numpy);
    Py_DECREF(m);
    return NULL;
}
