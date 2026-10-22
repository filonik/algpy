/*
 * Multicomplex numbers x + x0*e0 + x1*e1 + x01*e0*e1 + x2*e2 + x02*e0*e2 +
 * x12*e1*e2 + x012*e0*e1*e2 (e0*e0 = k0, e1*e1 = k1, e2*e2 = k2) exposed to
 * Python as NumPy user-defined dtypes.
 *
 * The parameter k is not part of the data: `multicomplex3(k0, k1, k2)` creates a
 * new scalar type and NumPy dtype on every call (not cached; keep and reuse it).
 * The parameters live in the mutable class attribute `K` (a list). The arithmetic
 * itself lives in `multicomplex3.h`.
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

#include "../../norm.h"

#include "multicomplex3.h"

/* Per-k type information */

typedef struct {
    PyObject *type;
    int typenum;
    vec3 k;
} Info;

typedef struct {
    PyObject_HEAD
    multicomplex3 v;
} PyMC3;


static Info *
info_of_type(PyTypeObject *type)
{
    PyObject *cap = PyObject_GetAttrString((PyObject *)type, "_info");
    Info *info;
    if (cap == NULL) {
        return NULL;
    }
    info = (Info *)PyCapsule_GetPointer(cap, "algpy.multicomplex3.info");
    Py_DECREF(cap);
    return info;
}

static PyObject *
mc3_from_value(PyTypeObject *type, multicomplex3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyMC3 *p = (PyMC3 *)tp_alloc(type, 0);
    if (p) {
        p->v = v;
    }
    return (PyObject *)p;
}

/*
 * Convert a multicomplex3 or real number. Returns 1 on success, 0 if `o` is
 * not convertible and -1 with an error set. `info` is the expected type.
 */
static int
mc3_from_object(Info *info, PyObject *o, multicomplex3 *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)info->type)) {
        *out = ((PyMC3 *)o)->v;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = multicomplex3_make(d, 0, 0, 0, 0, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pymc3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12",
                             "x012", NULL};
    double c[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
                                     c, c + 1, c + 2, c + 3, c + 4, c + 5,
                                     c + 6, c + 7)) {
        return NULL;
    }
    return mc3_from_value(type, multicomplex3_make(c[0], c[1], c[2], c[3],
                                                   c[4], c[5], c[6], c[7]));
}

static PyObject *
pymc3_repr(PyObject *self)
{
    multicomplex3 v = ((PyMC3 *)self)->v;
    char *s[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("multicomplex3(%s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6],
                                  s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pymc3_hash(PyObject *self)
{
    multicomplex3 v = ((PyMC3 *)self)->v;
    PyObject *obj;
    Py_hash_t h;
    if (vec8_is_real(v)) {
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

/* Resolve the type info of whichever operand is a multicomplex3. */
static Info *
binary_info(PyObject *a, PyObject *b)
{
    PyObject *o = PyType_GetSlot(Py_TYPE(a), Py_tp_new) == (void *)pymc3_new
                  ? a : b;
    return info_of_type(Py_TYPE(o));
}

static PyObject *
pymc3_richcompare(PyObject *a, PyObject *b, int op)
{
    multicomplex3 x, y;
    Info *info;
    int r;
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    info = binary_info(a, b);
    if (info == NULL) {
        return NULL;
    }
    if ((r = mc3_from_object(info, a, &x)) <= 0 ||
            (r = mc3_from_object(info, b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return PyBool_FromLong(vec8_all_eq(x, y) == (op == Py_EQ));
}

#define MC3_BINARY(name, expr)                                             \
    static PyObject *                                                      \
    pymc3_##name(PyObject *a, PyObject *b)                                 \
    {                                                                      \
        multicomplex3 x, y;                                                \
        Info *info;                                                        \
        int r;                                                             \
        if ((info = binary_info(a, b)) == NULL) {                          \
            return NULL;                                                   \
        }                                                                  \
        if ((r = mc3_from_object(info, a, &x)) <= 0 ||                     \
                (r = mc3_from_object(info, b, &y)) <= 0) {                 \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return mc3_from_value((PyTypeObject *)info->type, expr);           \
    }

MC3_BINARY(add, vec8_add(x, y))
MC3_BINARY(sub, vec8_sub(x, y))
MC3_BINARY(mul, multicomplex3_mul(info->k, x, y))
MC3_BINARY(div, multicomplex3_div(info->k, x, y))

static PyObject *
pymc3_pos(PyObject *self)
{
    return mc3_from_value(Py_TYPE(self), vec8_pos(((PyMC3 *)self)->v));
}

static PyObject *
pymc3_neg(PyObject *self)
{
    return mc3_from_value(Py_TYPE(self), vec8_neg(((PyMC3 *)self)->v));
}

static int
pymc3_bool(PyObject *self)
{
    return vec8_any_nonzero(((PyMC3 *)self)->v);
}

static PyObject *
pymc3_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return mc3_from_value(Py_TYPE(self),
                          multicomplex3_conjugate(((PyMC3 *)self)->v));
}

static PyObject *
pymc3_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc3_from_value((PyTypeObject *)cls, vec8_zero());
}

static PyObject *
pymc3_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return mc3_from_value((PyTypeObject *)cls, multicomplex3_one());
}

static PyObject *rebuild_func = NULL;

static PyObject *
pymc3_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {9};
    PyObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex3_minpoly(info->k, ((PyMC3 *)self)->v,
                           (double *)PyArray_DATA((PyArrayObject *)arr));
    return arr;
}

static PyObject *
pymc3_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    npy_intp dims[1] = {8};
    complex ev[8];
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
    multicomplex3_eivals(info->k, ((PyMC3 *)self)->v, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 8; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pymc3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    multicomplex3 v = ((PyMC3 *)self)->v;
    if (info == NULL) {
        return NULL;
    }
    return Py_BuildValue("(O(ddddddddddd))", rebuild_func, info->k.data[0],
                         info->k.data[1], info->k.data[2], v.data[0],
                         v.data[1], v.data[2], v.data[3], v.data[4],
                         v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pymc3_basis(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[1] = {ALGPY_MULTICOMPLEX3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex3_basis((multicomplex3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc3_coeffs(PyObject *cls, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type((PyTypeObject *)cls);
    npy_intp dims[2] = {ALGPY_MULTICOMPLEX3_N, ALGPY_MULTICOMPLEX3_N};
    PyArrayObject *arr;
    if (info == NULL) {
        return NULL;
    }
    arr = (PyArrayObject *)PyArray_SimpleNew(2, dims, info->typenum);
    if (arr == NULL) {
        return NULL;
    }
    multicomplex3_coeffs(
        info->k, (multicomplex3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pymc3_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    Info *info = info_of_type(Py_TYPE(self));
    if (info == NULL) {
        return NULL;
    }
    return algpy_norm_to_py(complex_root_real(multicomplex3_norm8(info->k, ((PyMC3 *)self)->v), 8));
}

static PyMethodDef pymc3_methods[] = {
    {"zero", pymc3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pymc3_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pymc3_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pymc3_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm", pymc3_norm, METH_NOARGS, "Norm (a float, or a complex1 if imaginary)."},
    {"conjugate", pymc3_conjugate, METH_NOARGS, "Conjugate."},
    {"minpoly", pymc3_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial."},
    {"eivals", pymc3_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"__reduce__", pymc3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pymc3_length(PyObject *self)
{
    return 8;
}

static PyObject *
pymc3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyMC3 *)self)->v.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pymc3_subscript(PyObject *self, PyObject *key)
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
    return pymc3_item(self, i);
}

static int
pymc3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pymc3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pymc3_slots[] = {
    {Py_tp_dealloc, pymc3_dealloc},
    {Py_tp_traverse, pymc3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pymc3_repr},
    {Py_tp_str, pymc3_repr},
    {Py_tp_hash, pymc3_hash},
    {Py_tp_doc, (void *)"Multicomplex number with float64 components"},
    {Py_tp_richcompare, pymc3_richcompare},
    {Py_tp_methods, pymc3_methods},
    {Py_sq_length, pymc3_length},
    {Py_sq_item, pymc3_item},
    {Py_mp_subscript, pymc3_subscript},
    {Py_tp_new, pymc3_new},
    {Py_nb_add, pymc3_add},
    {Py_nb_subtract, pymc3_sub},
    {Py_nb_multiply, pymc3_mul},
    {Py_nb_true_divide, pymc3_div},
    {Py_nb_positive, pymc3_pos},
    {Py_nb_negative, pymc3_neg},
    {Py_nb_bool, pymc3_bool},
    {0, NULL}
};

static PyType_Spec pymc3_spec = {
    .name = "algpy.algebras.multicomplex3",
    .basicsize = sizeof(PyMC3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pymc3_slots,
};

/* NumPy dtype support */

static PyObject *
npymc3_getitem(void *data, void *arr)
{
    multicomplex3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    memcpy(&v, data, sizeof(v));
    return mc3_from_value(descr->typeobj, v);
}

static int
npymc3_setitem(PyObject *item, void *data, void *arr)
{
    multicomplex3 v;
    PyArray_Descr *descr = PyArray_DESCR((PyArrayObject *)arr);
    Info *info = info_of_type(descr->typeobj);
    int r;
    if (info == NULL) {
        return -1;
    }
    r = mc3_from_object(info, item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to multicomplex3");
        double d[8];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "multicomplex3 needs exactly 8 components");
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
        v = multicomplex3_make(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
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
npymc3_copyswapn(void *dst, npy_intp dstride, void *src, npy_intp sstride,
                 npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(multicomplex3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npymc3_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(multicomplex3));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npymc3_nonzero(void *data, void *arr)
{
    multicomplex3 v;
    memcpy(&v, data, sizeof(v));
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npymc3_arrfuncs;

typedef struct { char c; multicomplex3 v; } align_test;

/* Casts from the builtin real types (target type is found via `toarr`) */

#define MC3_CAST_FROM(name, type)                                          \
    static void                                                            \
    npycast_##name##_to_mc3(void *from, void *to, npy_intp n,              \
                            void *fromarr, void *toarr)                    \
    {                                                                      \
        const type *src = (const type *)from;                              \
        multicomplex3 *dst = (multicomplex3 *)to;                          \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = multicomplex3_make((double)src[i], 0, 0, 0, 0, 0, 0, 0);                \
        }                                                                  \
    }

#define MC3_REAL_TYPES(X)                                                  \
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

MC3_REAL_TYPES(MC3_CAST_FROM)

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
                npycast_##name##_to_mc3) < 0 ||                            \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    MC3_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops. The loop data is a heap-allocated Op, one per (type, ufunc). */

typedef struct Op {
    const vec3 *k;
    multicomplex3 (*unary)(multicomplex3);
    multicomplex3 (*binary)(vec3, multicomplex3, multicomplex3);
    int (*pred)(multicomplex3, multicomplex3);
} Op;

static multicomplex3 op_add(vec3 k, multicomplex3 a, multicomplex3 b)
{ return vec8_add(a, b); }
static multicomplex3 op_sub(vec3 k, multicomplex3 a, multicomplex3 b)
{ return vec8_sub(a, b); }
static int op_not_equal(multicomplex3 a, multicomplex3 b)
{ return !vec8_all_eq(a, b); }

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    Op *op = (Op *)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        multicomplex3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op->unary(a);
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
        multicomplex3 a, b, r;
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
        multicomplex3 a, b;
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
    Op none = {NULL, NULL, NULL, NULL};
    Op o;
#define BIN(name, f) (o = none, o.binary = f, register_loop(info, numpy, name, loop_binary, 3, 0, o))
#define UN(name, f) (o = none, o.unary = f, register_loop(info, numpy, name, loop_unary, 2, 0, o))
#define PRED(name, f) (o = none, o.pred = f, register_loop(info, numpy, name, loop_predicate, 3, 1, o))
    return
        BIN("add", op_add) || BIN("subtract", op_sub) ||
        BIN("multiply", multicomplex3_mul) ||
        BIN("divide", multicomplex3_div) ||
        UN("positive", vec8_pos) || UN("negative", vec8_neg) ||
        PRED("equal", vec8_all_eq) || PRED("not_equal", op_not_equal)
        ? -1 : 0;
#undef BIN
#undef UN
#undef PRED
}

static int
register_conversions(PyObject *type, int n)
{
    PyObject *module, *attach, *result, *size;
    module = PyImport_ImportModule("algpy.algebras._conversions");
    if (module == NULL) {
        return -1;
    }
    attach = PyObject_GetAttrString(module, "attach");
    Py_DECREF(module);
    if (attach == NULL) {
        return -1;
    }
    size = PyLong_FromLong(n);
    if (size == NULL) {
        Py_DECREF(attach);
        return -1;
    }
    result = PyObject_CallFunctionObjArgs(attach, type, size, NULL);
    Py_DECREF(size);
    Py_DECREF(attach);
    if (result == NULL) {
        return -1;
    }
    Py_DECREF(result);
    return 0;
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
    .tp_name = "algpy.multicomplex3.type",
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

    type = PyType_FromSpecWithBases(&pymc3_spec,
                                    (PyObject *)&PyGenericArrType_Type);
    if (type == NULL) {
        goto fail;
    }
    Py_SET_TYPE(type, &MetaType);
    Py_INCREF(type);  /* kept alive by the registered dtype */
    info->type = type;
    cap = PyCapsule_New(info, "algpy.multicomplex3.info", NULL);
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
    proto->elsize = sizeof(multicomplex3);
    proto->alignment = offsetof(align_test, v);
    proto->f = &npymc3_arrfuncs;
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

    if (register_conversions(type, 8) < 0) {
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
mc3_factory(PyObject *module, PyObject *args)
{
    double k0, k1, k2;
    if (!PyArg_ParseTuple(args, "ddd", &k0, &k1, &k2)) {
        return NULL;
    }
    return create_type(k0, k1, k2);
}

static PyObject *
mc3_rebuild(PyObject *module, PyObject *args)
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
    res = mc3_from_value((PyTypeObject *)type, multicomplex3_make(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]));
    Py_DECREF(type);
    return res;
}

static PyMethodDef module_methods[] = {
    {"multicomplex3", mc3_factory, METH_VARARGS,
     "multicomplex3(k0, k1, k2): scalar type / dtype of numbers with e0*e0 = k0, e1*e1 = k1, e2*e2 = k2."},
    {"_rebuild", mc3_rebuild, METH_VARARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_multicomplex3",
    "Multicomplex number types for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__multicomplex3(void)
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
    PyArray_InitArrFuncs(&npymc3_arrfuncs);
    npymc3_arrfuncs.getitem = npymc3_getitem;
    npymc3_arrfuncs.setitem = npymc3_setitem;
    npymc3_arrfuncs.copyswapn = npymc3_copyswapn;
    npymc3_arrfuncs.copyswap = npymc3_copyswap;
    npymc3_arrfuncs.nonzero = npymc3_nonzero;

    rebuild_func = PyObject_GetAttrString(m, "_rebuild");
    if (rebuild_func == NULL) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}
