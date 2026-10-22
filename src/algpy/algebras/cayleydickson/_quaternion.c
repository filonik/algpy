/*
 * Quaternions exposed to Python as a NumPy user-defined dtype.
 *
 * This file only contains the Python/NumPy boilerplate, modelled on
 * `_rational_tests.c` from the NumPy repository. The arithmetic itself lives
 * in `quaternion.h`.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stddef.h>
#include <string.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

#include "quaternion.h"
#include "../norm.h"

/* Python scalar type */

typedef struct {
    PyObject_HEAD
    quaternion q;
} PyQuaternion;

static PyObject *PyQuaternion_Type = NULL;
static int npy_quaternion = -1;

static PyObject *
PyQuaternion_FromQuaternion(PyTypeObject *type, quaternion q)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyQuaternion *p = (PyQuaternion *)tp_alloc(type, 0);
    if (p) {
        p->q = q;
    }
    return (PyObject *)p;
}

#define NEW_QUATERNION(q) \
    PyQuaternion_FromQuaternion((PyTypeObject *)PyQuaternion_Type, (q))

/*
 * Convert a quaternion or real number to a quaternion.
 * Returns 1 on success, 0 if `o` is not convertible (no error set) and
 * -1 with an error set on failure.
 */
static int
quaternion_from_object(PyObject *o, quaternion *out)
{
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyQuaternion_Type)) {
        *out = ((PyQuaternion *)o)->q;
        return 1;
    }
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            /* e.g. OverflowError for huge ints */
            return -1;
        }
        *out = quaternion_make(d, 0, 0, 0);
        return 1;
    }
    return 0;
}

static PyObject *
pyquaternion_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"w", "x", "y", "z", NULL};
    double w = 0, x = 0, y = 0, z = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                     &w, &x, &y, &z)) {
        return NULL;
    }
    return PyQuaternion_FromQuaternion(type, quaternion_make(w, x, y, z));
}

static PyObject *
pyquaternion_repr(PyObject *self)
{
    quaternion q = ((PyQuaternion *)self)->q;
    const double *v = q.data;
    char *s[4] = {NULL, NULL, NULL, NULL};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("quaternion(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyquaternion_hash(PyObject *self)
{
    quaternion q = ((PyQuaternion *)self)->q;
    PyObject *tuple;
    Py_hash_t h;
    if (q.data[1] == 0 && q.data[2] == 0 && q.data[3] == 0) {
        /* hash(quaternion(w)) == hash(w), like complex */
        PyObject *f = PyFloat_FromDouble(q.data[0]);
        if (f == NULL) {
            return -1;
        }
        h = PyObject_Hash(f);
        Py_DECREF(f);
        return h;
    }
    tuple = Py_BuildValue("(dddd)", q.data[0], q.data[1], q.data[2], q.data[3]);
    if (tuple == NULL) {
        return -1;
    }
    h = PyObject_Hash(tuple);
    Py_DECREF(tuple);
    return h;
}

static PyObject *
pyquaternion_richcompare(PyObject *a, PyObject *b, int op)
{
    quaternion x, y;
    int r;
    if ((r = quaternion_from_object(a, &x)) <= 0 ||
            (r = quaternion_from_object(b, &y)) <= 0) {
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

#define QUATERNION_BINARY(name, expr)                                      \
    static PyObject *                                                      \
    pyquaternion_##name(PyObject *a, PyObject *b)                          \
    {                                                                      \
        quaternion x, y;                                                   \
        int r;                                                             \
        if ((r = quaternion_from_object(a, &x)) <= 0 ||                    \
                (r = quaternion_from_object(b, &y)) <= 0) {                \
            if (r < 0) {                                                   \
                return NULL;                                               \
            }                                                              \
            Py_RETURN_NOTIMPLEMENTED;                                      \
        }                                                                  \
        return NEW_QUATERNION(expr(x, y));                                 \
    }

QUATERNION_BINARY(add, vec4_add)
QUATERNION_BINARY(sub, vec4_sub)
QUATERNION_BINARY(mul, quaternion_mul)
QUATERNION_BINARY(div, quaternion_div)

#define QUATERNION_UNARY(name, expr)                                       \
    static PyObject *                                                      \
    pyquaternion_##name(PyObject *self)                                    \
    {                                                                      \
        return NEW_QUATERNION(expr(((PyQuaternion *)self)->q));            \
    }

QUATERNION_UNARY(pos, vec4_pos)
QUATERNION_UNARY(neg, vec4_neg)

static PyObject *
pyquaternion_zero(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyQuaternion_FromQuaternion((PyTypeObject *)cls, vec4_zero());
}

static PyObject *
pyquaternion_one(PyObject *cls, PyObject *Py_UNUSED(args))
{
    return PyQuaternion_FromQuaternion((PyTypeObject *)cls, quaternion_one());
}

static int
pyquaternion_bool(PyObject *self)
{
    return vec4_any_nonzero(((PyQuaternion *)self)->q);
}

static PyObject *
pyquaternion_exp(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_QUATERNION(quaternion_exp(((PyQuaternion *)self)->q));
}

static PyObject *
pyquaternion_log(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_QUATERNION(quaternion_log(((PyQuaternion *)self)->q));
}

static PyObject *
pyquaternion_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(quaternion_norm2(((PyQuaternion *)self)->q));
}

static PyObject *
pyquaternion_conjugate(PyObject *self, PyObject *Py_UNUSED(args))
{
    return NEW_QUATERNION(quaternion_conjugate(((PyQuaternion *)self)->q));
}

static PyObject *
pyquaternion_minpoly(PyObject *self, PyObject *Py_UNUSED(args))
{
    quaternion q = ((PyQuaternion *)self)->q;
    npy_intp dims[1] = {3};
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double *d;
    if (arr == NULL) {
        return NULL;
    }
    d = (double *)PyArray_DATA((PyArrayObject *)arr);
    quaternion_minpoly(q, d);
    return arr;
}

static PyObject *
pyquaternion_eivals(PyObject *self, PyObject *Py_UNUSED(args))
{
    quaternion q = ((PyQuaternion *)self)->q;
    npy_intp dims[1] = {2};
    vec2 ev[2];
    npy_cdouble *d;
    PyObject *arr = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    int i;
    if (arr == NULL) {
        return NULL;
    }
    quaternion_eivals(q, ev);
    d = (npy_cdouble *)PyArray_DATA((PyArrayObject *)arr);
    for (i = 0; i < 2; i++) {
        npy_csetreal(&d[i], ev[i].data[0]);
        npy_csetimag(&d[i], ev[i].data[1]);
    }
    return arr;
}

static PyObject *
pyquaternion_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    quaternion q = ((PyQuaternion *)self)->q;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         q.data[0], q.data[1], q.data[2], q.data[3]);
}

static PyObject *
pyquaternion_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {ALGPY_QUATERNION_N};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_quaternion);
    if (arr == NULL) {
        return NULL;
    }
    quaternion_basis((quaternion *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyquaternion_coeffs(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[2] = {ALGPY_QUATERNION_N, ALGPY_QUATERNION_N};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        2, dims, npy_quaternion);
    if (arr == NULL) {
        return NULL;
    }
    quaternion_coeffs((quaternion *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyquaternion_norm(PyObject *self, PyObject *Py_UNUSED(args))
{
    return algpy_norm_to_py(complex_root_real(quaternion_norm2(((PyQuaternion *)self)->q), 2));
}

static PyMethodDef pyquaternion_methods[] = {
    {"zero", pyquaternion_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"one", pyquaternion_one, METH_NOARGS | METH_CLASS, "Multiplicative identity."},
    {"basis", pyquaternion_basis, METH_NOARGS | METH_CLASS, "Algebra basis."},
    {"coeffs", pyquaternion_coeffs, METH_NOARGS | METH_CLASS, "Basis multiplication table."},
    {"norm2", pyquaternion_norm2, METH_NOARGS, "Squared norm."},
    {"norm", pyquaternion_norm, METH_NOARGS, "Norm (a float, or a complex1 if imaginary)."},
    {"conjugate", pyquaternion_conjugate, METH_NOARGS, "Quaternion conjugate."},
    {"minpoly", pyquaternion_minpoly, METH_NOARGS, "Coefficients of the minimal polynomial."},
    {"eivals", pyquaternion_eivals, METH_NOARGS, "Eigenvalues as an array of Python complex."},
    {"exp", pyquaternion_exp, METH_NOARGS, "Quaternion exponential."},
    {"log", pyquaternion_log, METH_NOARGS, "Principal quaternion logarithm."},
    {"__reduce__", pyquaternion_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyquaternion_length(PyObject *self)
{
    return 4;
}

static PyObject *
pyquaternion_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyQuaternion *)self)->q.data[i]);
}

/* numpy's generic scalar defines its own subscript, so override it. */
static PyObject *
pyquaternion_subscript(PyObject *self, PyObject *key)
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
    return pyquaternion_item(self, i);
}

static int
pyquaternion_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyquaternion_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

/* `tp_base` is only known after `import_array()`, so it is passed in init. */
static PyType_Slot pyquaternion_slots[] = {
    {Py_tp_dealloc, pyquaternion_dealloc},
    {Py_tp_traverse, pyquaternion_traverse},
    /* Must be set: `np.generic`'s allocator is not GC-aware. */
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyquaternion_repr},
    {Py_tp_str, pyquaternion_repr},
    {Py_tp_hash, pyquaternion_hash},
    {Py_tp_doc, (void *)"Quaternion with float64 components"},
    {Py_tp_richcompare, pyquaternion_richcompare},
    {Py_tp_methods, pyquaternion_methods},
    {Py_sq_length, pyquaternion_length},
    {Py_sq_item, pyquaternion_item},
    {Py_mp_subscript, pyquaternion_subscript},
    {Py_tp_new, pyquaternion_new},
    {Py_nb_add, pyquaternion_add},
    {Py_nb_subtract, pyquaternion_sub},
    {Py_nb_multiply, pyquaternion_mul},
    {Py_nb_true_divide, pyquaternion_div},
    {Py_nb_positive, pyquaternion_pos},
    {Py_nb_negative, pyquaternion_neg},
    {Py_nb_bool, pyquaternion_bool},
    {0, NULL}
};

static PyType_Spec pyquaternion_spec = {
    .name = "algpy.algebras.cayleydickson.quaternion",
    .basicsize = sizeof(PyQuaternion),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyquaternion_slots,
};

/* NumPy dtype support */

static PyObject *
npyquaternion_getitem(void *data, void *arr)
{
    quaternion q;
    memcpy(&q, data, sizeof(q));
    return NEW_QUATERNION(q);
}

static int
npyquaternion_setitem(PyObject *item, void *data, void *arr)
{
    quaternion q;
    int r = quaternion_from_object(item, &q);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to quaternion");
        double v[4];
        int i;
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "quaternion needs exactly 4 components");
            return -1;
        }
        for (i = 0; i < 4; i++) {
            v[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (v[i] == -1.0 && PyErr_Occurred()) {
                Py_DECREF(seq);
                return -1;
            }
        }
        Py_DECREF(seq);
        q = quaternion_make(v[0], v[1], v[2], v[3]);
    }
    memcpy(data, &q, sizeof(q));
    return 0;
}

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static void
npyquaternion_copyswapn(void *dst, npy_intp dstride, void *src,
                        npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(quaternion));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
npyquaternion_copyswap(void *dst, void *src, int swap, void *arr)
{
    if (src != NULL) {
        memcpy(dst, src, sizeof(quaternion));
    }
    if (swap) {
        swap_bytes((char *)dst);
    }
}

static npy_bool
npyquaternion_nonzero(void *data, void *arr)
{
    quaternion q;
    memcpy(&q, data, sizeof(q));
    return vec4_any_nonzero(q) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs npyquaternion_arrfuncs;

typedef struct { char c; quaternion q; } align_test;

PyArray_DescrProto npyquaternion_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL,                   /* typeobj, set to the `quaternion` type below */
    'V',                    /* kind */
    'q',                    /* type */
    '=',                    /* byteorder */
    NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM, /* flags */
    0,                      /* type_num */
    sizeof(quaternion),     /* elsize */
    offsetof(align_test, q),/* alignment */
    0,                      /* subarray */
    0,                      /* fields */
    0,                      /* names */
    &npyquaternion_arrfuncs,/* f */
};

/* Casts from the builtin real types */

#define QUATERNION_CAST_FROM(name, type)                                   \
    static void                                                            \
    npycast_##name##_to_quaternion(void *from, void *to, npy_intp n,       \
                                   void *fromarr, void *toarr)             \
    {                                                                      \
        const type *src = (const type *)from;                              \
        quaternion *dst = (quaternion *)to;                                \
        npy_intp i;                                                        \
        for (i = 0; i < n; i++) {                                          \
            dst[i] = quaternion_make((double)src[i], 0, 0, 0);             \
        }                                                                  \
    }

#define QUATERNION_REAL_TYPES(X)                                           \
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

QUATERNION_REAL_TYPES(QUATERNION_CAST_FROM)

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
                npycast_##name##_to_quaternion) < 0 ||                     \
                PyArray_RegisterCanCast(from, typenum, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                               \
            return -1;                                                     \
        }                                                                  \
        Py_DECREF(from);                                                   \
    }
    QUATERNION_REAL_TYPES(REGISTER_CAST)
#undef REGISTER_CAST
    return 0;
}

/* Ufunc loops */

typedef quaternion (*unary_op)(quaternion);
typedef quaternion (*binary_op)(quaternion, quaternion);
typedef int (*binary_pred)(quaternion, quaternion);

static void
loop_unary(char **args, npy_intp const *dims, npy_intp const *steps,
           void *data)
{
    unary_op op = (unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        quaternion a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_binary(char **args, npy_intp const *dims, npy_intp const *steps,
            void *data)
{
    binary_op op = (binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        quaternion a, b, r;
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
loop_predicate(char **args, npy_intp const *dims, npy_intp const *steps,
               void *data)
{
    binary_pred op = (binary_pred)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i, n = dims[0];
    for (i = 0; i < n; i++) {
        quaternion a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static int
quaternion_not_equal(quaternion a, quaternion b)
{
    return !vec4_all_eq(a, b);
}

static int
register_loop(PyObject *numpy, const char *name, PyUFuncGenericFunction loop,
              int nargs, int in_out_bool, void *data)
{
    int types[3];
    int i;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    int res;
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = npy_quaternion;
    }
    if (in_out_bool) {
        types[nargs - 1] = NPY_BOOL;
    }
    res = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_quaternion,
                                      loop, types, data);
    Py_DECREF(ufunc);
    return res;
}

static int
register_ufuncs(PyObject *numpy)
{
    return
        register_loop(numpy, "add", loop_binary, 3, 0, vec4_add) ||
        register_loop(numpy, "subtract", loop_binary, 3, 0,
                      vec4_sub) ||
        register_loop(numpy, "multiply", loop_binary, 3, 0,
                      quaternion_mul) ||
        register_loop(numpy, "divide", loop_binary, 3, 0, quaternion_div) ||
        register_loop(numpy, "positive", loop_unary, 2, 0,
                      vec4_pos) ||
        register_loop(numpy, "negative", loop_unary, 2, 0,
                      vec4_neg) ||
        register_loop(numpy, "exp", loop_unary, 2, 0, quaternion_exp) ||
        register_loop(numpy, "log", loop_unary, 2, 0, quaternion_log) ||
        register_loop(numpy, "equal", loop_predicate, 3, 1, vec4_all_eq) ||
        register_loop(numpy, "not_equal", loop_predicate, 3, 1,
                      quaternion_not_equal)
        ? -1 : 0;
}

/* Module */

static PyMethodDef module_methods[] = {
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_quaternion",
    "Quaternion number type for NumPy", -1, module_methods,
};

PyMODINIT_FUNC
PyInit__quaternion(void)
{
    PyObject *m = NULL;
    PyObject *numpy = NULL;
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

    PyQuaternion_Type = PyType_FromModuleAndSpec(
            m, &pyquaternion_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyQuaternion_Type == NULL) {
        goto fail;
    }

    PyArray_InitArrFuncs(&npyquaternion_arrfuncs);
    npyquaternion_arrfuncs.getitem = npyquaternion_getitem;
    npyquaternion_arrfuncs.setitem = npyquaternion_setitem;
    npyquaternion_arrfuncs.copyswapn = npyquaternion_copyswapn;
    npyquaternion_arrfuncs.copyswap = npyquaternion_copyswap;
    npyquaternion_arrfuncs.nonzero = npyquaternion_nonzero;

    npyquaternion_descr_proto.typeobj = (PyTypeObject *)PyQuaternion_Type;
    Py_SET_TYPE((PyObject *)&npyquaternion_descr_proto, &PyArrayDescr_Type);
    npy_quaternion = PyArray_RegisterDataType(&npyquaternion_descr_proto);
    if (npy_quaternion < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_quaternion);
    if (descr == NULL) {
        goto fail;
    }
    /* Support dtype(quaternion) syntax */
    if (PyObject_SetAttrString(PyQuaternion_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);

    if (register_casts(npy_quaternion) < 0 || register_ufuncs(numpy) < 0) {
        goto fail;
    }

    if (PyModule_AddObjectRef(m, "quaternion", PyQuaternion_Type) < 0) {
        goto fail;
    }
    Py_DECREF(numpy);
    return m;

fail:
    Py_XDECREF(numpy);
    Py_XDECREF(m);
    return NULL;
}
