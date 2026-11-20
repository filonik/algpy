/*
 * Exterior algebra on 1 generator exposed as a NumPy user-defined dtype with a
 * lattice interface: join (|), meet (&), ldual and rdual (also registered on
 * the shared algpy.lattices ufuncs) (and ~).
 * The arithmetic itself lives in exterior1.h.
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

#include "exterior1.h"

typedef struct {
    PyObject_HEAD
    exterior1 v;
} PyExterior1;

static PyObject *PyExterior1_Type = NULL;
static int npy_exterior1 = -1;

static PyObject *
PyExterior1_FromValue(PyTypeObject *type, exterior1 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyExterior1 *p = (PyExterior1 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_EXTERIOR1(v) \
    PyExterior1_FromValue((PyTypeObject *)PyExterior1_Type, (v))

static int
real_from_object(PyObject *o, double *out)
{
    if (PyFloat_Check(o) || PyLong_Check(o) ||
            PyArray_IsScalar(o, Integer) || PyArray_IsScalar(o, Floating) ||
            PyArray_IsScalar(o, Bool)) {
        double d = PyFloat_AsDouble(o);
        if (d == -1.0 && PyErr_Occurred()) {
            return -1;
        }
        *out = d;
        return 1;
    }
    return 0;
}

static int
exterior1_from_object(PyObject *o, exterior1 *out)
{
    double d;
    int r;
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyExterior1_Type)) {
        *out = ((PyExterior1 *)o)->v;
        return 1;
    }
    r = real_from_object(o, &d);
    if (r > 0) {
        *out = exterior1_from_re(d);
    }
    return r;
}

static PyObject *
pyexterior1_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[2] = {0};
    static char *kwlist[] = {"x", "x0", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dd", kwlist,
            &d[0], &d[1])) {
        return NULL;
    }
    return PyExterior1_FromValue(type, exterior1_make(d));
}

static PyObject *
pyexterior1_repr(PyObject *self)
{
    exterior1 v = ((PyExterior1 *)self)->v;
    char *s[2] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 2; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("exterior1(%s, %s)",
                                  s[0], s[1]);
done:
    for (i = 0; i < 2; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyexterior1_hash(PyObject *self)
{
    exterior1 v = ((PyExterior1 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_EXTERIOR1_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (exterior1_is_real(v)) {
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
pyexterior1_richcompare(PyObject *a, PyObject *b, int op)
{
    exterior1 x, y;
    int r;
    if ((r = exterior1_from_object(a, &x)) <= 0 ||
            (r = exterior1_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(exterior1_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!exterior1_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define EXTERIOR1_BINARY(name, expr)                                            \
    static PyObject *                                                        \
    pyexterior1_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        exterior1 x, y;                                                       \
        int r;                                                               \
        if ((r = exterior1_from_object(a, &x)) <= 0 ||                        \
                (r = exterior1_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_EXTERIOR1(expr(x, y));                                     \
    }

EXTERIOR1_BINARY(add, exterior1_add)
EXTERIOR1_BINARY(sub, exterior1_sub)
EXTERIOR1_BINARY(join, exterior1_join)
EXTERIOR1_BINARY(meet, exterior1_meet)

/* Multiplication and division only scale by real numbers. */
static PyObject *
pyexterior1_mul(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (PyObject_TypeCheck(a, (PyTypeObject *)PyExterior1_Type)) {
        if ((r = real_from_object(b, &d)) <= 0) {
            if (r < 0) {
                return NULL;
            }
            Py_RETURN_NOTIMPLEMENTED;
        }
        return NEW_EXTERIOR1(exterior1_lmul(((PyExterior1 *)a)->v, d));
    }
    if ((r = real_from_object(a, &d)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_EXTERIOR1(exterior1_rmul(d, ((PyExterior1 *)b)->v));
}

static PyObject *
pyexterior1_div(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (!PyObject_TypeCheck(a, (PyTypeObject *)PyExterior1_Type) ||
            (r = real_from_object(b, &d)) == 0) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (r < 0) {
        return NULL;
    }
    return NEW_EXTERIOR1(exterior1_ldiv(((PyExterior1 *)a)->v, d));
}

#define EXTERIOR1_UNARY(name, expr)                                             \
    static PyObject *                                                        \
    pyexterior1_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_EXTERIOR1(expr(((PyExterior1 *)self)->v));                  \
    }

EXTERIOR1_UNARY(pos, exterior1_pos)
EXTERIOR1_UNARY(neg, exterior1_neg)
EXTERIOR1_UNARY(invert, exterior1_dual)
static int
pyexterior1_bool(PyObject *self)
{
    exterior1 v = ((PyExterior1 *)self)->v;
    return vec2_any_nonzero(v);
}

#define EXTERIOR1_METHOD(name, op)                                              \
    static PyObject *                                                        \
    pyexterior1_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return PyExterior1_FromValue(Py_TYPE(self), op(((PyExterior1 *)self)->v));    \
    }

EXTERIOR1_METHOD(ldual, exterior1_ldual)
EXTERIOR1_METHOD(rdual, exterior1_rdual)
EXTERIOR1_METHOD(dual, exterior1_dual)
#define EXTERIOR1_CONSTANT(name, value)                                         \
    static PyObject *                                                        \
    pyexterior1_##name(PyObject *cls, PyObject *Py_UNUSED(args))              \
    {                                                                        \
        return PyExterior1_FromValue((PyTypeObject *)cls, value());                \
    }

EXTERIOR1_CONSTANT(zero, exterior1_zero)
EXTERIOR1_CONSTANT(bottom, exterior1_bottom)
EXTERIOR1_CONSTANT(top, exterior1_top)

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyexterior1_grade(PyObject *self, Py_ssize_t k)
{
    exterior1 v = ((PyExterior1 *)self)->v;
    switch (k) {
        case 0:
            return NEW_EXTERIOR1(exterior1_grade0(v));
        case 1:
            return NEW_EXTERIOR1(exterior1_grade1(v));
        default:
            return NEW_EXTERIOR1(exterior1_zero());
    }
}

static PyObject *
pyexterior1_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyexterior1_grade, 2);
}

static PyGetSetDef pyexterior1_getset[] = {
    {"grade", pyexterior1_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyexterior1_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(exterior1_norm2(((PyExterior1 *)self)->v));
}

static PyObject *
pyexterior1_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior1);
    if (arr == NULL) {
        return NULL;
    }
    exterior1_basis((exterior1 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior1_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior1);
    if (arr == NULL) {
        return NULL;
    }
    exterior1_gens((exterior1 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior1_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    exterior1 v = ((PyExterior1 *)self)->v;
    return Py_BuildValue("(O(dd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1]);
}

static PyObject *
pyexterior1_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_exterior1, 2);
}

static PyObject *
pyexterior1_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_exterior1, 2);
}

static PyMethodDef pyexterior1_methods[] = {
    {"to_float", (PyCFunction)pyexterior1_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert lattice values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyexterior1_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the lattice dtype."},
    {"zero", pyexterior1_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"bottom", pyexterior1_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pyexterior1_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"basis", pyexterior1_basis, METH_NOARGS | METH_CLASS, "Coefficient basis."},
    {"gens", pyexterior1_gens, METH_NOARGS | METH_CLASS, "Generators in basis order; excludes bottom."},
    {"ldual", pyexterior1_ldual, METH_NOARGS, "Left dual: a | a.ldual() == top for basis blades."},
    {"rdual", pyexterior1_rdual, METH_NOARGS, "Right dual: a.rdual() | a == top for basis blades."},
    {"dual", pyexterior1_dual, METH_NOARGS, "Dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"norm2", pyexterior1_norm2, METH_NOARGS, "Sum of squared coefficients."},
    {"__reduce__", pyexterior1_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyexterior1_length(PyObject *Py_UNUSED(self))
{
    return 2;
}

static PyObject *
pyexterior1_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 2) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyExterior1 *)self)->v.data[i]);
}

static PyObject *
pyexterior1_subscript(PyObject *self, PyObject *key)
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
    return pyexterior1_item(self, i);
}

static int
pyexterior1_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyexterior1_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyexterior1_slots[] = {
    {Py_tp_dealloc, pyexterior1_dealloc},
    {Py_tp_traverse, pyexterior1_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyexterior1_repr},
    {Py_tp_str, pyexterior1_repr},
    {Py_tp_hash, pyexterior1_hash},
    {Py_tp_doc, (void *)"Exterior algebra on 1 generator as a lattice with 2 float64 coefficients"},
    {Py_tp_richcompare, pyexterior1_richcompare},
    {Py_tp_methods, pyexterior1_methods},
    {Py_tp_getset, pyexterior1_getset},
    {Py_sq_length, pyexterior1_length},
    {Py_sq_item, pyexterior1_item},
    {Py_mp_subscript, pyexterior1_subscript},
    {Py_tp_new, pyexterior1_new},
    {Py_nb_add, pyexterior1_add},
    {Py_nb_subtract, pyexterior1_sub},
    {Py_nb_multiply, pyexterior1_mul},
    {Py_nb_true_divide, pyexterior1_div},
    {Py_nb_or, pyexterior1_join},
    {Py_nb_and, pyexterior1_meet},
    {Py_nb_positive, pyexterior1_pos},
    {Py_nb_negative, pyexterior1_neg},
    {Py_nb_invert, pyexterior1_invert},
    {Py_nb_bool, pyexterior1_bool},
    {0, NULL}
};

static PyType_Spec pyexterior1_spec = {
    .name = "algpy.lattices.clifford.exterior1",
    .basicsize = sizeof(PyExterior1),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyexterior1_slots,
};

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

static PyObject *
npyexterior1_getitem(void *data, void *arr)
{
    exterior1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_EXTERIOR1(v);
}

static int
npyexterior1_setitem(PyObject *item, void *data, void *arr)
{
    exterior1 v;
    int r = exterior1_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to exterior1");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 2) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "exterior1 needs exactly 2 components");
            return -1;
        }
        for (int i = 0; i < 2; i++) {
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
pyexterior1_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(exterior1));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyexterior1_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyexterior1_copyswapn(dst, sizeof(exterior1), src, sizeof(exterior1), 1,
                         swap, arr);
}

static npy_bool
pyexterior1_nonzero(void *data, void *arr)
{
    exterior1 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec2_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyexterior1_arrfuncs;
typedef struct { char c; exterior1 v; } pyexterior1_align_test;

static PyArray_DescrProto pyexterior1_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(exterior1), offsetof(pyexterior1_align_test, v), 0, 0, 0,
    &pyexterior1_arrfuncs,
};

#define EXTERIOR1_CAST_FROM(name, type, convert)                                \
    static void                                                              \
    pyexterior1_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        exterior1 *dst = (exterior1 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = exterior1_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define EXTERIOR1_REAL_TYPES(X) \
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

EXTERIOR1_REAL_TYPES(EXTERIOR1_CAST_FROM)

static int
pyexterior1_register_casts(void)
{
#define REGISTER_EXTERIOR1_CAST(name, type, convert)                            \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_exterior1,                    \
                                     pyexterior1_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_exterior1, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    EXTERIOR1_REAL_TYPES(REGISTER_EXTERIOR1_CAST)
#undef REGISTER_EXTERIOR1_CAST
    return 0;
}

typedef exterior1 (*exterior1_unary_op)(exterior1);
typedef exterior1 (*exterior1_binary_op)(exterior1, exterior1);
typedef int (*exterior1_predicate)(exterior1, exterior1);

static void
loop_exterior1_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    exterior1_unary_op op = (exterior1_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior1 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_exterior1_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior1_binary_op op = (exterior1_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior1 a, b, r;
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
loop_exterior1_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    exterior1_predicate op = (exterior1_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior1 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_exterior1_rmul(char **args, npy_intp const *dims,
                   npy_intp const *steps, void *Py_UNUSED(data))
{
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        double a;
        exterior1 b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = exterior1_rmul(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

typedef exterior1 (*exterior1_scalar_op)(exterior1, double);

static void
loop_exterior1_scalar(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior1_scalar_op op = (exterior1_scalar_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior1 a, r;
        double b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = op(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

/* Register a loop for name; a negative entry in types means the lattice dtype. */
static int
pyexterior1_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs,
                         const int *signature, void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = signature[i] < 0 ? npy_exterior1 : signature[i];
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_exterior1,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyexterior1_register_ufuncs(PyObject *numpy)
{
    static const int unary[2] = {-1, -1};
    static const int binary[3] = {-1, -1, -1};
    static const int predicate[3] = {-1, -1, NPY_BOOL};
    static const int rscalar[3] = {NPY_DOUBLE, -1, -1};
    static const int lscalar[3] = {-1, NPY_DOUBLE, -1};
    return
        pyexterior1_register_loop(numpy, "add", loop_exterior1_binary, 3, binary, exterior1_add) ||
        pyexterior1_register_loop(numpy, "subtract", loop_exterior1_binary, 3, binary, exterior1_sub) ||
        pyexterior1_register_loop(numpy, "multiply", loop_exterior1_rmul, 3, rscalar, NULL) ||
        pyexterior1_register_loop(numpy, "multiply", loop_exterior1_scalar, 3, lscalar, exterior1_lmul) ||
        pyexterior1_register_loop(numpy, "divide", loop_exterior1_scalar, 3, lscalar, exterior1_ldiv) ||
        pyexterior1_register_loop(numpy, "bitwise_or", loop_exterior1_binary, 3, binary, exterior1_join) ||
        pyexterior1_register_loop(numpy, "bitwise_and", loop_exterior1_binary, 3, binary, exterior1_meet) ||
        pyexterior1_register_loop(numpy, "invert", loop_exterior1_unary, 2, unary, exterior1_dual) ||
        pyexterior1_register_loop(numpy, "positive", loop_exterior1_unary, 2, unary, exterior1_pos) ||
        pyexterior1_register_loop(numpy, "negative", loop_exterior1_unary, 2, unary, exterior1_neg) ||
        pyexterior1_register_loop(numpy, "equal", loop_exterior1_predicate, 3, predicate, exterior1_equal) ||
        pyexterior1_register_loop(numpy, "not_equal", loop_exterior1_predicate, 3, predicate, exterior1_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pyexterior1_register_lattice_ufuncs(PyObject *ufuncs)
{
    static const int unary[2] = {-1, -1};
    return
        pyexterior1_register_loop(ufuncs, "ldual", loop_exterior1_unary, 2, unary, exterior1_ldual) ||
        pyexterior1_register_loop(ufuncs, "rdual", loop_exterior1_unary, 2, unary, exterior1_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_exterior1", "exterior1 lattice type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__exterior1(void)
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
    PyExterior1_Type = PyType_FromModuleAndSpec(
        m, &pyexterior1_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyExterior1_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyexterior1_arrfuncs);
    pyexterior1_arrfuncs.getitem = npyexterior1_getitem;
    pyexterior1_arrfuncs.setitem = npyexterior1_setitem;
    pyexterior1_arrfuncs.copyswapn = pyexterior1_copyswapn;
    pyexterior1_arrfuncs.copyswap = pyexterior1_copyswap;
    pyexterior1_arrfuncs.nonzero = pyexterior1_nonzero;
    pyexterior1_descr_proto.typeobj = (PyTypeObject *)PyExterior1_Type;
    Py_SET_TYPE((PyObject *)&pyexterior1_descr_proto, &PyArrayDescr_Type);
    npy_exterior1 = PyArray_RegisterDataType(&pyexterior1_descr_proto);
    if (npy_exterior1 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_exterior1);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyExterior1_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyExterior1_Type, 2) < 0) {
        goto fail;
    }
    if (pyexterior1_register_casts() < 0 ||
            pyexterior1_register_ufuncs(numpy) < 0 ||
            pyexterior1_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "exterior1", PyExterior1_Type) < 0) {
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
