/*
 * Exterior algebra on 2 generators exposed as a NumPy user-defined dtype with a
 * lattice interface: join (|), meet (&), ldual and rdual (also registered on
 * the shared algpy.lattices ufuncs).
 * The arithmetic itself lives in exterior2.h.
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

#include "exterior2.h"

typedef struct {
    PyObject_HEAD
    exterior2 v;
} PyExterior2;

static PyObject *PyExterior2_Type = NULL;
static int npy_exterior2 = -1;

static PyObject *
PyExterior2_FromValue(PyTypeObject *type, exterior2 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyExterior2 *p = (PyExterior2 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_EXTERIOR2(v) \
    PyExterior2_FromValue((PyTypeObject *)PyExterior2_Type, (v))

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
exterior2_from_object(PyObject *o, exterior2 *out)
{
    double d;
    int r;
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyExterior2_Type)) {
        *out = ((PyExterior2 *)o)->v;
        return 1;
    }
    r = real_from_object(o, &d);
    if (r > 0) {
        *out = exterior2_from_re(d);
    }
    return r;
}

static PyObject *
pyexterior2_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[4] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
            &d[0], &d[1], &d[2], &d[3])) {
        return NULL;
    }
    return PyExterior2_FromValue(type, exterior2_make(d));
}

static PyObject *
pyexterior2_repr(PyObject *self)
{
    exterior2 v = ((PyExterior2 *)self)->v;
    char *s[4] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 4; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("exterior2(%s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3]);
done:
    for (i = 0; i < 4; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyexterior2_hash(PyObject *self)
{
    exterior2 v = ((PyExterior2 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_EXTERIOR2_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (exterior2_is_real(v)) {
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
pyexterior2_richcompare(PyObject *a, PyObject *b, int op)
{
    exterior2 x, y;
    int r;
    if ((r = exterior2_from_object(a, &x)) <= 0 ||
            (r = exterior2_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(exterior2_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!exterior2_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define EXTERIOR2_BINARY(name, expr)                                            \
    static PyObject *                                                        \
    pyexterior2_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        exterior2 x, y;                                                       \
        int r;                                                               \
        if ((r = exterior2_from_object(a, &x)) <= 0 ||                        \
                (r = exterior2_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_EXTERIOR2(expr(x, y));                                     \
    }

EXTERIOR2_BINARY(add, exterior2_add)
EXTERIOR2_BINARY(sub, exterior2_sub)
EXTERIOR2_BINARY(join, exterior2_join)
EXTERIOR2_BINARY(meet, exterior2_meet)

/* Multiplication and division only scale by real numbers. */
static PyObject *
pyexterior2_mul(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (PyObject_TypeCheck(a, (PyTypeObject *)PyExterior2_Type)) {
        if ((r = real_from_object(b, &d)) <= 0) {
            if (r < 0) {
                return NULL;
            }
            Py_RETURN_NOTIMPLEMENTED;
        }
        return NEW_EXTERIOR2(exterior2_lmul(((PyExterior2 *)a)->v, d));
    }
    if ((r = real_from_object(a, &d)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_EXTERIOR2(exterior2_rmul(d, ((PyExterior2 *)b)->v));
}

static PyObject *
pyexterior2_div(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (!PyObject_TypeCheck(a, (PyTypeObject *)PyExterior2_Type) ||
            (r = real_from_object(b, &d)) == 0) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (r < 0) {
        return NULL;
    }
    return NEW_EXTERIOR2(exterior2_ldiv(((PyExterior2 *)a)->v, d));
}

#define EXTERIOR2_UNARY(name, expr)                                             \
    static PyObject *                                                        \
    pyexterior2_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_EXTERIOR2(expr(((PyExterior2 *)self)->v));                  \
    }

EXTERIOR2_UNARY(pos, exterior2_pos)
EXTERIOR2_UNARY(neg, exterior2_neg)

/* ldual and rdual differ in even dimensions, so ~ is deliberately undefined. */
static PyObject *
pyexterior2_invert(PyObject *Py_UNUSED(self))
{
    PyErr_SetString(PyExc_TypeError,
                    "~ is undefined for exterior2: ldual and rdual differ in "
                    "even dimensions; use ldual() or rdual()");
    return NULL;
}

static int
pyexterior2_bool(PyObject *self)
{
    exterior2 v = ((PyExterior2 *)self)->v;
    return vec4_any_nonzero(v);
}

#define EXTERIOR2_METHOD(name, op)                                              \
    static PyObject *                                                        \
    pyexterior2_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return PyExterior2_FromValue(Py_TYPE(self), op(((PyExterior2 *)self)->v));    \
    }

EXTERIOR2_METHOD(ldual, exterior2_ldual)
EXTERIOR2_METHOD(rdual, exterior2_rdual)

#define EXTERIOR2_CONSTANT(name, value)                                         \
    static PyObject *                                                        \
    pyexterior2_##name(PyObject *cls, PyObject *Py_UNUSED(args))              \
    {                                                                        \
        return PyExterior2_FromValue((PyTypeObject *)cls, value());                \
    }

EXTERIOR2_CONSTANT(zero, exterior2_zero)
EXTERIOR2_CONSTANT(bottom, exterior2_bottom)
EXTERIOR2_CONSTANT(top, exterior2_top)

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyexterior2_grade(PyObject *self, Py_ssize_t k)
{
    exterior2 v = ((PyExterior2 *)self)->v;
    switch (k) {
        case 0:
            return NEW_EXTERIOR2(exterior2_grade0(v));
        case 1:
            return NEW_EXTERIOR2(exterior2_grade1(v));
        case 2:
            return NEW_EXTERIOR2(exterior2_grade2(v));
        default:
            return NEW_EXTERIOR2(exterior2_zero());
    }
}

static PyObject *
pyexterior2_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyexterior2_grade, 3);
}

static PyGetSetDef pyexterior2_getset[] = {
    {"grade", pyexterior2_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyexterior2_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(exterior2_norm2(((PyExterior2 *)self)->v));
}

static PyObject *
pyexterior2_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior2);
    if (arr == NULL) {
        return NULL;
    }
    exterior2_basis((exterior2 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior2_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {2};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior2);
    if (arr == NULL) {
        return NULL;
    }
    exterior2_gens((exterior2 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior2_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    exterior2 v = ((PyExterior2 *)self)->v;
    return Py_BuildValue("(O(dddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3]);
}

static PyObject *
pyexterior2_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_exterior2, 4);
}

static PyObject *
pyexterior2_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_exterior2, 4);
}

static PyMethodDef pyexterior2_methods[] = {
    {"to_float", (PyCFunction)pyexterior2_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert lattice values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyexterior2_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the lattice dtype."},
    {"zero", pyexterior2_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"bottom", pyexterior2_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pyexterior2_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"basis", pyexterior2_basis, METH_NOARGS | METH_CLASS, "Coefficient basis."},
    {"gens", pyexterior2_gens, METH_NOARGS | METH_CLASS, "Generators in basis order; excludes bottom."},
    {"ldual", pyexterior2_ldual, METH_NOARGS, "Left dual: a | a.ldual() == top for basis blades."},
    {"rdual", pyexterior2_rdual, METH_NOARGS, "Right dual: a.rdual() | a == top for basis blades."},
    {"norm2", pyexterior2_norm2, METH_NOARGS, "Sum of squared coefficients."},
    {"__reduce__", pyexterior2_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyexterior2_length(PyObject *Py_UNUSED(self))
{
    return 4;
}

static PyObject *
pyexterior2_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 4) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyExterior2 *)self)->v.data[i]);
}

static PyObject *
pyexterior2_subscript(PyObject *self, PyObject *key)
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
    return pyexterior2_item(self, i);
}

static int
pyexterior2_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyexterior2_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyexterior2_slots[] = {
    {Py_tp_dealloc, pyexterior2_dealloc},
    {Py_tp_traverse, pyexterior2_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyexterior2_repr},
    {Py_tp_str, pyexterior2_repr},
    {Py_tp_hash, pyexterior2_hash},
    {Py_tp_doc, (void *)"Exterior algebra on 2 generators as a lattice with 4 float64 coefficients"},
    {Py_tp_richcompare, pyexterior2_richcompare},
    {Py_tp_methods, pyexterior2_methods},
    {Py_tp_getset, pyexterior2_getset},
    {Py_sq_length, pyexterior2_length},
    {Py_sq_item, pyexterior2_item},
    {Py_mp_subscript, pyexterior2_subscript},
    {Py_tp_new, pyexterior2_new},
    {Py_nb_add, pyexterior2_add},
    {Py_nb_subtract, pyexterior2_sub},
    {Py_nb_multiply, pyexterior2_mul},
    {Py_nb_true_divide, pyexterior2_div},
    {Py_nb_or, pyexterior2_join},
    {Py_nb_and, pyexterior2_meet},
    {Py_nb_positive, pyexterior2_pos},
    {Py_nb_negative, pyexterior2_neg},
    {Py_nb_invert, pyexterior2_invert},
    {Py_nb_bool, pyexterior2_bool},
    {0, NULL}
};

static PyType_Spec pyexterior2_spec = {
    .name = "algpy.lattices.clifford.exterior2",
    .basicsize = sizeof(PyExterior2),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyexterior2_slots,
};

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

static PyObject *
npyexterior2_getitem(void *data, void *arr)
{
    exterior2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_EXTERIOR2(v);
}

static int
npyexterior2_setitem(PyObject *item, void *data, void *arr)
{
    exterior2 v;
    int r = exterior2_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to exterior2");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 4) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "exterior2 needs exactly 4 components");
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
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    memcpy(data, &v, sizeof(v));
    return 0;
}

static void
pyexterior2_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(exterior2));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyexterior2_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyexterior2_copyswapn(dst, sizeof(exterior2), src, sizeof(exterior2), 1,
                         swap, arr);
}

static npy_bool
pyexterior2_nonzero(void *data, void *arr)
{
    exterior2 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec4_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyexterior2_arrfuncs;
typedef struct { char c; exterior2 v; } pyexterior2_align_test;

static PyArray_DescrProto pyexterior2_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(exterior2), offsetof(pyexterior2_align_test, v), 0, 0, 0,
    &pyexterior2_arrfuncs,
};

#define EXTERIOR2_CAST_FROM(name, type, convert)                                \
    static void                                                              \
    pyexterior2_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        exterior2 *dst = (exterior2 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = exterior2_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define EXTERIOR2_REAL_TYPES(X) \
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

EXTERIOR2_REAL_TYPES(EXTERIOR2_CAST_FROM)

static int
pyexterior2_register_casts(void)
{
#define REGISTER_EXTERIOR2_CAST(name, type, convert)                            \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_exterior2,                    \
                                     pyexterior2_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_exterior2, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    EXTERIOR2_REAL_TYPES(REGISTER_EXTERIOR2_CAST)
#undef REGISTER_EXTERIOR2_CAST
    return 0;
}

typedef exterior2 (*exterior2_unary_op)(exterior2);
typedef exterior2 (*exterior2_binary_op)(exterior2, exterior2);
typedef int (*exterior2_predicate)(exterior2, exterior2);

static void
loop_exterior2_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    exterior2_unary_op op = (exterior2_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior2 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_exterior2_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior2_binary_op op = (exterior2_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior2 a, b, r;
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
loop_exterior2_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    exterior2_predicate op = (exterior2_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior2 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_exterior2_rmul(char **args, npy_intp const *dims,
                   npy_intp const *steps, void *Py_UNUSED(data))
{
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        double a;
        exterior2 b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = exterior2_rmul(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

typedef exterior2 (*exterior2_scalar_op)(exterior2, double);

static void
loop_exterior2_scalar(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior2_scalar_op op = (exterior2_scalar_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior2 a, r;
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
pyexterior2_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs,
                         const int *signature, void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = signature[i] < 0 ? npy_exterior2 : signature[i];
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_exterior2,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyexterior2_register_ufuncs(PyObject *numpy)
{
    static const int unary[2] = {-1, -1};
    static const int binary[3] = {-1, -1, -1};
    static const int predicate[3] = {-1, -1, NPY_BOOL};
    static const int rscalar[3] = {NPY_DOUBLE, -1, -1};
    static const int lscalar[3] = {-1, NPY_DOUBLE, -1};
    return
        pyexterior2_register_loop(numpy, "add", loop_exterior2_binary, 3, binary, exterior2_add) ||
        pyexterior2_register_loop(numpy, "subtract", loop_exterior2_binary, 3, binary, exterior2_sub) ||
        pyexterior2_register_loop(numpy, "multiply", loop_exterior2_rmul, 3, rscalar, NULL) ||
        pyexterior2_register_loop(numpy, "multiply", loop_exterior2_scalar, 3, lscalar, exterior2_lmul) ||
        pyexterior2_register_loop(numpy, "divide", loop_exterior2_scalar, 3, lscalar, exterior2_ldiv) ||
        pyexterior2_register_loop(numpy, "bitwise_or", loop_exterior2_binary, 3, binary, exterior2_join) ||
        pyexterior2_register_loop(numpy, "bitwise_and", loop_exterior2_binary, 3, binary, exterior2_meet) ||
        pyexterior2_register_loop(numpy, "positive", loop_exterior2_unary, 2, unary, exterior2_pos) ||
        pyexterior2_register_loop(numpy, "negative", loop_exterior2_unary, 2, unary, exterior2_neg) ||
        pyexterior2_register_loop(numpy, "equal", loop_exterior2_predicate, 3, predicate, exterior2_equal) ||
        pyexterior2_register_loop(numpy, "not_equal", loop_exterior2_predicate, 3, predicate, exterior2_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pyexterior2_register_lattice_ufuncs(PyObject *ufuncs)
{
    static const int unary[2] = {-1, -1};
    return
        pyexterior2_register_loop(ufuncs, "ldual", loop_exterior2_unary, 2, unary, exterior2_ldual) ||
        pyexterior2_register_loop(ufuncs, "rdual", loop_exterior2_unary, 2, unary, exterior2_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_exterior2", "exterior2 lattice type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__exterior2(void)
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
    PyExterior2_Type = PyType_FromModuleAndSpec(
        m, &pyexterior2_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyExterior2_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyexterior2_arrfuncs);
    pyexterior2_arrfuncs.getitem = npyexterior2_getitem;
    pyexterior2_arrfuncs.setitem = npyexterior2_setitem;
    pyexterior2_arrfuncs.copyswapn = pyexterior2_copyswapn;
    pyexterior2_arrfuncs.copyswap = pyexterior2_copyswap;
    pyexterior2_arrfuncs.nonzero = pyexterior2_nonzero;
    pyexterior2_descr_proto.typeobj = (PyTypeObject *)PyExterior2_Type;
    Py_SET_TYPE((PyObject *)&pyexterior2_descr_proto, &PyArrayDescr_Type);
    npy_exterior2 = PyArray_RegisterDataType(&pyexterior2_descr_proto);
    if (npy_exterior2 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_exterior2);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyExterior2_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyExterior2_Type, 4) < 0) {
        goto fail;
    }
    if (pyexterior2_register_casts() < 0 ||
            pyexterior2_register_ufuncs(numpy) < 0 ||
            pyexterior2_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "exterior2", PyExterior2_Type) < 0) {
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
