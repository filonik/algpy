/*
 * Exterior algebra on 4 generators exposed as a NumPy user-defined dtype with a
 * lattice interface: join (|), meet (&), ldual and rdual (also registered on
 * the shared algpy.lattices ufuncs).
 * The arithmetic itself lives in exterior4.h.
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

#include "exterior4.h"

typedef struct {
    PyObject_HEAD
    exterior4 v;
} PyExterior4;

static PyObject *PyExterior4_Type = NULL;
static int npy_exterior4 = -1;

static PyObject *
PyExterior4_FromValue(PyTypeObject *type, exterior4 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyExterior4 *p = (PyExterior4 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_EXTERIOR4(v) \
    PyExterior4_FromValue((PyTypeObject *)PyExterior4_Type, (v))

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
exterior4_from_object(PyObject *o, exterior4 *out)
{
    double d;
    int r;
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyExterior4_Type)) {
        *out = ((PyExterior4 *)o)->v;
        return 1;
    }
    r = real_from_object(o, &d);
    if (r > 0) {
        *out = exterior4_from_re(d);
    }
    return r;
}

static PyObject *
pyexterior4_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[16] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15])) {
        return NULL;
    }
    return PyExterior4_FromValue(type, exterior4_make(d));
}

static PyObject *
pyexterior4_repr(PyObject *self)
{
    exterior4 v = ((PyExterior4 *)self)->v;
    char *s[16] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 16; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("exterior4(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15]);
done:
    for (i = 0; i < 16; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyexterior4_hash(PyObject *self)
{
    exterior4 v = ((PyExterior4 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_EXTERIOR4_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (exterior4_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyexterior4_richcompare(PyObject *a, PyObject *b, int op)
{
    exterior4 x, y;
    int r;
    if ((r = exterior4_from_object(a, &x)) <= 0 ||
            (r = exterior4_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(exterior4_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!exterior4_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define EXTERIOR4_BINARY(name, expr)                                            \
    static PyObject *                                                        \
    pyexterior4_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        exterior4 x, y;                                                       \
        int r;                                                               \
        if ((r = exterior4_from_object(a, &x)) <= 0 ||                        \
                (r = exterior4_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_EXTERIOR4(expr(x, y));                                     \
    }

EXTERIOR4_BINARY(add, exterior4_add)
EXTERIOR4_BINARY(sub, exterior4_sub)
EXTERIOR4_BINARY(join, exterior4_join)
EXTERIOR4_BINARY(meet, exterior4_meet)

/* Multiplication and division only scale by real numbers. */
static PyObject *
pyexterior4_mul(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (PyObject_TypeCheck(a, (PyTypeObject *)PyExterior4_Type)) {
        if ((r = real_from_object(b, &d)) <= 0) {
            if (r < 0) {
                return NULL;
            }
            Py_RETURN_NOTIMPLEMENTED;
        }
        return NEW_EXTERIOR4(exterior4_lmul(((PyExterior4 *)a)->v, d));
    }
    if ((r = real_from_object(a, &d)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_EXTERIOR4(exterior4_rmul(d, ((PyExterior4 *)b)->v));
}

static PyObject *
pyexterior4_div(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (!PyObject_TypeCheck(a, (PyTypeObject *)PyExterior4_Type) ||
            (r = real_from_object(b, &d)) == 0) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (r < 0) {
        return NULL;
    }
    return NEW_EXTERIOR4(exterior4_ldiv(((PyExterior4 *)a)->v, d));
}

#define EXTERIOR4_UNARY(name, expr)                                             \
    static PyObject *                                                        \
    pyexterior4_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_EXTERIOR4(expr(((PyExterior4 *)self)->v));                  \
    }

EXTERIOR4_UNARY(pos, exterior4_pos)
EXTERIOR4_UNARY(neg, exterior4_neg)

/* ldual and rdual differ in even dimensions, so ~ is deliberately undefined. */
static PyObject *
pyexterior4_invert(PyObject *Py_UNUSED(self))
{
    PyErr_SetString(PyExc_TypeError,
                    "~ is undefined for exterior4: ldual and rdual differ in "
                    "even dimensions; use ldual() or rdual()");
    return NULL;
}

static int
pyexterior4_bool(PyObject *self)
{
    exterior4 v = ((PyExterior4 *)self)->v;
    return vec16_any_nonzero(v);
}

#define EXTERIOR4_METHOD(name, op)                                              \
    static PyObject *                                                        \
    pyexterior4_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return PyExterior4_FromValue(Py_TYPE(self), op(((PyExterior4 *)self)->v));    \
    }

EXTERIOR4_METHOD(ldual, exterior4_ldual)
EXTERIOR4_METHOD(rdual, exterior4_rdual)

#define EXTERIOR4_CONSTANT(name, value)                                         \
    static PyObject *                                                        \
    pyexterior4_##name(PyObject *cls, PyObject *Py_UNUSED(args))              \
    {                                                                        \
        return PyExterior4_FromValue((PyTypeObject *)cls, value());                \
    }

EXTERIOR4_CONSTANT(zero, exterior4_zero)
EXTERIOR4_CONSTANT(bottom, exterior4_bottom)
EXTERIOR4_CONSTANT(top, exterior4_top)

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyexterior4_grade(PyObject *self, Py_ssize_t k)
{
    exterior4 v = ((PyExterior4 *)self)->v;
    switch (k) {
        case 0:
            return NEW_EXTERIOR4(exterior4_grade0(v));
        case 1:
            return NEW_EXTERIOR4(exterior4_grade1(v));
        case 2:
            return NEW_EXTERIOR4(exterior4_grade2(v));
        case 3:
            return NEW_EXTERIOR4(exterior4_grade3(v));
        case 4:
            return NEW_EXTERIOR4(exterior4_grade4(v));
        default:
            return NEW_EXTERIOR4(exterior4_zero());
    }
}

static PyObject *
pyexterior4_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyexterior4_grade, 5);
}

static PyGetSetDef pyexterior4_getset[] = {
    {"grade", pyexterior4_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyexterior4_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(exterior4_norm2(((PyExterior4 *)self)->v));
}

static PyObject *
pyexterior4_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {16};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior4);
    if (arr == NULL) {
        return NULL;
    }
    exterior4_basis((exterior4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior4_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {4};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior4);
    if (arr == NULL) {
        return NULL;
    }
    exterior4_gens((exterior4 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior4_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    exterior4 v = ((PyExterior4 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15]);
}

static PyObject *
pyexterior4_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_exterior4, 16);
}

static PyObject *
pyexterior4_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_exterior4, 16);
}

static PyMethodDef pyexterior4_methods[] = {
    {"to_float", (PyCFunction)pyexterior4_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert lattice values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyexterior4_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the lattice dtype."},
    {"zero", pyexterior4_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"bottom", pyexterior4_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pyexterior4_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"basis", pyexterior4_basis, METH_NOARGS | METH_CLASS, "Coefficient basis."},
    {"gens", pyexterior4_gens, METH_NOARGS | METH_CLASS, "Generators in basis order; excludes bottom."},
    {"ldual", pyexterior4_ldual, METH_NOARGS, "Left dual: a | a.ldual() == top for basis blades."},
    {"rdual", pyexterior4_rdual, METH_NOARGS, "Right dual: a.rdual() | a == top for basis blades."},
    {"norm2", pyexterior4_norm2, METH_NOARGS, "Sum of squared coefficients."},
    {"__reduce__", pyexterior4_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyexterior4_length(PyObject *Py_UNUSED(self))
{
    return 16;
}

static PyObject *
pyexterior4_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 16) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyExterior4 *)self)->v.data[i]);
}

static PyObject *
pyexterior4_subscript(PyObject *self, PyObject *key)
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
        i += 16;
    }
    return pyexterior4_item(self, i);
}

static int
pyexterior4_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyexterior4_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyexterior4_slots[] = {
    {Py_tp_dealloc, pyexterior4_dealloc},
    {Py_tp_traverse, pyexterior4_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyexterior4_repr},
    {Py_tp_str, pyexterior4_repr},
    {Py_tp_hash, pyexterior4_hash},
    {Py_tp_doc, (void *)"Exterior algebra on 4 generators as a lattice with 16 float64 coefficients"},
    {Py_tp_richcompare, pyexterior4_richcompare},
    {Py_tp_methods, pyexterior4_methods},
    {Py_tp_getset, pyexterior4_getset},
    {Py_sq_length, pyexterior4_length},
    {Py_sq_item, pyexterior4_item},
    {Py_mp_subscript, pyexterior4_subscript},
    {Py_tp_new, pyexterior4_new},
    {Py_nb_add, pyexterior4_add},
    {Py_nb_subtract, pyexterior4_sub},
    {Py_nb_multiply, pyexterior4_mul},
    {Py_nb_true_divide, pyexterior4_div},
    {Py_nb_or, pyexterior4_join},
    {Py_nb_and, pyexterior4_meet},
    {Py_nb_positive, pyexterior4_pos},
    {Py_nb_negative, pyexterior4_neg},
    {Py_nb_invert, pyexterior4_invert},
    {Py_nb_bool, pyexterior4_bool},
    {0, NULL}
};

static PyType_Spec pyexterior4_spec = {
    .name = "algpy.lattices.clifford.exterior4",
    .basicsize = sizeof(PyExterior4),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyexterior4_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyexterior4_getitem(void *data, void *arr)
{
    exterior4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_EXTERIOR4(v);
}

static int
npyexterior4_setitem(PyObject *item, void *data, void *arr)
{
    exterior4 v;
    int r = exterior4_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to exterior4");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 16) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "exterior4 needs exactly 16 components");
            return -1;
        }
        for (int i = 0; i < 16; i++) {
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
pyexterior4_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(exterior4));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyexterior4_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyexterior4_copyswapn(dst, sizeof(exterior4), src, sizeof(exterior4), 1,
                         swap, arr);
}

static npy_bool
pyexterior4_nonzero(void *data, void *arr)
{
    exterior4 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec16_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyexterior4_arrfuncs;
typedef struct { char c; exterior4 v; } pyexterior4_align_test;

static PyArray_DescrProto pyexterior4_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(exterior4), offsetof(pyexterior4_align_test, v), 0, 0, 0,
    &pyexterior4_arrfuncs,
};

#define EXTERIOR4_CAST_FROM(name, type, convert)                                \
    static void                                                              \
    pyexterior4_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        exterior4 *dst = (exterior4 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = exterior4_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define EXTERIOR4_REAL_TYPES(X) \
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

EXTERIOR4_REAL_TYPES(EXTERIOR4_CAST_FROM)

static int
pyexterior4_register_casts(void)
{
#define REGISTER_EXTERIOR4_CAST(name, type, convert)                            \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_exterior4,                    \
                                     pyexterior4_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_exterior4, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    EXTERIOR4_REAL_TYPES(REGISTER_EXTERIOR4_CAST)
#undef REGISTER_EXTERIOR4_CAST
    return 0;
}

typedef exterior4 (*exterior4_unary_op)(exterior4);
typedef exterior4 (*exterior4_binary_op)(exterior4, exterior4);
typedef int (*exterior4_predicate)(exterior4, exterior4);

static void
loop_exterior4_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    exterior4_unary_op op = (exterior4_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior4 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_exterior4_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior4_binary_op op = (exterior4_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior4 a, b, r;
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
loop_exterior4_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    exterior4_predicate op = (exterior4_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior4 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_exterior4_rmul(char **args, npy_intp const *dims,
                   npy_intp const *steps, void *Py_UNUSED(data))
{
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        double a;
        exterior4 b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = exterior4_rmul(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

typedef exterior4 (*exterior4_scalar_op)(exterior4, double);

static void
loop_exterior4_scalar(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior4_scalar_op op = (exterior4_scalar_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior4 a, r;
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
pyexterior4_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs,
                         const int *signature, void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = signature[i] < 0 ? npy_exterior4 : signature[i];
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_exterior4,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyexterior4_register_ufuncs(PyObject *numpy)
{
    static const int unary[2] = {-1, -1};
    static const int binary[3] = {-1, -1, -1};
    static const int predicate[3] = {-1, -1, NPY_BOOL};
    static const int rscalar[3] = {NPY_DOUBLE, -1, -1};
    static const int lscalar[3] = {-1, NPY_DOUBLE, -1};
    return
        pyexterior4_register_loop(numpy, "add", loop_exterior4_binary, 3, binary, exterior4_add) ||
        pyexterior4_register_loop(numpy, "subtract", loop_exterior4_binary, 3, binary, exterior4_sub) ||
        pyexterior4_register_loop(numpy, "multiply", loop_exterior4_rmul, 3, rscalar, NULL) ||
        pyexterior4_register_loop(numpy, "multiply", loop_exterior4_scalar, 3, lscalar, exterior4_lmul) ||
        pyexterior4_register_loop(numpy, "divide", loop_exterior4_scalar, 3, lscalar, exterior4_ldiv) ||
        pyexterior4_register_loop(numpy, "bitwise_or", loop_exterior4_binary, 3, binary, exterior4_join) ||
        pyexterior4_register_loop(numpy, "bitwise_and", loop_exterior4_binary, 3, binary, exterior4_meet) ||
        pyexterior4_register_loop(numpy, "positive", loop_exterior4_unary, 2, unary, exterior4_pos) ||
        pyexterior4_register_loop(numpy, "negative", loop_exterior4_unary, 2, unary, exterior4_neg) ||
        pyexterior4_register_loop(numpy, "equal", loop_exterior4_predicate, 3, predicate, exterior4_equal) ||
        pyexterior4_register_loop(numpy, "not_equal", loop_exterior4_predicate, 3, predicate, exterior4_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pyexterior4_register_lattice_ufuncs(PyObject *ufuncs)
{
    static const int unary[2] = {-1, -1};
    return
        pyexterior4_register_loop(ufuncs, "ldual", loop_exterior4_unary, 2, unary, exterior4_ldual) ||
        pyexterior4_register_loop(ufuncs, "rdual", loop_exterior4_unary, 2, unary, exterior4_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_exterior4", "exterior4 lattice type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__exterior4(void)
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
    PyExterior4_Type = PyType_FromModuleAndSpec(
        m, &pyexterior4_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyExterior4_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyexterior4_arrfuncs);
    pyexterior4_arrfuncs.getitem = npyexterior4_getitem;
    pyexterior4_arrfuncs.setitem = npyexterior4_setitem;
    pyexterior4_arrfuncs.copyswapn = pyexterior4_copyswapn;
    pyexterior4_arrfuncs.copyswap = pyexterior4_copyswap;
    pyexterior4_arrfuncs.nonzero = pyexterior4_nonzero;
    pyexterior4_descr_proto.typeobj = (PyTypeObject *)PyExterior4_Type;
    Py_SET_TYPE((PyObject *)&pyexterior4_descr_proto, &PyArrayDescr_Type);
    npy_exterior4 = PyArray_RegisterDataType(&pyexterior4_descr_proto);
    if (npy_exterior4 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_exterior4);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyExterior4_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyExterior4_Type, 16) < 0) {
        goto fail;
    }
    if (pyexterior4_register_casts() < 0 ||
            pyexterior4_register_ufuncs(numpy) < 0 ||
            pyexterior4_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "exterior4", PyExterior4_Type) < 0) {
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
