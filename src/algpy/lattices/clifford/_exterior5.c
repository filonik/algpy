/*
 * Exterior algebra on 5 generators exposed as a NumPy user-defined dtype with a
 * lattice interface: join (|), meet (&), ldual and rdual (also registered on
 * the shared algpy.lattices ufuncs) (and ~).
 * The arithmetic itself lives in exterior5.h.
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

#include "exterior5.h"

typedef struct {
    PyObject_HEAD
    exterior5 v;
} PyExterior5;

static PyObject *PyExterior5_Type = NULL;
static int npy_exterior5 = -1;

static PyObject *
PyExterior5_FromValue(PyTypeObject *type, exterior5 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyExterior5 *p = (PyExterior5 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_EXTERIOR5(v) \
    PyExterior5_FromValue((PyTypeObject *)PyExterior5_Type, (v))

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
exterior5_from_object(PyObject *o, exterior5 *out)
{
    double d;
    int r;
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyExterior5_Type)) {
        *out = ((PyExterior5 *)o)->v;
        return 1;
    }
    r = real_from_object(o, &d);
    if (r > 0) {
        *out = exterior5_from_re(d);
    }
    return r;
}

static PyObject *
pyexterior5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[32] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", "x3", "x03", "x13", "x013", "x23", "x023", "x123", "x0123", "x4", "x04", "x14", "x014", "x24", "x024", "x124", "x0124", "x34", "x034", "x134", "x0134", "x234", "x0234", "x1234", "x01234", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddddddddddddddddddddddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23], &d[24], &d[25], &d[26], &d[27], &d[28], &d[29], &d[30], &d[31])) {
        return NULL;
    }
    return PyExterior5_FromValue(type, exterior5_make(d));
}

static PyObject *
pyexterior5_repr(PyObject *self)
{
    exterior5 v = ((PyExterior5 *)self)->v;
    char *s[32] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 32; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("exterior5(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11], s[12], s[13], s[14], s[15], s[16], s[17], s[18], s[19], s[20], s[21], s[22], s[23], s[24], s[25], s[26], s[27], s[28], s[29], s[30], s[31]);
done:
    for (i = 0; i < 32; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyexterior5_hash(PyObject *self)
{
    exterior5 v = ((PyExterior5 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_EXTERIOR5_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (exterior5_is_real(v)) {
        obj = PyFloat_FromDouble(v.data[0]);
    }
    else {
        obj = Py_BuildValue("(dddddddddddddddddddddddddddddddd)", v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
    }
    if (obj == NULL) {
        return -1;
    }
    result = PyObject_Hash(obj);
    Py_DECREF(obj);
    return result;
}

static PyObject *
pyexterior5_richcompare(PyObject *a, PyObject *b, int op)
{
    exterior5 x, y;
    int r;
    if ((r = exterior5_from_object(a, &x)) <= 0 ||
            (r = exterior5_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(exterior5_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!exterior5_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define EXTERIOR5_BINARY(name, expr)                                            \
    static PyObject *                                                        \
    pyexterior5_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        exterior5 x, y;                                                       \
        int r;                                                               \
        if ((r = exterior5_from_object(a, &x)) <= 0 ||                        \
                (r = exterior5_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_EXTERIOR5(expr(x, y));                                     \
    }

EXTERIOR5_BINARY(add, exterior5_add)
EXTERIOR5_BINARY(sub, exterior5_sub)
EXTERIOR5_BINARY(join, exterior5_join)
EXTERIOR5_BINARY(meet, exterior5_meet)

/* Multiplication and division only scale by real numbers. */
static PyObject *
pyexterior5_mul(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (PyObject_TypeCheck(a, (PyTypeObject *)PyExterior5_Type)) {
        if ((r = real_from_object(b, &d)) <= 0) {
            if (r < 0) {
                return NULL;
            }
            Py_RETURN_NOTIMPLEMENTED;
        }
        return NEW_EXTERIOR5(exterior5_lmul(((PyExterior5 *)a)->v, d));
    }
    if ((r = real_from_object(a, &d)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_EXTERIOR5(exterior5_rmul(d, ((PyExterior5 *)b)->v));
}

static PyObject *
pyexterior5_div(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (!PyObject_TypeCheck(a, (PyTypeObject *)PyExterior5_Type) ||
            (r = real_from_object(b, &d)) == 0) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (r < 0) {
        return NULL;
    }
    return NEW_EXTERIOR5(exterior5_ldiv(((PyExterior5 *)a)->v, d));
}

#define EXTERIOR5_UNARY(name, expr)                                             \
    static PyObject *                                                        \
    pyexterior5_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_EXTERIOR5(expr(((PyExterior5 *)self)->v));                  \
    }

EXTERIOR5_UNARY(pos, exterior5_pos)
EXTERIOR5_UNARY(neg, exterior5_neg)
EXTERIOR5_UNARY(invert, exterior5_dual)
static int
pyexterior5_bool(PyObject *self)
{
    exterior5 v = ((PyExterior5 *)self)->v;
    return vec32_any_nonzero(v);
}

#define EXTERIOR5_METHOD(name, op)                                              \
    static PyObject *                                                        \
    pyexterior5_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return PyExterior5_FromValue(Py_TYPE(self), op(((PyExterior5 *)self)->v));    \
    }

EXTERIOR5_METHOD(ldual, exterior5_ldual)
EXTERIOR5_METHOD(rdual, exterior5_rdual)
EXTERIOR5_METHOD(dual, exterior5_dual)
#define EXTERIOR5_CONSTANT(name, value)                                         \
    static PyObject *                                                        \
    pyexterior5_##name(PyObject *cls, PyObject *Py_UNUSED(args))              \
    {                                                                        \
        return PyExterior5_FromValue((PyTypeObject *)cls, value());                \
    }

EXTERIOR5_CONSTANT(zero, exterior5_zero)
EXTERIOR5_CONSTANT(bottom, exterior5_bottom)
EXTERIOR5_CONSTANT(top, exterior5_top)

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyexterior5_grade(PyObject *self, Py_ssize_t k)
{
    exterior5 v = ((PyExterior5 *)self)->v;
    switch (k) {
        case 0:
            return NEW_EXTERIOR5(exterior5_grade0(v));
        case 1:
            return NEW_EXTERIOR5(exterior5_grade1(v));
        case 2:
            return NEW_EXTERIOR5(exterior5_grade2(v));
        case 3:
            return NEW_EXTERIOR5(exterior5_grade3(v));
        case 4:
            return NEW_EXTERIOR5(exterior5_grade4(v));
        case 5:
            return NEW_EXTERIOR5(exterior5_grade5(v));
        default:
            return NEW_EXTERIOR5(exterior5_zero());
    }
}

static PyObject *
pyexterior5_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyexterior5_grade, 6);
}

static PyGetSetDef pyexterior5_getset[] = {
    {"grade", pyexterior5_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyexterior5_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(exterior5_norm2(((PyExterior5 *)self)->v));
}

static PyObject *
pyexterior5_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {32};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior5);
    if (arr == NULL) {
        return NULL;
    }
    exterior5_basis((exterior5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior5_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {5};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior5);
    if (arr == NULL) {
        return NULL;
    }
    exterior5_gens((exterior5 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior5_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    exterior5 v = ((PyExterior5 *)self)->v;
    return Py_BuildValue("(O(dddddddddddddddddddddddddddddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7], v.data[8], v.data[9], v.data[10], v.data[11], v.data[12], v.data[13], v.data[14], v.data[15], v.data[16], v.data[17], v.data[18], v.data[19], v.data[20], v.data[21], v.data[22], v.data[23], v.data[24], v.data[25], v.data[26], v.data[27], v.data[28], v.data[29], v.data[30], v.data[31]);
}

static PyObject *
pyexterior5_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_exterior5, 32);
}

static PyObject *
pyexterior5_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_exterior5, 32);
}

static PyMethodDef pyexterior5_methods[] = {
    {"to_float", (PyCFunction)pyexterior5_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert lattice values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyexterior5_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the lattice dtype."},
    {"zero", pyexterior5_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"bottom", pyexterior5_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pyexterior5_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"basis", pyexterior5_basis, METH_NOARGS | METH_CLASS, "Coefficient basis."},
    {"gens", pyexterior5_gens, METH_NOARGS | METH_CLASS, "Generators in basis order; excludes bottom."},
    {"ldual", pyexterior5_ldual, METH_NOARGS, "Left dual: a | a.ldual() == top for basis blades."},
    {"rdual", pyexterior5_rdual, METH_NOARGS, "Right dual: a.rdual() | a == top for basis blades."},
    {"dual", pyexterior5_dual, METH_NOARGS, "Dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"norm2", pyexterior5_norm2, METH_NOARGS, "Sum of squared coefficients."},
    {"__reduce__", pyexterior5_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyexterior5_length(PyObject *Py_UNUSED(self))
{
    return 32;
}

static PyObject *
pyexterior5_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 32) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyExterior5 *)self)->v.data[i]);
}

static PyObject *
pyexterior5_subscript(PyObject *self, PyObject *key)
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
        i += 32;
    }
    return pyexterior5_item(self, i);
}

static int
pyexterior5_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyexterior5_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyexterior5_slots[] = {
    {Py_tp_dealloc, pyexterior5_dealloc},
    {Py_tp_traverse, pyexterior5_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyexterior5_repr},
    {Py_tp_str, pyexterior5_repr},
    {Py_tp_hash, pyexterior5_hash},
    {Py_tp_doc, (void *)"Exterior algebra on 5 generators as a lattice with 32 float64 coefficients"},
    {Py_tp_richcompare, pyexterior5_richcompare},
    {Py_tp_methods, pyexterior5_methods},
    {Py_tp_getset, pyexterior5_getset},
    {Py_sq_length, pyexterior5_length},
    {Py_sq_item, pyexterior5_item},
    {Py_mp_subscript, pyexterior5_subscript},
    {Py_tp_new, pyexterior5_new},
    {Py_nb_add, pyexterior5_add},
    {Py_nb_subtract, pyexterior5_sub},
    {Py_nb_multiply, pyexterior5_mul},
    {Py_nb_true_divide, pyexterior5_div},
    {Py_nb_or, pyexterior5_join},
    {Py_nb_and, pyexterior5_meet},
    {Py_nb_positive, pyexterior5_pos},
    {Py_nb_negative, pyexterior5_neg},
    {Py_nb_invert, pyexterior5_invert},
    {Py_nb_bool, pyexterior5_bool},
    {0, NULL}
};

static PyType_Spec pyexterior5_spec = {
    .name = "algpy.lattices.clifford.exterior5",
    .basicsize = sizeof(PyExterior5),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyexterior5_slots,
};

static void
swap_bytes(char *p)
{
    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 4; j++) {
            char t = p[8 * i + j];
            p[8 * i + j] = p[8 * i + 7 - j];
            p[8 * i + 7 - j] = t;
        }
    }
}

static PyObject *
npyexterior5_getitem(void *data, void *arr)
{
    exterior5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_EXTERIOR5(v);
}

static int
npyexterior5_setitem(PyObject *item, void *data, void *arr)
{
    exterior5 v;
    int r = exterior5_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to exterior5");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 32) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "exterior5 needs exactly 32 components");
            return -1;
        }
        for (int i = 0; i < 32; i++) {
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
pyexterior5_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(exterior5));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyexterior5_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyexterior5_copyswapn(dst, sizeof(exterior5), src, sizeof(exterior5), 1,
                         swap, arr);
}

static npy_bool
pyexterior5_nonzero(void *data, void *arr)
{
    exterior5 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec32_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyexterior5_arrfuncs;
typedef struct { char c; exterior5 v; } pyexterior5_align_test;

static PyArray_DescrProto pyexterior5_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(exterior5), offsetof(pyexterior5_align_test, v), 0, 0, 0,
    &pyexterior5_arrfuncs,
};

#define EXTERIOR5_CAST_FROM(name, type, convert)                                \
    static void                                                              \
    pyexterior5_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        exterior5 *dst = (exterior5 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = exterior5_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define EXTERIOR5_REAL_TYPES(X) \
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

EXTERIOR5_REAL_TYPES(EXTERIOR5_CAST_FROM)

static int
pyexterior5_register_casts(void)
{
#define REGISTER_EXTERIOR5_CAST(name, type, convert)                            \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_exterior5,                    \
                                     pyexterior5_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_exterior5, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    EXTERIOR5_REAL_TYPES(REGISTER_EXTERIOR5_CAST)
#undef REGISTER_EXTERIOR5_CAST
    return 0;
}

typedef exterior5 (*exterior5_unary_op)(exterior5);
typedef exterior5 (*exterior5_binary_op)(exterior5, exterior5);
typedef int (*exterior5_predicate)(exterior5, exterior5);

static void
loop_exterior5_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    exterior5_unary_op op = (exterior5_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior5 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_exterior5_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior5_binary_op op = (exterior5_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior5 a, b, r;
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
loop_exterior5_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    exterior5_predicate op = (exterior5_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior5 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_exterior5_rmul(char **args, npy_intp const *dims,
                   npy_intp const *steps, void *Py_UNUSED(data))
{
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        double a;
        exterior5 b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = exterior5_rmul(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

typedef exterior5 (*exterior5_scalar_op)(exterior5, double);

static void
loop_exterior5_scalar(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior5_scalar_op op = (exterior5_scalar_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior5 a, r;
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
pyexterior5_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs,
                         const int *signature, void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = signature[i] < 0 ? npy_exterior5 : signature[i];
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_exterior5,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyexterior5_register_ufuncs(PyObject *numpy)
{
    static const int unary[2] = {-1, -1};
    static const int binary[3] = {-1, -1, -1};
    static const int predicate[3] = {-1, -1, NPY_BOOL};
    static const int rscalar[3] = {NPY_DOUBLE, -1, -1};
    static const int lscalar[3] = {-1, NPY_DOUBLE, -1};
    return
        pyexterior5_register_loop(numpy, "add", loop_exterior5_binary, 3, binary, exterior5_add) ||
        pyexterior5_register_loop(numpy, "subtract", loop_exterior5_binary, 3, binary, exterior5_sub) ||
        pyexterior5_register_loop(numpy, "multiply", loop_exterior5_rmul, 3, rscalar, NULL) ||
        pyexterior5_register_loop(numpy, "multiply", loop_exterior5_scalar, 3, lscalar, exterior5_lmul) ||
        pyexterior5_register_loop(numpy, "divide", loop_exterior5_scalar, 3, lscalar, exterior5_ldiv) ||
        pyexterior5_register_loop(numpy, "bitwise_or", loop_exterior5_binary, 3, binary, exterior5_join) ||
        pyexterior5_register_loop(numpy, "bitwise_and", loop_exterior5_binary, 3, binary, exterior5_meet) ||
        pyexterior5_register_loop(numpy, "invert", loop_exterior5_unary, 2, unary, exterior5_dual) ||
        pyexterior5_register_loop(numpy, "positive", loop_exterior5_unary, 2, unary, exterior5_pos) ||
        pyexterior5_register_loop(numpy, "negative", loop_exterior5_unary, 2, unary, exterior5_neg) ||
        pyexterior5_register_loop(numpy, "equal", loop_exterior5_predicate, 3, predicate, exterior5_equal) ||
        pyexterior5_register_loop(numpy, "not_equal", loop_exterior5_predicate, 3, predicate, exterior5_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pyexterior5_register_lattice_ufuncs(PyObject *ufuncs)
{
    static const int unary[2] = {-1, -1};
    return
        pyexterior5_register_loop(ufuncs, "ldual", loop_exterior5_unary, 2, unary, exterior5_ldual) ||
        pyexterior5_register_loop(ufuncs, "rdual", loop_exterior5_unary, 2, unary, exterior5_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_exterior5", "exterior5 lattice type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__exterior5(void)
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
    PyExterior5_Type = PyType_FromModuleAndSpec(
        m, &pyexterior5_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyExterior5_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyexterior5_arrfuncs);
    pyexterior5_arrfuncs.getitem = npyexterior5_getitem;
    pyexterior5_arrfuncs.setitem = npyexterior5_setitem;
    pyexterior5_arrfuncs.copyswapn = pyexterior5_copyswapn;
    pyexterior5_arrfuncs.copyswap = pyexterior5_copyswap;
    pyexterior5_arrfuncs.nonzero = pyexterior5_nonzero;
    pyexterior5_descr_proto.typeobj = (PyTypeObject *)PyExterior5_Type;
    Py_SET_TYPE((PyObject *)&pyexterior5_descr_proto, &PyArrayDescr_Type);
    npy_exterior5 = PyArray_RegisterDataType(&pyexterior5_descr_proto);
    if (npy_exterior5 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_exterior5);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyExterior5_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyExterior5_Type, 32) < 0) {
        goto fail;
    }
    if (pyexterior5_register_casts() < 0 ||
            pyexterior5_register_ufuncs(numpy) < 0 ||
            pyexterior5_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "exterior5", PyExterior5_Type) < 0) {
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
