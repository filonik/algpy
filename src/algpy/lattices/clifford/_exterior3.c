/*
 * Exterior algebra on 3 generators exposed as a NumPy user-defined dtype with a
 * lattice interface: join (|), meet (&), ldual and rdual (also registered on
 * the shared algpy.lattices ufuncs) (and ~).
 * The arithmetic itself lives in exterior3.h.
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

#include "exterior3.h"

typedef struct {
    PyObject_HEAD
    exterior3 v;
} PyExterior3;

static PyObject *PyExterior3_Type = NULL;
static int npy_exterior3 = -1;

static PyObject *
PyExterior3_FromValue(PyTypeObject *type, exterior3 v)
{
    allocfunc tp_alloc = (allocfunc)PyType_GetSlot(type, Py_tp_alloc);
    PyExterior3 *p = (PyExterior3 *)tp_alloc(type, 0);
    if (p != NULL) {
        p->v = v;
    }
    return (PyObject *)p;
}

#define NEW_EXTERIOR3(v) \
    PyExterior3_FromValue((PyTypeObject *)PyExterior3_Type, (v))

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
exterior3_from_object(PyObject *o, exterior3 *out)
{
    double d;
    int r;
    if (PyObject_TypeCheck(o, (PyTypeObject *)PyExterior3_Type)) {
        *out = ((PyExterior3 *)o)->v;
        return 1;
    }
    r = real_from_object(o, &d);
    if (r > 0) {
        *out = exterior3_from_re(d);
    }
    return r;
}

static PyObject *
pyexterior3_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    double d[8] = {0};
    static char *kwlist[] = {"x", "x0", "x1", "x01", "x2", "x02", "x12", "x012", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddddddd", kwlist,
            &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7])) {
        return NULL;
    }
    return PyExterior3_FromValue(type, exterior3_make(d));
}

static PyObject *
pyexterior3_repr(PyObject *self)
{
    exterior3 v = ((PyExterior3 *)self)->v;
    char *s[8] = {0};
    PyObject *result = NULL;
    int i;
    for (i = 0; i < 8; i++) {
        s[i] = PyOS_double_to_string(v.data[i], 'g', 6, 0, NULL);
        if (s[i] == NULL) {
            goto done;
        }
    }
    result = PyUnicode_FromFormat("exterior3(%s, %s, %s, %s, %s, %s, %s, %s)",
                                  s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
done:
    for (i = 0; i < 8; i++) {
        PyMem_Free(s[i]);
    }
    return result;
}

static Py_hash_t
pyexterior3_hash(PyObject *self)
{
    exterior3 v = ((PyExterior3 *)self)->v;
    PyObject *obj;
    Py_hash_t result;
    int i;
    /* NaN hashes must not depend on temporary Python float allocations. */
    for (i = 0; i < ALGPY_EXTERIOR3_N; i++) {
        if (isnan(v.data[i])) {
            return 0;
        }
    }
    if (exterior3_is_real(v)) {
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
pyexterior3_richcompare(PyObject *a, PyObject *b, int op)
{
    exterior3 x, y;
    int r;
    if ((r = exterior3_from_object(a, &x)) <= 0 ||
            (r = exterior3_from_object(b, &y)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    switch (op) {
        case Py_EQ:
            return PyBool_FromLong(exterior3_equal(x, y));
        case Py_NE:
            return PyBool_FromLong(!exterior3_equal(x, y));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

#define EXTERIOR3_BINARY(name, expr)                                            \
    static PyObject *                                                        \
    pyexterior3_##name(PyObject *a, PyObject *b)                              \
    {                                                                        \
        exterior3 x, y;                                                       \
        int r;                                                               \
        if ((r = exterior3_from_object(a, &x)) <= 0 ||                        \
                (r = exterior3_from_object(b, &y)) <= 0) {                    \
            if (r < 0) {                                                     \
                return NULL;                                                 \
            }                                                                \
            Py_RETURN_NOTIMPLEMENTED;                                        \
        }                                                                    \
        return NEW_EXTERIOR3(expr(x, y));                                     \
    }

EXTERIOR3_BINARY(add, exterior3_add)
EXTERIOR3_BINARY(sub, exterior3_sub)
EXTERIOR3_BINARY(join, exterior3_join)
EXTERIOR3_BINARY(meet, exterior3_meet)

/* Multiplication and division only scale by real numbers. */
static PyObject *
pyexterior3_mul(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (PyObject_TypeCheck(a, (PyTypeObject *)PyExterior3_Type)) {
        if ((r = real_from_object(b, &d)) <= 0) {
            if (r < 0) {
                return NULL;
            }
            Py_RETURN_NOTIMPLEMENTED;
        }
        return NEW_EXTERIOR3(exterior3_lmul(((PyExterior3 *)a)->v, d));
    }
    if ((r = real_from_object(a, &d)) <= 0) {
        if (r < 0) {
            return NULL;
        }
        Py_RETURN_NOTIMPLEMENTED;
    }
    return NEW_EXTERIOR3(exterior3_rmul(d, ((PyExterior3 *)b)->v));
}

static PyObject *
pyexterior3_div(PyObject *a, PyObject *b)
{
    double d;
    int r;
    if (!PyObject_TypeCheck(a, (PyTypeObject *)PyExterior3_Type) ||
            (r = real_from_object(b, &d)) == 0) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (r < 0) {
        return NULL;
    }
    return NEW_EXTERIOR3(exterior3_ldiv(((PyExterior3 *)a)->v, d));
}

#define EXTERIOR3_UNARY(name, expr)                                             \
    static PyObject *                                                        \
    pyexterior3_##name(PyObject *self)                                        \
    {                                                                        \
        return NEW_EXTERIOR3(expr(((PyExterior3 *)self)->v));                  \
    }

EXTERIOR3_UNARY(pos, exterior3_pos)
EXTERIOR3_UNARY(neg, exterior3_neg)
EXTERIOR3_UNARY(invert, exterior3_dual)
static int
pyexterior3_bool(PyObject *self)
{
    exterior3 v = ((PyExterior3 *)self)->v;
    return vec8_any_nonzero(v);
}

#define EXTERIOR3_METHOD(name, op)                                              \
    static PyObject *                                                        \
    pyexterior3_##name(PyObject *self, PyObject *Py_UNUSED(args))             \
    {                                                                        \
        return PyExterior3_FromValue(Py_TYPE(self), op(((PyExterior3 *)self)->v));    \
    }

EXTERIOR3_METHOD(ldual, exterior3_ldual)
EXTERIOR3_METHOD(rdual, exterior3_rdual)
EXTERIOR3_METHOD(dual, exterior3_dual)
#define EXTERIOR3_CONSTANT(name, value)                                         \
    static PyObject *                                                        \
    pyexterior3_##name(PyObject *cls, PyObject *Py_UNUSED(args))              \
    {                                                                        \
        return PyExterior3_FromValue((PyTypeObject *)cls, value());                \
    }

EXTERIOR3_CONSTANT(zero, exterior3_zero)
EXTERIOR3_CONSTANT(bottom, exterior3_bottom)
EXTERIOR3_CONSTANT(top, exterior3_top)

/* Projection onto grade k, for 0 <= k <= generators; used by x.grade. */
static PyObject *
pyexterior3_grade(PyObject *self, Py_ssize_t k)
{
    exterior3 v = ((PyExterior3 *)self)->v;
    switch (k) {
        case 0:
            return NEW_EXTERIOR3(exterior3_grade0(v));
        case 1:
            return NEW_EXTERIOR3(exterior3_grade1(v));
        case 2:
            return NEW_EXTERIOR3(exterior3_grade2(v));
        case 3:
            return NEW_EXTERIOR3(exterior3_grade3(v));
        default:
            return NEW_EXTERIOR3(exterior3_zero());
    }
}

static PyObject *
pyexterior3_grades(PyObject *self, void *Py_UNUSED(closure))
{
    return algpy_grades_view(self, pyexterior3_grade, 4);
}

static PyGetSetDef pyexterior3_getset[] = {
    {"grade", pyexterior3_grades, NULL, "Read-only grade projections: x.grade[k] for 0 <= k < len(x.grade).", NULL},
    {NULL, NULL, NULL, NULL, NULL}
};

static PyObject *
pyexterior3_norm2(PyObject *self, PyObject *Py_UNUSED(args))
{
    return PyFloat_FromDouble(exterior3_norm2(((PyExterior3 *)self)->v));
}

static PyObject *
pyexterior3_basis(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {8};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior3);
    if (arr == NULL) {
        return NULL;
    }
    exterior3_basis((exterior3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior3_gens(PyObject *Py_UNUSED(cls), PyObject *Py_UNUSED(args))
{
    npy_intp dims[1] = {3};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(
        1, dims, npy_exterior3);
    if (arr == NULL) {
        return NULL;
    }
    exterior3_gens((exterior3 *)PyArray_DATA(arr));
    return (PyObject *)arr;
}

static PyObject *
pyexterior3_reduce(PyObject *self, PyObject *Py_UNUSED(args))
{
    exterior3 v = ((PyExterior3 *)self)->v;
    return Py_BuildValue("(O(dddddddd))", (PyObject *)Py_TYPE(self),
                         v.data[0], v.data[1], v.data[2], v.data[3], v.data[4], v.data[5], v.data[6], v.data[7]);
}

static PyObject *
pyexterior3_to_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_to_float(args, kwargs, npy_exterior3, 8);
}

static PyObject *
pyexterior3_from_float(PyObject *Py_UNUSED(cls), PyObject *args, PyObject *kwargs)
{
    return algpy_from_float(args, kwargs, npy_exterior3, 8);
}

static PyMethodDef pyexterior3_methods[] = {
    {"to_float", (PyCFunction)pyexterior3_to_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert lattice values to real coefficient arrays, with optional dtype."},
    {"from_float", (PyCFunction)pyexterior3_from_float, METH_VARARGS | METH_KEYWORDS | METH_CLASS, "Convert arrays with a final coefficient axis to the lattice dtype."},
    {"zero", pyexterior3_zero, METH_NOARGS | METH_CLASS, "Additive identity."},
    {"bottom", pyexterior3_bottom, METH_NOARGS | METH_CLASS, "Identity of join (|): the scalar 1."},
    {"top", pyexterior3_top, METH_NOARGS | METH_CLASS, "Identity of meet (&): the pseudoscalar."},
    {"basis", pyexterior3_basis, METH_NOARGS | METH_CLASS, "Coefficient basis."},
    {"gens", pyexterior3_gens, METH_NOARGS | METH_CLASS, "Generators in basis order; excludes bottom."},
    {"ldual", pyexterior3_ldual, METH_NOARGS, "Left dual: a | a.ldual() == top for basis blades."},
    {"rdual", pyexterior3_rdual, METH_NOARGS, "Right dual: a.rdual() | a == top for basis blades."},
    {"dual", pyexterior3_dual, METH_NOARGS, "Dual; ldual and rdual coincide in odd dimensions. Also ~a."},
    {"norm2", pyexterior3_norm2, METH_NOARGS, "Sum of squared coefficients."},
    {"__reduce__", pyexterior3_reduce, METH_NOARGS, "Support for pickle."},
    {NULL, NULL, 0, NULL}
};

static Py_ssize_t
pyexterior3_length(PyObject *Py_UNUSED(self))
{
    return 8;
}

static PyObject *
pyexterior3_item(PyObject *self, Py_ssize_t i)
{
    if (i < 0 || i >= 8) {
        PyErr_SetString(PyExc_IndexError, "index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(((PyExterior3 *)self)->v.data[i]);
}

static PyObject *
pyexterior3_subscript(PyObject *self, PyObject *key)
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
    return pyexterior3_item(self, i);
}

static int
pyexterior3_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    return 0;
}

static void
pyexterior3_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);
    PyObject_GC_UnTrack(self);
    tp_free(self);
    Py_DECREF(tp);
}

static PyType_Slot pyexterior3_slots[] = {
    {Py_tp_dealloc, pyexterior3_dealloc},
    {Py_tp_traverse, pyexterior3_traverse},
    {Py_tp_alloc, PyType_GenericAlloc},
    {Py_tp_free, PyObject_GC_Del},
    {Py_tp_repr, pyexterior3_repr},
    {Py_tp_str, pyexterior3_repr},
    {Py_tp_hash, pyexterior3_hash},
    {Py_tp_doc, (void *)"Exterior algebra on 3 generators as a lattice with 8 float64 coefficients"},
    {Py_tp_richcompare, pyexterior3_richcompare},
    {Py_tp_methods, pyexterior3_methods},
    {Py_tp_getset, pyexterior3_getset},
    {Py_sq_length, pyexterior3_length},
    {Py_sq_item, pyexterior3_item},
    {Py_mp_subscript, pyexterior3_subscript},
    {Py_tp_new, pyexterior3_new},
    {Py_nb_add, pyexterior3_add},
    {Py_nb_subtract, pyexterior3_sub},
    {Py_nb_multiply, pyexterior3_mul},
    {Py_nb_true_divide, pyexterior3_div},
    {Py_nb_or, pyexterior3_join},
    {Py_nb_and, pyexterior3_meet},
    {Py_nb_positive, pyexterior3_pos},
    {Py_nb_negative, pyexterior3_neg},
    {Py_nb_invert, pyexterior3_invert},
    {Py_nb_bool, pyexterior3_bool},
    {0, NULL}
};

static PyType_Spec pyexterior3_spec = {
    .name = "algpy.lattices.clifford.exterior3",
    .basicsize = sizeof(PyExterior3),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .slots = pyexterior3_slots,
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
npyexterior3_getitem(void *data, void *arr)
{
    exterior3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return NEW_EXTERIOR3(v);
}

static int
npyexterior3_setitem(PyObject *item, void *data, void *arr)
{
    exterior3 v;
    int r = exterior3_from_object(item, &v);
    if (r < 0) {
        return -1;
    }
    if (r == 0) {
        PyObject *seq = PySequence_Fast(item, "cannot convert to exterior3");
        if (seq == NULL) {
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(seq) != 8) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError,
                            "exterior3 needs exactly 8 components");
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
pyexterior3_copyswapn(void *dst, npy_intp dstride, void *src,
                     npy_intp sstride, npy_intp n, int swap, void *arr)
{
    npy_intp i;
    if (src != NULL) {
        for (i = 0; i < n; i++) {
            memcpy((char *)dst + i * dstride, (char *)src + i * sstride,
                   sizeof(exterior3));
        }
    }
    if (swap) {
        for (i = 0; i < n; i++) {
            swap_bytes((char *)dst + i * dstride);
        }
    }
}

static void
pyexterior3_copyswap(void *dst, void *src, int swap, void *arr)
{
    pyexterior3_copyswapn(dst, sizeof(exterior3), src, sizeof(exterior3), 1,
                         swap, arr);
}

static npy_bool
pyexterior3_nonzero(void *data, void *arr)
{
    exterior3 v;
    memcpy(&v, data, sizeof(v));
    if (arr != NULL && PyArray_ISBYTESWAPPED((PyArrayObject *)arr)) {
        swap_bytes((char *)&v);
    }
    return vec8_any_nonzero(v) ? NPY_TRUE : NPY_FALSE;
}

static PyArray_ArrFuncs pyexterior3_arrfuncs;
typedef struct { char c; exterior3 v; } pyexterior3_align_test;

static PyArray_DescrProto pyexterior3_descr_proto = {
    PyObject_HEAD_INIT(0)
    NULL, 'V', 'c', '=', NPY_NEEDS_PYAPI | NPY_USE_GETITEM | NPY_USE_SETITEM,
    0, sizeof(exterior3), offsetof(pyexterior3_align_test, v), 0, 0, 0,
    &pyexterior3_arrfuncs,
};

#define EXTERIOR3_CAST_FROM(name, type, convert)                                \
    static void                                                              \
    pyexterior3_cast_##name(void *from, void *to, npy_intp n,                 \
                           void *fromarr, void *toarr)                       \
    {                                                                        \
        const type *src = (const type *)from;                               \
        exterior3 *dst = (exterior3 *)to;                                     \
        npy_intp i;                                                          \
        for (i = 0; i < n; i++) {                                            \
            dst[i] = exterior3_from_re(convert(src[i]));                      \
        }                                                                    \
    }

#define EXTERIOR3_REAL_TYPES(X) \
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

EXTERIOR3_REAL_TYPES(EXTERIOR3_CAST_FROM)

static int
pyexterior3_register_casts(void)
{
#define REGISTER_EXTERIOR3_CAST(name, type, convert)                            \
    {                                                                        \
        PyArray_Descr *from = PyArray_DescrFromType(NPY_##name);             \
        if (from == NULL) {                                                  \
            return -1;                                                       \
        }                                                                    \
        if (PyArray_RegisterCastFunc(from, npy_exterior3,                    \
                                     pyexterior3_cast_##name) < 0 ||          \
                PyArray_RegisterCanCast(from, npy_exterior3, NPY_NOSCALAR) < 0) {\
            Py_DECREF(from);                                                 \
            return -1;                                                       \
        }                                                                    \
        Py_DECREF(from);                                                     \
    }
    EXTERIOR3_REAL_TYPES(REGISTER_EXTERIOR3_CAST)
#undef REGISTER_EXTERIOR3_CAST
    return 0;
}

typedef exterior3 (*exterior3_unary_op)(exterior3);
typedef exterior3 (*exterior3_binary_op)(exterior3, exterior3);
typedef int (*exterior3_predicate)(exterior3, exterior3);

static void
loop_exterior3_unary(char **args, npy_intp const *dims,
                    npy_intp const *steps, void *data)
{
    exterior3_unary_op op = (exterior3_unary_op)data;
    char *in = args[0], *out = args[1];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior3 a, r;
        memcpy(&a, in, sizeof(a));
        r = op(a);
        memcpy(out, &r, sizeof(r));
        in += steps[0];
        out += steps[1];
    }
}

static void
loop_exterior3_binary(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior3_binary_op op = (exterior3_binary_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior3 a, b, r;
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
loop_exterior3_predicate(char **args, npy_intp const *dims,
                        npy_intp const *steps, void *data)
{
    exterior3_predicate op = (exterior3_predicate)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior3 a, b;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        *(npy_bool *)out = op(a, b) ? NPY_TRUE : NPY_FALSE;
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

static void
loop_exterior3_rmul(char **args, npy_intp const *dims,
                   npy_intp const *steps, void *Py_UNUSED(data))
{
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        double a;
        exterior3 b, r;
        memcpy(&a, in1, sizeof(a));
        memcpy(&b, in2, sizeof(b));
        r = exterior3_rmul(a, b);
        memcpy(out, &r, sizeof(r));
        in1 += steps[0];
        in2 += steps[1];
        out += steps[2];
    }
}

typedef exterior3 (*exterior3_scalar_op)(exterior3, double);

static void
loop_exterior3_scalar(char **args, npy_intp const *dims,
                     npy_intp const *steps, void *data)
{
    exterior3_scalar_op op = (exterior3_scalar_op)data;
    char *in1 = args[0], *in2 = args[1], *out = args[2];
    npy_intp i;
    for (i = 0; i < dims[0]; i++) {
        exterior3 a, r;
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
pyexterior3_register_loop(PyObject *numpy, const char *name,
                         PyUFuncGenericFunction loop, int nargs,
                         const int *signature, void *data)
{
    int types[3], i, result;
    PyObject *ufunc = PyObject_GetAttrString(numpy, name);
    if (ufunc == NULL) {
        return -1;
    }
    for (i = 0; i < nargs; i++) {
        types[i] = signature[i] < 0 ? npy_exterior3 : signature[i];
    }
    result = PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, npy_exterior3,
                                         loop, types, data);
    Py_DECREF(ufunc);
    return result;
}

static int
pyexterior3_register_ufuncs(PyObject *numpy)
{
    static const int unary[2] = {-1, -1};
    static const int binary[3] = {-1, -1, -1};
    static const int predicate[3] = {-1, -1, NPY_BOOL};
    static const int rscalar[3] = {NPY_DOUBLE, -1, -1};
    static const int lscalar[3] = {-1, NPY_DOUBLE, -1};
    return
        pyexterior3_register_loop(numpy, "add", loop_exterior3_binary, 3, binary, exterior3_add) ||
        pyexterior3_register_loop(numpy, "subtract", loop_exterior3_binary, 3, binary, exterior3_sub) ||
        pyexterior3_register_loop(numpy, "multiply", loop_exterior3_rmul, 3, rscalar, NULL) ||
        pyexterior3_register_loop(numpy, "multiply", loop_exterior3_scalar, 3, lscalar, exterior3_lmul) ||
        pyexterior3_register_loop(numpy, "divide", loop_exterior3_scalar, 3, lscalar, exterior3_ldiv) ||
        pyexterior3_register_loop(numpy, "bitwise_or", loop_exterior3_binary, 3, binary, exterior3_join) ||
        pyexterior3_register_loop(numpy, "bitwise_and", loop_exterior3_binary, 3, binary, exterior3_meet) ||
        pyexterior3_register_loop(numpy, "invert", loop_exterior3_unary, 2, unary, exterior3_dual) ||
        pyexterior3_register_loop(numpy, "positive", loop_exterior3_unary, 2, unary, exterior3_pos) ||
        pyexterior3_register_loop(numpy, "negative", loop_exterior3_unary, 2, unary, exterior3_neg) ||
        pyexterior3_register_loop(numpy, "equal", loop_exterior3_predicate, 3, predicate, exterior3_equal) ||
        pyexterior3_register_loop(numpy, "not_equal", loop_exterior3_predicate, 3, predicate, exterior3_not_equal)
        ? -1 : 0;
}

/* Loops for the shared algpy.lattices ldual and rdual ufuncs. */
static int
pyexterior3_register_lattice_ufuncs(PyObject *ufuncs)
{
    static const int unary[2] = {-1, -1};
    return
        pyexterior3_register_loop(ufuncs, "ldual", loop_exterior3_unary, 2, unary, exterior3_ldual) ||
        pyexterior3_register_loop(ufuncs, "rdual", loop_exterior3_unary, 2, unary, exterior3_rdual)
        ? -1 : 0;
}

static PyMethodDef module_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_exterior3", "exterior3 lattice type for NumPy",
    -1, module_methods,
};

PyMODINIT_FUNC
PyInit__exterior3(void)
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
    PyExterior3_Type = PyType_FromModuleAndSpec(
        m, &pyexterior3_spec, (PyObject *)&PyGenericArrType_Type);
    if (PyExterior3_Type == NULL) {
        goto fail;
    }
    PyArray_InitArrFuncs(&pyexterior3_arrfuncs);
    pyexterior3_arrfuncs.getitem = npyexterior3_getitem;
    pyexterior3_arrfuncs.setitem = npyexterior3_setitem;
    pyexterior3_arrfuncs.copyswapn = pyexterior3_copyswapn;
    pyexterior3_arrfuncs.copyswap = pyexterior3_copyswap;
    pyexterior3_arrfuncs.nonzero = pyexterior3_nonzero;
    pyexterior3_descr_proto.typeobj = (PyTypeObject *)PyExterior3_Type;
    Py_SET_TYPE((PyObject *)&pyexterior3_descr_proto, &PyArrayDescr_Type);
    npy_exterior3 = PyArray_RegisterDataType(&pyexterior3_descr_proto);
    if (npy_exterior3 < 0) {
        goto fail;
    }
    descr = PyArray_DescrFromType(npy_exterior3);
    if (descr == NULL) {
        goto fail;
    }
    if (PyObject_SetAttrString(PyExterior3_Type, "dtype",
                               (PyObject *)descr) < 0) {
        Py_DECREF(descr);
        goto fail;
    }
    Py_DECREF(descr);
    if (algpy_set_dimension(PyExterior3_Type, 8) < 0) {
        goto fail;
    }
    if (pyexterior3_register_casts() < 0 ||
            pyexterior3_register_ufuncs(numpy) < 0 ||
            pyexterior3_register_lattice_ufuncs(ufuncs) < 0 ||
            PyModule_AddObjectRef(m, "exterior3", PyExterior3_Type) < 0) {
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
