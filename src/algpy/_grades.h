#ifndef ALGPY_GRADES_H
#define ALGPY_GRADES_H

/* Read-only grade access shared by the algebra and lattice wrappers.
 *
 * x.grade is a sequence view: x.grade[k] is the projection onto grade k, for
 * 0 <= k < len(x.grade), and negative indices count from the top grade.
 * Algebra scalars are immutable, so the view rejects item assignment.
 *
 * Each wrapper supplies a function returning the grade-k projection of its
 * scalar for 0 <= k < count. Each extension module creates its own view
 * type on first use. */

#include <Python.h>

typedef PyObject *(*algpy_grade_func)(PyObject *owner, Py_ssize_t k);

typedef struct {
    PyObject_HEAD
    PyObject *owner;
    algpy_grade_func grade;
    Py_ssize_t count;
} AlgpyGrades;

static PyObject *algpy_grades_type = NULL;

static Py_ssize_t
algpy_grades_length(PyObject *self)
{
    return ((AlgpyGrades *)self)->count;
}

static PyObject *
algpy_grades_item(PyObject *self, Py_ssize_t k)
{
    AlgpyGrades *view = (AlgpyGrades *)self;
    if (k < 0 || k >= view->count) {
        PyErr_SetString(PyExc_IndexError, "grade index out of range");
        return NULL;
    }
    return view->grade(view->owner, k);
}

static PyObject *
algpy_grades_subscript(PyObject *self, PyObject *key)
{
    Py_ssize_t k;
    if (!PyIndex_Check(key)) {
        PyErr_SetString(PyExc_TypeError, "grade indices must be integers");
        return NULL;
    }
    k = PyNumber_AsSsize_t(key, PyExc_IndexError);
    if (k == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (k < 0) {
        k += ((AlgpyGrades *)self)->count;
    }
    return algpy_grades_item(self, k);
}

static int
algpy_grades_ass_subscript(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(key),
                           PyObject *Py_UNUSED(value))
{
    PyErr_SetString(PyExc_TypeError,
                    "algebra values are immutable; grades cannot be assigned");
    return -1;
}

static PyObject *
algpy_grades_repr(PyObject *self)
{
    return PyUnicode_FromFormat("<grades of %R>", ((AlgpyGrades *)self)->owner);
}

static int
algpy_grades_traverse(PyObject *self, visitproc visit, void *arg)
{
    Py_VISIT(Py_TYPE(self));
    Py_VISIT(((AlgpyGrades *)self)->owner);
    return 0;
}

static int
algpy_grades_clear(PyObject *self)
{
    Py_CLEAR(((AlgpyGrades *)self)->owner);
    return 0;
}

static void
algpy_grades_dealloc(PyObject *self)
{
    PyTypeObject *tp = Py_TYPE(self);
    PyObject_GC_UnTrack(self);
    algpy_grades_clear(self);
    PyObject_GC_Del(self);
    Py_DECREF(tp);
}

static PyType_Slot algpy_grades_slots[] = {
    {Py_tp_dealloc, algpy_grades_dealloc},
    {Py_tp_traverse, algpy_grades_traverse},
    {Py_tp_clear, algpy_grades_clear},
    {Py_tp_repr, algpy_grades_repr},
    {Py_tp_doc, (void *)"Read-only grade projections: x.grade[k]."},
    {Py_sq_length, algpy_grades_length},
    {Py_sq_item, algpy_grades_item},
    {Py_mp_length, algpy_grades_length},
    {Py_mp_subscript, algpy_grades_subscript},
    {Py_mp_ass_subscript, algpy_grades_ass_subscript},
    {0, NULL}
};

static PyType_Spec algpy_grades_spec = {
    .name = "algpy.grades",
    .basicsize = sizeof(AlgpyGrades),
    .flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,
    .slots = algpy_grades_slots,
};

/* The x.grade view of owner, which has count grades. */
static inline PyObject *
algpy_grades_view(PyObject *owner, algpy_grade_func grade, Py_ssize_t count)
{
    AlgpyGrades *view;
    if (algpy_grades_type == NULL) {
        algpy_grades_type = PyType_FromSpec(&algpy_grades_spec);
        if (algpy_grades_type == NULL) {
            return NULL;
        }
    }
    view = PyObject_GC_New(AlgpyGrades, (PyTypeObject *)algpy_grades_type);
    if (view == NULL) {
        return NULL;
    }
    view->owner = Py_NewRef(owner);
    view->grade = grade;
    view->count = count;
    PyObject_GC_Track((PyObject *)view);
    return (PyObject *)view;
}

#endif /* ALGPY_GRADES_H */
