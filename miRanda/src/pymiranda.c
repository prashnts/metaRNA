#include "Python.h"

static PyObject* libpymiranda_find_targets(PyObject *self, PyObject *args, PyObject *keywds) {
    int voltage;
    char *state = "a stiff";
    char *action = "voom";
    char *type = "Norwegian Blue";

    static char *kwlist[] = {"voltage", "state", "action", "type", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "i|sss", kwlist,
                                     &voltage, &state, &action, &type))
        return NULL;

    printf("-- This parrot wouldn't %s if you put %i Volts through it.\n",
           action, voltage);
    printf("-- Lovely plumage, the %s -- It's %s!\n", type, state);

    Py_RETURN_NONE;
}

static PyMethodDef libpymiranda_methods[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and libpymiranda_parrot() takes
     * three.
     */
    {"parrot", (PyCFunction)libpymiranda_parrot, METH_VARARGS | METH_KEYWORDS,
     "Print a lovely skit to standard output."},
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef libpymirandamodule = {
    PyModuleDef_HEAD_INIT,
    "libpymiranda",
    NULL,
    -1,
    libpymiranda_methods
};

PyMODINIT_FUNC
PyInit_libpymiranda(void)
{
    return PyModule_Create(&libpymirandamodule);
}
