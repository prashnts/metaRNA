#include "Python.h"

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "miranda.h"

double scale = 4.0;
int strict = 0;
double gap_open = -9.0;
int temperature_override = 25;
double gap_extend = -4.0;
double score_threshold = 100.0;
double energy_threshold = 1.0;
int length_5p_for_weighting = 8;
int length_3p_for_weighting;

static PyObject* libpymiranda_find_targets(PyObject *self,
    PyObject *args, PyObject *keywds) {

  char* gene_seq;
  char* mirna_seq;

  initialize_bases();
  initialize_scores();

  static char *kwlist[] = {"gene_seq", "mirna_seq", "scale", "strict",
    "gap_open", "gap_extend", "score_threshold", "energy_threshold",
    "length_5p_for_weighting", "temperature", NULL};

  int parsed = PyArg_ParseTupleAndKeywords(args, keywds,
    "ss|diddddii",
    kwlist,
    &gene_seq, &mirna_seq, &scale, &strict, &gap_open, &gap_extend,
    &score_threshold, &energy_threshold, &length_5p_for_weighting,
    &temperature_override
  );

  ExpString *outjson = 0;
  create_ExpString(&outjson);
  find_targets(gene_seq, mirna_seq, outjson);

  PyObject* ret_val = Py_BuildValue("s", access_ExpString(outjson));

  destroy_ExpString(&outjson);

  if (parsed) {
    return ret_val;
  } else {
    return NULL;
  }
}

static PyMethodDef libpymiranda_methods[] = {
  {
    "find_targets",
    (PyCFunction)libpymiranda_find_targets,
    METH_VARARGS | METH_KEYWORDS,
    "Find the miRNA and Gene Sequence targets and Binding Energy."
  },
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef libpymirandamodule = {
  PyModuleDef_HEAD_INIT,
  "libpymiranda",
  NULL,
  -1,
  libpymiranda_methods
};

PyMODINIT_FUNC PyInit_libpymiranda(void) {
  return PyModule_Create(&libpymirandamodule);
}
