#include "Python.h"

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "miranda.h"

double scale = 4.0;
int strict = 0;
int debug = 0;
double gap_open = -9.0;
double gap_extend = -4.0;
double score_threshold = 100.0;
double energy_threshold = 1.0;
int length_5p_for_weighting = 8;
int length_3p_for_weighting;
int key_value_pairs = 1;
int no_energy = 0;
int verbosity = 0;
int truncated = 0;
int restricted = 0;

static PyObject* libpymiranda_find_targets(PyObject *self, PyObject *args, PyObject *keywds) {

  char* gene_seq = "GCTACAGTTTTTATTTAGCATGGGGATTGCAGAGTGACCAGCACACTGGACTCC"
    "GAGGTGGTTCAGACAAGACAGAGGGGAGCAGTGGCCATCATCCTCCCGCCAGGAGCTTCTTCGTTCCTG"
    "CGCATATAGACTGTACATTATGAAGAATACCCAGGAAGACTTTGTGACTGTCACTTGCTGCTTTTTCTG"
    "CGCTTCAGTAACAAGTGTTGGCAAACGAGACTTTCTCCTGGCCCCTGCCTGCTGGAGATCAGCATGCCT"
    "GTCCTTTCAGTCTGATCCATCCATCTCTCTCTTGCCTGAGGGGAAAGAGAGATGGGCCAGGCAGAGAAC"
    "AGAACTGGAGGCAGTCCATCTA";

  char* mirna_seq = "UAGCAGCACGUAAAUAUUGGCG";

  initialize_bases(); /* Prepare the generic base lookup array*/
  initialize_scores();
    static char *kwlist[] = {"gene_seq", "mirna_seq", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|ss", kwlist,
                                     &gene_seq, &mirna_seq))
        return NULL;

    char* obj = find_targets(gene_seq, mirna_seq);

    return Py_BuildValue("s", obj);
}

static PyMethodDef libpymiranda_methods[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and libpymiranda_find_targets() takes
     * three.
     */
    {"find_targets", (PyCFunction)libpymiranda_find_targets, METH_VARARGS | METH_KEYWORDS,
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
