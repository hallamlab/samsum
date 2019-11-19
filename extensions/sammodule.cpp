#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

/*
Some function nonspecific boilerplate follows this section
*/

struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int module_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int module_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}
#define INITERROR return NULL
// End of boilerplate


// Function signatures go here
static PyObject *get_mapped_reads(PyObject *self, PyObject *args);

static PyObject *get_alignment_strings(PyObject *self, PyObject *args);
// End function signatures


// Docstrings for functions go here
static char get_mapped_reads_docstring[] =
        "Parses a SAM file and returns the read names of every read that was mapped to a reference sequence.\n";

static char get_alignment_strings_docstring[] =
        "Parses a SAM file and returns a string representing the first eight fields for every alignment made.\n";
// End of docstrings

// Define all of the module methods in this:
static PyMethodDef module_methods[] = {
        {"get_mapped_reads",
        get_mapped_reads,
        METH_VARARGS,
        get_mapped_reads_docstring},
        {NULL, NULL, 0, NULL},
        {"get_alignment_strings",
        get_alignment_strings,
        METH_VARARGS,
        get_alignment_strings_docstring},
        {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_sam_module",
    "Suite of functions for parsing SAM and BAM files.",
    sizeof(struct module_state),
    module_methods,    /* m_methods */
    NULL,              /* m_reload */
    module_traverse,   /* m_traverse */
    module_clear,      /* m_clear */
    NULL,              /* m_free */
};
//"get_mapped_reads",
//"Parses a SAM file and returns the read names of every read that was mapped to a reference sequence.",


PyMODINIT_FUNC PyInit__sam_module(void) {
    PyObject *m = PyModule_Create(&module_def);

    if (m == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(m);

    st->error = PyErr_NewException("_SAM.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

    return m;
}

static PyObject *get_mapped_reads(PyObject *self, PyObject *args) {
    char * aln_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sii", &aln_file, &min_length, &min_map_qual)) {
        return NULL;
    }

    PyObject *all_reads;
    /*
    *
    */
    std::cout << "Parsing alignment file " << aln_file << std::endl;
    return all_reads;
}

static PyObject *get_alignment_strings(PyObject *self, PyObject *args) {
    char * aln_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sii", &aln_file, &min_length, &min_map_qual)) {
        return NULL;
    }

    PyObject *all_reads;
    /*
    *
    */
    std::cout << "Parsing alignment file " << aln_file << std::endl;
    return all_reads;
}
