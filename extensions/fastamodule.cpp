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
static PyObject *get_lengths(PyObject *self, PyObject *args);

static char get_lengths_docstring[] =
        "Reads the FASTA file, returning a list of interleaved record headers and sequences.";

static PyMethodDef module_methods[] = {
        {"get_lengths",
        get_lengths,
        METH_VARARGS,
        get_lengths_docstring},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_fasta_module",
    "This module rapidly loads and a FASTA file using C++. There is no position-level formatting or validation.",
    sizeof(struct module_state),
    module_methods,    /* m_methods */
    NULL,              /* m_reload */
    module_traverse,   /* m_traverse */
    module_clear,      /* m_clear */
    NULL,              /* m_free */
};

PyMODINIT_FUNC PyInit__fasta_module(void) {
    PyObject *m = PyModule_Create(&module_def);

    if (m == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(m);

    st->error = PyErr_NewException("_FASTA_module.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

    return m;
}


static PyObject *get_lengths(PyObject *self, PyObject *args) {
    char * fasta_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    if (!PyArg_ParseTuple(args, "si", &fasta_file, &min_length)) {
        return NULL;
    }

    PyObject *seq_lens;
    /*
    *
    */
    std::cout << "Parsing FASTA file " << fasta_file << std::endl;
    return seq_lens;
}