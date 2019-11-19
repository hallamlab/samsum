#ifndef __FASTAREADER
#define __FASTAREADER
#include "utilities.h"
#include <string.h>
#include <Python.h>

using namespace std;

static PyObject *get_lengths(PyObject *self, PyObject *args);

static char get_lengths_docstring[] =
        "Reads the FASTA file and formats it (checking duplicate headers, ambiguity characters, etc.) for TreeSAPP";

static PyMethodDef module_methods[] = {
        {"_get_lengths",
        get_lengths,
        METH_VARARGS,
        get_lengths_docstring},
        {NULL, NULL, 0, NULL}
};

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

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_get_lengths",
    "This module rapidly loads and formats sequences from a FASTA file using C from within TreeSAPP",
    sizeof(struct module_state),
    module_methods,    /* m_methods */
    NULL,              /* m_reload */
    module_traverse,   /* m_traverse */
    module_clear,      /* m_clear */
    NULL,              /* m_free */
};

#define INITERROR return NULL

PyMODINIT_FUNC PyInit__fasta_reader(void) {
    PyObject *m = PyModule_Create(&module_def);

    if (m == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(m);

    st->error = PyErr_NewException("_get_lengths.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

    return m;
}

class FastaReader {

private:
     string contigs_file;
public:
     FastaReader(const string & contigs_file);
     void get_fasta_sequence_info(map<string, unsigned long> &contigs_dictionary) ;
     std::string extract_sequence_name(const std::string &name);
     string getContigsFileName();
};

#endif // __FASTAREADER
