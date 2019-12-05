#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cstdlib>
#include <iostream>
#include "sambamparser.h"
#include "helper.h"

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
    /*
      * Create a new SamFileParser instance
      * Read the alignments using SamFileParser::consume_sam()
      * Identify the reads with multiple alignments (mutlireads)
      * Redistribute the weights of these reads based on its alignment multiplicity
      * Return a list of interleaved `read_name`s and `reference_name, weight, left-most position, CIGAR`
    */
    char * aln_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sii", &aln_file, &min_length, &min_map_qual)) {
        return NULL;
    }

    PyObject *mapping_info_py = PyList_New(0);
    std::cout << "Parsing alignment file " << aln_file << std::endl;

    bool verbose = true;
    vector<MATCH> mapped_reads;
    mapped_reads.reserve(80000000);
    map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > reads_dict;
    map<std::string, float > multireads;

    SamFileParser sam_file(aln_file, "sam");
    int status = sam_file.consume_sam(mapped_reads, reads_dict, verbose);
    if ( status > 0 )
        return mapping_info_py;

    // Identify multireads with and count the number of secondary adn supplementary alignments
    long num_secondary_hits = identify_multireads(reads_dict, multireads,
                                                  sam_file.num_multireads, sam_file.num_singletons);

    // Redistribute read weights using multiple alignment information in reads_dict
    assign_read_weights(mapped_reads, reads_dict);

    // Set the SamFileParser values
    sam_file.unique_queries = reads_dict.size();
    sam_file.secondary_alns = num_secondary_hits;
    sam_file.num_distinct_reads_mapped = sam_file.num_mapped - num_secondary_hits;

    // Print the various SAM alignment stats
    if ( verbose )
        std::cout << sam_file.summarise();

    // Reformat the MATCH objects into the strings required
    if ( verbose )
        cout << "Formatting alignment strings... ";
    vector<std::string> query_info = format_matches_for_service(mapped_reads);
    if ( verbose )
        cout << "done." << endl;

    if ( verbose )
        cout << "Converting strings to Python objects... ";
    long x = 0;
    vector<std::string>::iterator qi_it;
    std::string str;
    for (qi_it = query_info.begin(); qi_it != query_info.end(); qi_it++ ) {
        str = *qi_it;
        if (PyList_Append(mapping_info_py, Py_BuildValue("s", str.c_str())) == -1)
            x++;
    }
    if ( verbose )
        cout << "done." << endl;
    if (x > 0) {
        sprintf(sam_file.buf, "WARNING: Failed to append %ld/%ld items into mapped reads list.", x, query_info.size());
        cerr << sam_file.buf << endl;
    }

    return mapping_info_py;
}

static PyObject *get_alignment_strings(PyObject *self, PyObject *args) {
    char * aln_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sii", &aln_file, &min_length, &min_map_qual)) {
        return NULL;
    }

    PyObject *mapping_info_py;
    std::cout << "Parsing alignment file " << aln_file << std::endl;
    return mapping_info_py;
}