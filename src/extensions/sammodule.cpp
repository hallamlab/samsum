#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cstdlib>
#include <iostream>
#include "sambamparser.h"
#include <string.h>
//#include "helper.h"

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

    if(PyType_Ready(&MatchType) < 0)
        return NULL;
    
    Py_INCREF((PyObject *) &MatchType);
    PyModule_AddObject(m, "MATCH", (PyObject *) &MatchType);

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
    char * index;
    bool all_alignments;  // A flag indicating whether secondary and supplementary alignments should be used (True)
    int aln_percent;  // The minimum alignment length - this currently isn't used here
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sbiis", &aln_file, &all_alignments, &aln_percent, &min_map_qual, &index)) {
        return NULL;
    }

    PyObject *mapping_info_py = PyList_New(0);
    std::cout << "Parsing alignment file " << aln_file << std::endl;

    bool verbose = true;
    vector<MATCH *> mapped_reads;
    float unmapped_weight_sum;
    cout << "Reserving space for mapped reads... " << std::flush;
    mapped_reads.reserve(8000000); // still required if 
    cout << "done." << endl;
    map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > reads_dict;
    map<std::string, float > multireads;

    SamFileParser sam_file(aln_file, "sam");
    int status = sam_file.consume_sam(mapped_reads, all_alignments, verbose);
    if ( status > 0 )
        return mapping_info_py;

    if (check_reads_paired(mapped_reads))
        unmapped_weight_sum = (sam_file.num_unmapped*0.5);
    else
        unmapped_weight_sum = sam_file.num_unmapped;

    sam_file.alignment_multiplicity_audit(mapped_reads, reads_dict);

    // Identify multireads with and count the number of secondary adn supplementary alignments
    long num_secondary_hits = identify_multireads(reads_dict, multireads,
                                                  sam_file.num_multireads, sam_file.num_singletons);

    // Redistribute read weights using multiple alignment information in reads_dict
    assign_read_weights(mapped_reads, reads_dict);
    remove_low_quality_matches(mapped_reads, min_map_qual, unmapped_weight_sum);

    // Set the SamFileParser values
    sam_file.unique_queries = reads_dict.size();
    sam_file.secondary_alns = num_secondary_hits;
    sam_file.num_distinct_reads_mapped = sam_file.num_mapped - num_secondary_hits;

    // Add a match object that stores the number of unmapped reads
    
    MATCH *unmapped = Match_cnew();
    unmapped->w = unmapped_weight_sum;
    unmapped->query = "NA";
    unmapped->subject = "UNMAPPED";
    unmapped->parity = 0;
    mapped_reads.push_back(unmapped);

    // Print the various SAM alignment stats
    if ( verbose )
        std::cout << sam_file.summarise();

    if ( verbose )
        cout << "Calculating alignment positions... " << std::flush;

    add_alignment_positions(mapped_reads, index); //update match end and read_length

    if ( verbose )
        cout << "done." << endl << std::flush;

    if ( verbose )
        cout << "Building alignment list... " <<std::flush;

    long x = 0;
    vector<MATCH *>::iterator qi_it;
    for (qi_it = mapped_reads.begin(); qi_it != mapped_reads.end(); ++qi_it ) {
        MATCH *mt = (*qi_it);
        if (PyList_Append(mapping_info_py, Py_BuildValue("O", (PyObject *)mt)) == -1)
            x++;
    }

    if ( verbose )
        cout << "done." << endl << std::flush;

    if (x > 0) {
        sprintf(sam_file.buf, "WARNING: Failed to append %ld/%zu items into mapped reads list.", x, mapped_reads.size());
        cerr << sam_file.buf << endl;
    }
    return mapping_info_py;
}

static PyObject *get_alignment_strings(PyObject *self, PyObject *args) {
    /* Parameters:
      * args: A list of arguments received from the Python call that includes the SAM/BAM file, the minimum alignment
        length, and the minimum mapping quality for filtering.
     * Functionality:
      * IN PROGRESS: This will be developed if a need is identified.
      * Returns the alignment strings directly from the SAM/BAM file with minimal filtering of poor alignments.
    */
    char * aln_file;  // This could either be a SAM or BAM file
    int min_length;  // The minimum alignment length
    int min_map_qual;  // The minimum mapping quality
    if (!PyArg_ParseTuple(args, "sii", &aln_file, &min_length, &min_map_qual)) {
        return NULL;
    }

    PyObject *mapping_info_py = PyList_New(0);
    std::cout << "Parsing alignment file " << aln_file << std::endl;
    return mapping_info_py;
}
