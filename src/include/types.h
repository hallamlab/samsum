#ifndef __RPKM_TYPE
#define __RPKM_TYPE
#include <Python.h>
#include <map>
#include <vector>
#include <iostream>
#include <ctype.h>
#include <stdlib.h>
#include "structmember.h"
using namespace std;


typedef struct {
    /*
      * paired is true if the read has a mate pair, else false (from a single-end run)
      * parity is 0 if it is the forward read, 1 if it is reverse in the pair
      * mapped is 1 if the read was not unmapped
      * orphan is 1 if the mate was not successfully aligned
      * multi is 1 if the read is not a primary alignment
      * chimeric is 1 if parts of the read aligned to different loci
      * singleton is 1 if the mate was not successfully aligned
     */
    PyObject_HEAD
    char * query;
    char *subject;
    char *cigar;
    /*unsigned int start, end, mq; */
    unsigned int start, end, mq, read_length;
    bool paired;
    bool parity; // Forward or reverse
    bool mapped; // Did it map to a reference sequence
    bool orphan;
    bool multi;
    bool chimeric;  // Whether part of the read aligned to multiple different reference sequences
    bool singleton; // Whether its mate was aligned or not
    float w, percent_id; // The weight of that read, based on the number of alignments
    //_MATCH(): w(0) { } 
} MATCH;

extern PyTypeObject MatchType;


MATCH *Match_cnew(PyTypeObject *type = &MatchType);

void update_end_and_read_length(MATCH * self);

template< typename A, typename B, typename C, typename D>
struct QUADRUPLE {
     A first;
     B second;
     C third;
     D fourth;
};

#endif //__RPKM_TYPE
