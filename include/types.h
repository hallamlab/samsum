#ifndef __RPKM_TYPE
#define __RPKM_TYPE
#include <map>
#include <vector>
#include <iostream>

using namespace std;

typedef struct _READ_DATUM {
    unsigned int start,
            end;
    bool multi;
    std::string name;
} READ_DATUM;

typedef vector<READ_DATUM> READ_DATA;

typedef struct _CONTIG {
    /*
     * L is the length of the contig
     * M is the list of all reads aligned to the contig with start, end, boolean for multiread and its name
     */
    unsigned long L;
    int hits;
    READ_DATA M;
} CONTIG;


typedef struct _MATCH {
    /*
      * parity is 1 if it is second in the pair
      * mapped is 1 if the read was not unmapped
      * orphan is 1 if the mate was not successfully aligned
      * multi is 1 if the read is not a primary alignment
      * chimeric is 1 if parts of the read aligned to different loci
      * singleton is 1 if the mate was not successfully aligned
     */
    std::string query, subject;
    unsigned int start, end;
    bool parity; 
    bool mapped;
    bool orphan;
//    bool multi;
    bool chimeric;
    bool singleton;
    float  w; //no idea what this does...
    _MATCH(): w(0) { } 
} MATCH;


template< typename A, typename B, typename C, typename D>
struct QUADRUPLE {
     A first;
     B second;
     C third;
     D fourth;
};


typedef struct _COVERAGE {
    float coverage;
    float numreads;
    unsigned int sequence_length, uncovered_length;
} COVERAGE;
#endif //__RPKM_TYPE
