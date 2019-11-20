#ifndef _MATHOUTPUTPARSER
#define _MATHOUTPUTPARSER
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "utilities.h"
#include "types.h"

using namespace std;

class MatchOutputParser {
    protected:
        unsigned long num_unmapped_reads;
    public:
        /* class variables */
        std::string filename;
        std::string format;
        std::ifstream input;
        char buf[1000];
        vector<char *> tempv;
        vector<char *> fields;
        // The following variables are general file parsing stats
        unsigned long num_alignments;
        unsigned long num_mapped;
        unsigned long num_unmapped;
        unsigned long num_fwd;
        unsigned long num_rev;
        /* class functions */
        MatchOutputParser(const std::string &filename, const std::string &format);
        virtual ~MatchOutputParser() = 0;
        unsigned long get_Num_Unmapped_Reads();
        virtual bool nextline(MATCH &match)=0;
};

//subclass of the MatchOutputParser
class SamFileParser: virtual public MatchOutputParser {
    public:
        SamFileParser(const std::string &filename, const std::string &format);
        void consume_sam(vector<MATCH> &all_reads, map<std::string, float > &multireads);
        virtual bool nextline(MATCH &match);
        bool getMateInfo(unsigned int i, MATCH &match);
        ~SamFileParser();
};

#endif //_MATHOUTPUTPARSER
