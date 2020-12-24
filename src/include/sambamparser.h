#ifndef _MATHOUTPUTPARSER
#define _MATHOUTPUTPARSER
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "utilities.h"
#include "helper.h"
#include "types.h"

using namespace std;

class MatchOutputParser {
    protected:
        // The following variables are general file parsing stats
        unsigned long num_lines;
        unsigned long num_fwd;
        unsigned long num_rev;
        unsigned long num_unpaired;
    public:
        /* Class Variables */
        unsigned long unique_queries;
        unsigned long num_mapped;
        unsigned long num_unmapped;
        unsigned long num_multireads;
        unsigned long secondary_alns;
        unsigned long num_singletons;
        unsigned long num_distinct_reads_mapped;
        std::string filename;
        std::string format;
        std::ifstream input;
        char buf[1000];
        vector<char *> fields;
        /* Class Functions */
        MatchOutputParser(const std::string &filename, const std::string &format);
        virtual ~MatchOutputParser() = 0;
        std::string summarise();
        unsigned long get_Num_Unmapped_Reads();
        virtual bool nextline(MATCH *match)=0;
};

//subclass of the MatchOutputParser
class SamFileParser: virtual public MatchOutputParser {
    public:
        /* Class Variables */
        std::string header_pattern;
        std::string unmapped_pattern;
        /* Class Functions */
        SamFileParser(const std::string &filename, const std::string &format);
        int parse_header(map<std::string, int> &ref_dict);
        int consume_sam(vector<MATCH *> &all_reads, bool multireads, bool verbose);
        int alignment_multiplicity_audit(vector<MATCH *> &all_reads,
                                         map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict);
        virtual bool nextline(MATCH *match);
        bool getMateInfo(unsigned int i, MATCH *match);
        ~SamFileParser();
};

long identify_multireads(map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict,
                         map<std::string, float > &multireads, unsigned long &multi, unsigned long &num_singleton_reads);

float calculate_weight(int parity, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> &pair);

void assign_read_weights(vector<MATCH *> &all_reads,
                         map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict);

#endif //_MATHOUTPUTPARSER
