#ifndef __FASTAREADER
#define __FASTAREADER
#include "types.h"
#include "utilities.h"
#include <string.h>

using namespace std;


class FastaReader {
    private:
        string fasta_file;
    public:
        /* class variables */
        int num_records;
        map<string, unsigned long> seq_lengths;
        /* class functions */
        FastaReader(const string & fasta_file);
        void get_sequence_lengths() ;
//        std::string extract_sequence_name(const std::string &name);
        string getContigsFileName();
};

unsigned long create_contigs_dictionary(std::string fasta_file, std::map<std::string, CONTIG> &contigs_dictionary);

#endif // __FASTAREADER
