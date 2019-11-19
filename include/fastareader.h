#ifndef __FASTAREADER
#define __FASTAREADER
#include "types.h"
#include "utilities.h"
#include <string.h>

using namespace std;


class FastaReader {

private:
     string contigs_file;
public:
     FastaReader(const string & contigs_file);
     void get_fasta_sequence_info(map<string, unsigned long> &contigs_dictionary) ;
     std::string extract_sequence_name(const std::string &name);
     string getContigsFileName();
};

unsigned long create_contigs_dictionary(std::string contigs_file, std::map<std::string, CONTIG> &contigs_dictionary);

#endif // __FASTAREADER
