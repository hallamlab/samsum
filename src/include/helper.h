#ifndef _HELPER
#define _HELPER
#include <map>
#include <string>
#include <ostream>
#include <iterator>
#include <assert.h>
#include "types.h"
#include "sambamparser.h"

using namespace std;

unsigned long create_contigs_dictionary(std::string contigs_file, std::map<std::string, CONTIG> &contigs_dictionary);

void process_SAM(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,
                 const std::string &reads_map_file_format,
                 vector<MATCH> &all_reads,
                 map<std::string, float > &multireads,
                 bool show_status= false);

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,
                        unsigned long start, unsigned long end,
                        COVERAGE &coverage, unsigned int maxReadLength,
                        map<std::string, float > &multireads, bool multi_reads);

vector<std::string> format_matches_for_service(vector<MATCH> &all_reads);

void remove_low_quality_matches(vector<MATCH> &mapped_reads, unsigned int min_map_qual, float &unmapped_weight_sum);

#endif //_HELPER
