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

void add_alignment_positions(vector<MATCH *> &all_reads, char* &index);
void remove_low_quality_matches(vector<MATCH *> &mapped_reads, unsigned int min_map_qual, float &unmapped_weight_sum);
bool check_reads_paired(vector<MATCH *> &mapped_reads);

#endif //_HELPER
