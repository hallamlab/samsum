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

vector<std::string> format_matches_for_service(vector<MATCH> &all_reads);

void remove_low_quality_matches(vector<MATCH> &mapped_reads, unsigned int min_map_qual, float &unmapped_weight_sum);

#endif //_HELPER
