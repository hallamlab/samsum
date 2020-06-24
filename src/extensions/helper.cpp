#include "helper.h"
#include "types.h"
using namespace std;


// vector<std::string> format_matches_for_service(vector<MATCH *> &all_reads, char* &index) {
void format_matches_for_service(vector<MATCH *> &all_reads, char* &index) {
    /* Parameters:
      * all_reads: A vector populated by MATCH instances that contain information to be converted into strings
     * Functionality:
      * Iterates through all MATCH instances in all_reads and generates two strings for each: the name of the query
      and the relevant fields that are to be returned into Python.
    */


    for ( vector<MATCH *>::iterator it = all_reads.begin(); it != all_reads.end(); ++it)  {
        update_end_and_read_length(*it);
    }
}


void remove_low_quality_matches(vector<MATCH *> &mapped_reads, unsigned int min_map_qual, float &unmapped_weight_sum) {
    /* Parameters:
      * mapped_reads: A vector of MATCH instances that is to be filtered
      * min_map_qual: An integer representing the minimum mapping quality for an alignment to be included
      * unmapped_weight_sum: Reference to a float that tracks the sum weight of fragments
     * Functionality:
      * A new filtered_matches vector is created and all MATCH instances that pass the min_map_qual are appended to it.
      * Their weights are added to the unmapped_weight_sum float since these are no longer returned as MATCHes
    */
    vector<MATCH *> filtered_matches;
    filtered_matches.reserve(mapped_reads.size());
    for ( vector<MATCH *>::iterator it = mapped_reads.begin(); it != mapped_reads.end(); ++it)  {
        if ((*it)->mq < min_map_qual)
            unmapped_weight_sum += (*it)->w;
        else
            filtered_matches.push_back((*it));
    }
    mapped_reads.clear();
    mapped_reads = filtered_matches;
    filtered_matches.clear();
}
