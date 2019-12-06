#include "helper.h"
using namespace std;


#define _MAX 1000000000

void process_SAM(const std::string & reads_map_file, std::map<string, CONTIG> &contigs_dictionary,
                 const std::string &reads_map_file_format,
                 std::vector<MATCH> &all_reads,
                 map<std::string, float > &multireads,
                 bool show_status) {
    /* Parameters:
      * multireads: A map indexed by read names and is otherwise empty
     * Functionality:
      *
     */
    MATCH match;

    READ_DATA read_data;
    READ_DATUM read_datum;
    
    int i =0;
    
    if( show_status ) std::cout << "Number of hits processed : " ;
    // iterate through individual hits/alignments
    for(vector<MATCH>::iterator it=all_reads.begin();  it!= all_reads.end(); it++ )  {

        if( i >=_MAX ) break;

        if (show_status && i%10000==0) {
           std::cout << "\n\033[F\033[J";
           std::cout << i ;
        }
        i++;

        if ( contigs_dictionary.find(it->subject)==contigs_dictionary.end() ) {
            std::cout << "Missing contig " << it->subject << std::endl;
            std::cerr << "ERROR : Could not find the matched contig in the contig file " << std::endl;
            exit(1);
        }

        // Update the READ_DATUM instance with mutli (multiread) boolean value
        read_datum.name = it->query;
        if ( multireads.find(it->query) != multireads.end() ) {
            multireads[it->query] += 1;
            read_datum.multi = true;
        }
        else
            read_datum.multi = false;

        // Where is this data coming from?
        if (  it->start < it->end ) {
            read_datum.start = it->start;
            read_datum.end = it->end;
        }
        else {
            read_datum.start = it->end;
            read_datum.end = it->start;
        }

        contigs_dictionary[it->subject].M.push_back(read_datum);
    }
}

unsigned int getMaxReadSize( std::map<string, vector<MATCH> > &orf_dictionary,
                             std::map<string, CONTIG> &contigs_dictionary) {
    unsigned int size = 0;
    std::map<string, vector<MATCH> >::iterator itcont;

    for(itcont= orf_dictionary.begin(); itcont != orf_dictionary.end(); itcont++)  {
       for(std::vector<READ_DATUM>::iterator it = contigs_dictionary[itcont->first].M.begin(); it != contigs_dictionary[itcont->first].M.end(); it++) {
          if( size < it->end  - it->start ) size = it->end  - it->start ;
       }
    }

    return size;
} 


std::vector<READ_DATUM>::iterator binary_search(std::vector<READ_DATUM> &A, int seekValue) {
  // continually narrow search until just one element remains
  unsigned long imin, imax;
  imin = 0;
  imax = A.size();

  while (imin < imax)
    {
      unsigned int imid = (imin+imax)/2;
 
      // code must guarantee the interval is reduced at each iteration
      assert(imid < imax);

      // note: 0 <= imin < imax implies imid will always be less than imax
 
      // reduce the search
      if (A[imid].start < static_cast<unsigned int>(seekValue) )
        imin = imid + 1;
      else
        imax = imid;
    }

    std::vector<READ_DATUM>::iterator it = A.begin() + imin;
    return it ;
}

void substring_coverage(std::map<string, CONTIG> &contigs_dictionary, const std::string &contig,
                        unsigned long start, unsigned long end,
                        COVERAGE &coverage, unsigned int maxReadLength,
                        map<std::string, float > &multireads, bool multi_reads) {

    if ( contigs_dictionary.find(contig) == contigs_dictionary.end() || contigs_dictionary[contig].L == 0 ) {
        coverage.coverage = 0 ;
        coverage.numreads = 0 ;
        coverage.sequence_length = end - start ;
        coverage.uncovered_length = 0;
    }

    float numreads = 0;
    float _coverage = 0;
    unsigned long uncovered_length = 0;
    unsigned long p_end = start;

    int _seekValue =  maxReadLength == 0 ? 0 :  (start < maxReadLength ? 0 : start-maxReadLength );

    std::vector<READ_DATUM>::iterator it= contigs_dictionary[contig].M.begin();

    if ( _seekValue >0 )
        it =  binary_search(contigs_dictionary[contig].M, _seekValue);

    // iterate through every read that aligned to that contig
    for ( ; it != contigs_dictionary[contig].M.end(); it++) {

        uncovered_length  +=  ( p_end > it->start  || it->start > end) ? 0 : it->start - p_end;
        //make sure the read start and end are not going beyond the contig
        if( it->end > p_end )
            p_end = it->end;

        if( (start <= it->start && it->start <= end) ||  (start <= it->end && it->end <= end)  ) {
            numreads += 1;
        }

        // the subsequent reads are going off the end of the orf
        if( it->start > end )
            break;

        if (multi_reads && it->multi) {
            float read_multiplicity = multireads.find(it->name)->second;
            numreads += 1.0/read_multiplicity;
        }
    }


    uncovered_length += (p_end > end ) ? 0 : (end - p_end);

    unsigned long sequence_length = end - start;
    if( sequence_length > 0 )
        _coverage = ((float)(sequence_length - uncovered_length )/(float)sequence_length)*100;

    coverage.numreads = numreads;
    coverage.coverage = _coverage;
    coverage.sequence_length  = sequence_length;
    coverage.uncovered_length =  uncovered_length;
}

vector<std::string> format_matches_for_service(vector<MATCH> &all_reads) {
    /* Parameters:
      * all_reads: A vector populated by MATCH instances that contain information to be converted into strings
     * Functionality:
      * Iterates through all MATCH instances in all_reads and generates two strings for each: the name of the query
      and the relevant fields that are to be returned into Python.
    */
    vector<std::string> query_info;
    char buf[1000];

    for ( vector<MATCH>::iterator it = all_reads.begin(); it != all_reads.end(); it++)  {
        query_info.push_back(it->query);
        sprintf(buf, "%s\t%d\t%s\t%f", it->subject.c_str(), it->start, it->cigar.c_str(), it->w);
        query_info.push_back(buf);
    }

    return query_info;
}
