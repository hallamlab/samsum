#include <map>
#include <vector>
#include <algorithm>
#include "types.h"
#include "utilities.h"
#include "helper.h"
#include <assert.h>

using namespace std;
bool compare_triplets(const READ_DATUM &a, const READ_DATUM &b) {
   return a.start < b.start ? true : false; 
}
 
int main( int argc, char **argv ){
    // parse options
    Options options;

    if (!options.SetOptions(argc, argv)) {
        options.print_usage(argv[0]);
        exit(0);
    }
    map<string, CONTIG> contigs_dictionary;
    unsigned long total_contig_length = create_contigs_dictionary(options.contigs_file, contigs_dictionary);

    std::ostream *output;
    std::ofstream rpkm_output;
    rpkm_output.open(options.output_file.c_str(), std::ifstream::out);
    output = &rpkm_output;

    std::ostream *stats_out;
    std::ofstream Sout;
    if (!options.stats_file.empty()) {
        Sout.open(options.stats_file.c_str(), std::ifstream::out);
        stats_out = &Sout;
    }
 
    vector<MATCH> all_reads;
    all_reads.reserve(80000000);
    
    // Creating the read multiplicity counts data structure to track reads mapping to multiple reference sequences
    map<std::string, float > multireads;

    RUN_STATS stats;
    RUN_STATS _stats;
    for( vector<string>::iterator it = options.read_map_files.begin(); it != options.read_map_files.end(); it++)  {
        if( it->empty() )
            continue;
        all_reads.clear();
        if( options.show_status ) 
            std::cout << "\n" << "Parsing alignments from SAM file " << *it << std::endl;
        _stats = consume_sam(*it, options.reads_map_file_format, all_reads, multireads, options.show_status);

        stats = stats + _stats;
        process_SAM(*it, contigs_dictionary, options.reads_map_file_format, all_reads, multireads, options.show_status);
        stats.num_multireads += (int) multireads.size();

        if( options.show_status )
            _stats.printStats(&std::cout);

        if(!options.stats_file.empty()) {
            *stats_out << "\nStats for file :  " << *it << std::endl;
            _stats.printStats(stats_out);
        }
    }

    if ( options.show_status )
        std::cout << "\nSorting the read matches... ";
    for (map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
        std::sort(it->second.M.begin(), it->second.M.end(), compare_triplets);

    }
    if (options.show_status)
        std::cout << "done." << endl;

    // Calculate the number of reads present in the dataset (not number of lines in the SAM file!)
    if (options.num_reads == 0)
        stats.num_distinct_reads = stats.num_distinct_reads_mapped + stats.num_distinct_reads_unmapped;
    else
        stats.num_distinct_reads = options.num_reads;

    float total_covered_length = 0.0;
    float rpkm_sum = 0.0;
    COVERAGE coverage;
    if (options.show_status)
        std::cout << "Computing coverage and RPKM values for all contigs... ";

    // Write header to output csv
    *output << "Sample_name,Sequence_name,Reads_mapped,RPKM" << endl;
    // Write the number of unmapped reads to the output file - necessary for calculating proportion of reads mapped
    *output << options.output_file << ",UNMAPPED," << stats.num_distinct_reads_unmapped << ",NA" << endl;

    for (map<string, CONTIG>::iterator it = contigs_dictionary.begin(); it != contigs_dictionary.end(); it++) {
        substring_coverage(contigs_dictionary, it->first, 1, it->second.L, coverage, 0, multireads, options.multi_reads);
        total_covered_length += coverage.coverage*contigs_dictionary[it->first].L;
        it->second.rpkm = (1E9/ static_cast<float>(stats.num_distinct_reads))*
                              (coverage.numreads/ static_cast<float>(coverage.sequence_length));
        it->second.hits = static_cast<int>(coverage.numreads);
        rpkm_sum += it->second.rpkm;
        // Write the sample-name, seq_name, number of hits, and RPKM
        *output << options.output_file << ',' << it->first << ',' << it->second.hits <<  ',' << it->second.rpkm << endl;
    }
    rpkm_output.close();

    if (options.show_status)
        std::cout << "done." << endl;

}