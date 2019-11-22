#include "sambamparser.h"

using namespace std;

MatchOutputParser::MatchOutputParser(const std::string &filename, const std::string &format) {
     this->filename = filename;
     this->format = format;
     this->num_unmapped = 0;
};

unsigned long MatchOutputParser::get_Num_Unmapped_Reads() {
    return this->num_unmapped;
}

MatchOutputParser::~MatchOutputParser() {
}

std::string MatchOutputParser::summarise() {
    char buf[1000];
    std::string summary_str;
    summary_str.assign("Summary for " + this->filename + ":\n");
    sprintf(buf, "\tNumber of alignments:           %ld\n", this->num_alignments);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of mapped reads:         %ld\n", this->num_mapped);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of unmapped reads:       %ld\n", this->num_unmapped);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of forward reads:        %ld\n", this->num_fwd);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of reverse reads:        %ld\n", this->num_rev);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of multireads:           %ld\n", this->num_multireads);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of singleton alignments: %ld\n", this->num_singletons);
    summary_str.append(buf);

    return summary_str;
}

SamFileParser::SamFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
    /*
      * Constructor for SamFileParser class
      * Attempts to open the SAM file that was provided as the file name and throws an error, and returns, if unable to
      * Set all SamFileParser variables used for counting alignments while parsing to 0
    */
     this->filename = filename;
     this->input.open(filename.c_str(), std::ifstream::in);

     if(!this->input.good()){
         std::cerr << "ERROR: Unable to open '"<< filename <<"' for reading." << std::endl;
         return ;
     }
     this->num_alignments = 0;
     this->num_mapped = 0;
     this->num_unmapped = 0;
     this->num_fwd = 0;
     this->num_rev = 0;
     this->num_multireads = 0;
     this->num_singletons = 0;
     this->header_pattern.assign('@', 1);
     this->unmapped_pattern.assign('*', 1);
}

SamFileParser::~SamFileParser() {
   this->input.close();
}

bool SamFileParser::getMateInfo(unsigned int i, MATCH &match)  {
    /*
      * Details for the match object are outlined in types.h
    */

    unsigned int a = i;
    bool singleton = 0;
    a = a >> 2;
    match.mapped = !(a&1); 
    singleton = a&1;

    a = a >> 1;
    match.orphan = a&1; 
    singleton = singleton^(a&1);

    a = a >> 3;
    if ( a&1 )  {
         match.parity = 0; 
         a = a >> 1;
    }
    else {
         a = a >> 1;
         if ( a&1 )
             match.parity  = 1;
         else
             return false;
    }

    a = a >> 4;
    match.chimeric = a&1; 
    match.singleton = singleton;
    return true;
}

bool SamFileParser::nextline(MATCH &match) {
    /*
      * Function for iterating through lines in a SAM file
      * match instance provided has been instantiated
      * while loop iterates until the SAM header is passed then reads a single line if the line has >= 9 fields
    */
     string line;

     bool _success = false;
     while (std::getline(this->input, line).good()) {
         if (matchString(line, this->header_pattern, true) )
             continue;

         fields.clear();
         split(line, fields, this->buf,'\t');

         if (fields.size() < 9) continue;

         _success = true;
         break;
     }
    
     if ( _success )  {
         match.query =  fields[0];
         match.subject = std::string(fields[2]);
         match.start = atoi(fields[3]);
         match.end =  match.start + std::string(fields[9]).size();
         getMateInfo(static_cast<unsigned int>(atoi(fields[1])), match);

         return true;
     }
    return false;
}

void SamFileParser::consume_sam(vector<MATCH> &all_reads,
                                map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict,
                                bool show_status) {
    /*
      * Basic function for parsing a SAM file
      * All mapped reads are saved as a MATCH instance and these objects are stored in the all_reads vector
      * The number of mapped, unmapped, forward, and reverse reads are counted.
      * These are counts are non-unique so double counts could arise from reads with multiple alignments
    */
    MATCH match;

    if ( show_status )
        std::cout << "Number of SAM lines processed: " << std::endl;

    int i;
    struct QUADRUPLE <bool, bool, unsigned int, unsigned int> p;
    for ( i =0; ; i++ ) {
        if (show_status && i % 10000 == 0)
            std::cout << "\n\033[F\033[J" << i;

        if (!this->nextline(match))
            break;

        if (match.mapped)
            this->num_mapped++;
        else
            this->num_unmapped++;

        if (match.parity)
            this->num_rev++;
        else this->num_fwd++;

        if (reads_dict.find(match.query) == reads_dict.end()) {
            p.first = false;
            p.second = false;
            p.third = 0;
            p.fourth = 0;
            reads_dict[match.query] = p;
        }
        this->num_alignments++;

        // if it is not mapped then ignore it
        if (!match.mapped)
            continue;

        if (match.parity) {
            reads_dict[match.query].first = true;  // This is a forward read
            reads_dict[match.query].third++;
        }
        else {
            reads_dict[match.query].second = true;  // This is a reverse read
            reads_dict[match.query].fourth++;
        }

        // store it to process later by looking up the dictionary
        try {
            all_reads.push_back(match);
        }
        catch (...) {
            cout << "Failing " << match.query << "   " << all_reads.size() << endl;
        }
    }

    if ( show_status )
        std::cout << "\n\033[F\033[J" << i << std::endl;

    cout << "Unique reads parsed: " << reads_dict.size() << endl;
    cout << "Number of mapped reads: " << all_reads.size() << endl;

    return;
}


int identify_multireads(map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict,
                         map<std::string, float > &multireads) {
    /* Parameters:
      * reads_dict: A map indexed by read-names with QUADRUPLE values that store all reads in the SAM file
      * multireads: An empty map that is passed by reference of strings indexing floats
     * Functionality:
      * Count the number of singleton reads (reads that didn't map), multireads, and secondary hits
      * Counting singletons: iterate through all read names (keys) in reads_dict
      and if the first and second elements of QUADRUPLE are false, the number of singletons is incremented.
      * Counting multireads: If the third or fourth variable of the QUADRUPLE is greater than one,
      this indicates the read was aligned multiple times so multireads is incremented by 1.
      * Counting secondary hits: The number of alignments of a read is tracked by the third and fourth elements
      of QUADRUPLE so the number of secondary hits is 1-(QUADRUPLE.third|QUADRUPLE.fourth)
    */
    int num_secondary_hits = 0;
    int num_singleton_reads = 0;
    int num_multireads = 0;

    for ( map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> >::iterator it = reads_dict.begin();
          it != reads_dict.end();
          it++) {
        if( !(it->second.first && it->second.second) )
            num_singleton_reads++;
        if( it->second.third > 1) {
            num_multireads++;
            multireads[it->first] = 0.0;
            num_secondary_hits += it->second.third-1;
        }
        if( it->second.fourth  > 1) {
            num_multireads++;
            multireads[it->first] = 0.0;
            num_secondary_hits += it->second.fourth-1;
        }
    }

    return num_secondary_hits;
}


void assign_read_weights(vector<MATCH> &all_reads,
                         map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict) {
    /* Parameters:
      * all_reads: A complete list of MATCH instances, one for each mapped read
      * reads_dict: A map indexed by read-names with QUADRUPLE values that store all reads in the SAM file
     * Functionality:
      * Basically calculates the weights such that the sum of the paired reads is 1 - for FPKM calculation
      * Iterate through the all_reads vector and depending on whether the read was forward (parity == 0) or reverse
      (parity == 1) the respective alignment count variable (third|fourth) is used to calculate the weight such that
      the sum of forward and reverse (if applicable) equals one.
    */

    int n = 0;

    for ( vector<MATCH>::iterator it = all_reads.begin(); it != all_reads.end(); it++)  {
        if ( it->parity == 0  ) {
            if( reads_dict[it->query].first && reads_dict[it->query].second )
                it->w = 0.5/static_cast<float>(reads_dict[it->query].third);
            else
                it->w = 1/static_cast<float>(reads_dict[it->query].third);
        }
        else  {
            if( reads_dict[it->query].first && reads_dict[it->query].second )
                it->w = 0.5/static_cast<float>(reads_dict[it->query].fourth);
            else
                it->w = 1/static_cast<float>(reads_dict[it->query].fourth);
        }
    }

    if (n == 0)
        std::cerr << "ERROR: alignments were parsed incorrectly (none found)" << std::endl;
    return;
}