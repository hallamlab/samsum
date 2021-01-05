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
    sprintf(buf, "\tNumber of alignment lines:      %ld\n", this->num_lines);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of aligned reads:        %ld\n", this->num_mapped);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of unmapped reads:       %ld\n", this->num_unmapped);
    summary_str.append(buf);
    sprintf(buf, "\tUnique queries aligned          %ld\n", this->unique_queries);
    summary_str.append(buf);
    sprintf(buf, "\tForward reads aligned:          %ld\n", this->num_fwd);
    summary_str.append(buf);
    sprintf(buf, "\tReverse reads aligned:          %ld\n", this->num_rev);
    summary_str.append(buf);
    sprintf(buf, "\tUnpaired read alignments:       %ld\n", this->num_unpaired);
    summary_str.append(buf);
    sprintf(buf, "\tNumber of multireads:           %ld\n", this->num_multireads);
    summary_str.append(buf);
    sprintf(buf, "\tSecondary alignments:           %ld\n", this->secondary_alns);
    summary_str.append(buf);
    sprintf(buf, "\tOrphan alignments:              %ld\n", this->num_singletons);
    summary_str.append(buf);

    return summary_str;
}

SamFileParser::SamFileParser(const std::string &filename, const std::string &format):MatchOutputParser(filename, format) {
    /* Parameters:
      * filename: Name of the SAM file to be parsed
     * Functionality:
      * Constructor for SamFileParser class
      * Attempts to open the SAM file that was provided as the file name and throws an error, and returns, if unable to
      * Sets all SamFileParser variables used for counting alignments while parsing to 0
    */
     this->filename = filename;
     this->input.open(filename.c_str(), std::ifstream::in);
     this->num_lines = 0;
     this->unique_queries = 0;
     this->num_mapped = 0;
     this->num_unmapped = 0;
     this->num_fwd = 0;
     this->num_rev = 0;
     this->num_unpaired = 0;
     this->num_multireads = 0;
     this->secondary_alns = 0;
     this->num_singletons = 0;
     this->num_distinct_reads_mapped = 0;
     this->header_pattern.assign("@", 1);
     this->unmapped_pattern.assign("*", 1);
     return;
}

SamFileParser::~SamFileParser() {
   this->input.close();
}

bool SamFileParser::getMateInfo(unsigned int bitflag, MATCH *match)  {
    /* Parameters:
      * bitflag: The second column in a SAM file with bitwise encodings of mapping information
      * match: A MATCH instance
     * Functionality:
      * Perform bit-wise calculations to get information on the read's alignment, mate pairing, etc.
      * Returns `true` if either the first or second read in pair mapped, else `false`.
      * Details for the MATCH object are in types.h
    */

    unsigned int a = bitflag;
    bool orphan = 0;
    bool non_primary = 0;  // This is True if neither of the non-primary or supplementary alignment bits are set
    a = a >> 2;  // Skip the "read mapped" and "mapped in proper pair" bits
    match->mapped = !(a&1);  // mapped == 1 if the read was mapped b/c the third (0x4) bit isn't set
    orphan = a&1;  // orphan == 0 if the read was unmapped

    a = a >> 1;  // Move to the next, "mate unmapped" bit
    orphan = orphan^(a&1);  // `orphan` is 1 if the mate was unmapped (or doesn't exist) XOR `orphan` is set to 1

    a = a >> 3;  // Move to the sixth (0x64) "first in pair" bit
    if ( a&1 )  {
        match->parity = false;
        a = a >> 1;
    }
    else {
        a = a >> 1;  // Move to the seventh (0x128) "second in pair" bit
        if ( a&1 )
            match->parity  = true;
        else
            return false;
    }

    a = a >> 1;  // Move to the eighth, "not primary alignment" position
    non_primary = a&1;  // Should be 0 if it is the primary alignment
    a = a >> 3;  // Move to the eleventh (0x2048) "supplementary alignment" bit position
    match->multi = non_primary^(a&1);  // Is a multiread if it is either a non-primary XOR secondary alignment
    match->chimeric = a&1;  // This hints at a possible chimera, but it would have to be validated downstream
    match->singleton = orphan;
    return true;
}

bool SamFileParser::nextline(MATCH *match) {
    /* Parameters:
      * match: Reference to a MATCH instance that is to be populated with alignment information
     * Functionality:
      * Function for adding values of a short-read alignment in a SAM file to a match instance.
      * Specifically, the `paired`, `query`, `subject`, `start`, and `end` values are populated.
      * If the line of the SamFileParser matches the header_pattern, lines are skipped until they no longer match
      and the line has >= 9 tab-separated columns.
    */
    if (this->fields.size() < 9) return false;

    match->query = (char *)malloc(strlen(this->fields[0]) + 1);
    strcpy(match->query, this->fields[0]);
    match->subject = (char *)malloc(strlen(this->fields[2]) + 1);
    strcpy(match->subject, this->fields[2]);
    match->start = atoi(this->fields[3]);
    match->mq = atoi(this->fields[4]);
    match->cigar = (char *)malloc(strlen(this->fields[5]) + 1);
    strcpy(match->cigar, this->fields[5]);
    match->paired = getMateInfo(static_cast<unsigned int>(atoi(this->fields[1])), match);

    // TODO: test to ensure the end position is calculated correctly. It currently isn't.
    if ( match->parity ) // Read is second in pair and will be aligned right-to-left
    match->end =  match->start - std::string(this->fields[9]).size();
    else
    match->end =  match->start + std::string(this->fields[9]).size();

    return true;
}

int SamFileParser::parse_header(map<std::string, int> &ref_dict) {
    /* Parameters:
      * ref_dict: Pointer to a map of strings (to be reference names) as values and integers as keys
     * Functionality:
      * Iterates over the lines in a SAM file (SamFileParser.input attribute) while the lines match the
       SamFileParser.header_pattern attribute ('@').
      * Returns the line number that the header ends at.
    */
    string line;
    int line_no = 0;
    while (std::getline(this->input, line).good()) {
        if (match_string(line, this->header_pattern, true) ) {
            this->fields.clear();
            split(line, this->fields, this->buf, '\t');
            if (strcmp(this->fields[0], "@SQ") == 0) {
                ref_dict[lstrip(this->fields[1], ':')] = atoi(lstrip(this->fields[2], ':'));
            }
            else
                continue;
        }
        else {
            this->input.seekg(long(this->input.tellg())-(line.size()+1));
            return line_no;
        }
        line_no++;
    }
    return line_no;
}

int SamFileParser::consume_sam(vector<MATCH *> &all_alignments, bool multireads, bool show_status) {
    /* Parameters:
      * all_alignments: Pointer to a vector of MATCH objects that has yet to be populated
      * multireads: Boolean flag indicating whether reads that have multiple ambiguous mapping positions are used
      * show_stats: Boolean indicating whether the number of reads parsed should be printed to screen
     * Functionality:
      * Basic function for parsing a SAM file.
      * All mapped reads are saved as a MATCH instance and these objects are stored in all_alignments.
      * The number of mapped, unmapped, forward, and reverse reads are counted.
      * These are counts are non-unique so double counts could arise from reads with multiple alignments
    */
    string line;
    map<std::string, int> ref_dict;

     if(!this->input.good()) {
         std::cerr << "ERROR: Unable to open '"<< filename <<"' for reading." << std::endl;
         return 1;
     }

    this->parse_header(ref_dict);

    if ( show_status )
        std::cout << "Number of SAM alignment lines processed: " << std::endl;

    while (std::getline(this->input, line).good()) {
        this->num_lines++;
        if (show_status && this->num_lines % 10000 == 0)
            std::cout << "\n\033[F\033[J" << this->num_lines;
        this->fields.clear();
        split(line, this->fields, this->buf, '\t');
        if ( match_string(string(this->fields[2]), this->unmapped_pattern, true) ) {
            this->num_unmapped++;
            continue;
        }

        MATCH *match = Match_cnew();
        if (!this->nextline(match))
            break;

        this->num_mapped++;

        if (!match->paired)
            this->num_unpaired++;
        else {
            if (match->parity)
                this->num_rev++;
            else this->num_fwd++;
        }

        if (match->multi && !multireads)  // Drop secondary and supplementary alignments
            continue;

        // if it is not mapped then ignore it
        if (!match->mapped) {
            Py_DECREF((PyObject*)match);
            continue;
        }

        // store it to process later by looking up the dictionary
        try {
            all_alignments.push_back(match);
        }
        catch (...) {
            PyErr_Format(PyExc_RuntimeError, "Failing at %s.", match->query);
            return 1;
        }
    }
    this->fields.clear();

    if ( show_status )
        std::cout << "\n\033[F\033[J" << this->num_lines << std::endl;

    return 0;
}


int SamFileParser::alignment_multiplicity_audit(vector<MATCH *> &all_alignments,
                                                map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict) {
    /* Parameters:
      * all_alignments: Pointer to a vector of MATCH objects that has yet to be populated
      * reads_dict: Pointer to a map indexed by read-names with QUADRUPLE values that store all reads in the SAM file
     * Functionality:
      * 
    */
    struct QUADRUPLE <bool, bool, unsigned int, unsigned int> p;
    for ( vector<MATCH *>::iterator it = all_alignments.begin(); it != all_alignments.end(); ++it)  {
        if (reads_dict.find((*it)->query) == reads_dict.end()) {
            p.first = false;
            p.second = false;
            p.third = 0;
            p.fourth = 0;
            reads_dict[(*it)->query] = p;
        }

        if (!(*it)->parity) {
            reads_dict[(*it)->query].first = true;  // This is a forward read
            if ((*it)->mapped)
                reads_dict[(*it)->query].third++;
        }
        else {
            reads_dict[(*it)->query].second = true;  // This is a reverse read
            if ((*it)->mapped)
                reads_dict[(*it)->query].fourth++;
        }                       
    }
    return 0;
}


long identify_multireads(map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict,
                         map<std::string, float > &multireads, unsigned long &multi, unsigned long &num_singletons) {
    /* Parameters:
      * reads_dict: A map indexed by read-names with QUADRUPLE values that store all reads in the SAM file
      * multireads: An empty map that is passed by reference of strings indexing floats
     * Functionality:
      * Count the number of orphan reads (reads with mates that didn't map), multireads, and secondary hits
      * Counting singletons: iterate through all read names (keys) in reads_dict
      and if either the first or second elements of QUADRUPLE are false, the number of singletons is incremented.
      * Counting multireads: If the third or fourth variable of the QUADRUPLE is greater than one,
      this indicates the read was aligned multiple times so multireads is incremented by 1.
      * Counting secondary hits: The number of alignments of a read is tracked by the third and fourth elements
      of QUADRUPLE so the number of secondary hits is 1-(QUADRUPLE.third|QUADRUPLE.fourth)
    */
    long num_secondary_hits = 0;

    for ( map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> >::iterator it = reads_dict.begin();
          it != reads_dict.end();
          ++it) {
        if( !(it->second.first && it->second.second) )
            num_singletons++;
        if( it->second.third > 1) {
            multi++;
            multireads[it->first] = 0.0;
            num_secondary_hits += it->second.third-1;
        }
        if( it->second.fourth  > 1) {
            multi++;
            multireads[it->first] = 0.0;
            num_secondary_hits += it->second.fourth-1;
        }
    }

    return num_secondary_hits;
}


float calculate_weight(int parity, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> &pair) {
    float numerator = 1.0;
    // Is the read from a paired-end library AND did both of the reads map?
    if ( pair.first && pair.second )
        numerator = 0.5;
    // Calculate the read's weight based on the number of times it aligned
    if ( parity )  // This read is reverse
        return numerator/static_cast<float>(pair.fourth);
    else  // This read is forward
        return numerator/static_cast<float>(pair.third);
}


void assign_read_weights(vector<MATCH* > &all_reads,
                         map<std::string, struct QUADRUPLE<bool, bool, unsigned int, unsigned int> > &reads_dict) {
    /* Parameters:
      * all_reads: A complete list of MATCH instances, one for each mapped read
      * reads_dict: A map indexed by read-names with QUADRUPLE values that store all reads in the SAM file
     * Functionality:
      * Basically calculates the weights such that the sum of the paired reads is 1 - for FPKM calculation
      * Iterate through the all_reads vector and depending on whether the read was forward (parity == false) or reverse
      (parity == true) the respective alignment count variable (third|fourth) is used to calculate the weight such that
      the sum of forward and reverse (if applicable) equals one.
    */
    int n = 0;
    for ( vector<MATCH *>::iterator it = all_reads.begin(); it != all_reads.end(); ++it)  {
        (*it)->w = calculate_weight((*it)->parity, reads_dict[(*it)->query]);
        n++;
    }

    if (n == 0)
        PyErr_SetString(PyExc_TypeError, "Alignments were parsed incorrectly (none found)");
    return;
}
