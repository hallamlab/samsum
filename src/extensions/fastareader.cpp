#include "fastareader.h"


FastaReader::FastaReader( const string & fasta_file) {
    /* Instantiate the FastaReader instance */
    map<string, unsigned long> seq_len_map;
    this->fasta_file = fasta_file;
    this->seq_lengths = seq_len_map;
    return;
}

string FastaReader::getContigsFileName() {
    return this->fasta_file;
}

//void FastaReader::get_sequence_lengths(map< string, unsigned long> &contigs_dictionary) {
void FastaReader::get_sequence_lengths() {

    std::ifstream input(this->fasta_file.c_str());
    if(!input.good()){
        std::cerr << "ERROR: Unable to open '"<<this->fasta_file << "'for reading." << std::endl;
        return ;
    }

    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                name = extract_sequence_name(name);
                this->seq_lengths[name] = content.size() ;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        name = extract_sequence_name(name);
//        std::cout << name << " : " << content.size() << std::endl;
        this->seq_lengths[name] = content.size() ;
    }
    input.close();
    return ;
}

//std::string FastaReader::extract_sequence_name(const std::string &name) {
//     char  cstr[1000000];
//     strcpy(cstr, name.c_str());
//
//     char * cptr = cstr;
//
//     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++;
//     (*cptr) ='\0';
//
//     std::string sname(cstr);
//     return sname;
//}


unsigned long create_contigs_dictionary(std::string fasta_file, std::map<std::string, CONTIG> &contigs_dictionary) {

     FastaReader fastareader(fasta_file);
     map<string, unsigned long> contig_lengths;
     map<string, unsigned long>::iterator it_contig_lens;

//     fastareader.get_sequence_lengths(contig_lengths);
     unsigned long genome_length = 0;
//     CONTIG contig;
//     for(it_contig_lens = contig_lengths.begin(); it_contig_lens != contig_lengths.end(); it_contig_lens++ ) {
//        genome_length += it_contig_lens->second;
//        contig.L = it_contig_lens->second;
//        contigs_dictionary[it_contig_lens->first] = contig;
//     }

     return genome_length;
}