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
                content.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
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