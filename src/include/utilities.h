#ifndef _UTILITIES
#define _UTILITIES
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>

using namespace std;

void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');

bool match_string(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

char* lstrip(char *str, char c);

#endif //_UTILITIES

