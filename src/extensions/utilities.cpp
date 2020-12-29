#include <stdlib.h>
#include "utilities.h"

/*
CML -- this function has been made robust to lines that are >1000 characters
by allocating more space for buf if required
*/
void split(const string &strn, std::vector<char *> &v, char *buf, char d) {
  if (strn.length() > 1000 )
      buf = (char *) malloc ((strn.length() + 1) * sizeof(char));
  strcpy(buf, strn.c_str());
  v.clear();
  v.reserve(15);
  v.push_back(buf);
  while(*buf != '\0') {
     if(*buf==d) {
       *buf = '\0';
       v.push_back(buf+1);
     }
     buf++;
  }
}


char* lstrip(char *str, char c) {
    /* utilities::lstrip:
      Parameters:
        str: A string instance to subset, based on the first character found
        c: A character used to subset the string by
      Functionality:
        Finds the position of 'c' in 'str' and removes all characters to the left of 'c' from str
        Returns 1 if the character 'c' wasn't found in str
    */
    while (*str != '\0') {
        if (*str == c) {
            str++;
            return str;
        }
        str++;
    }
    std::cerr << "ERROR: character '" << c << "' was not found in string." << std::endl;
    std::exit(1);
}


bool match_string(const string &str, const string & string_to_match, bool from_start) {
    unsigned long pos = str.find(string_to_match);
    if ( from_start && pos == 0 ) return true;

    return !from_start && pos != string_to_match.npos;

}
