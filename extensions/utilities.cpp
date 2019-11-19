#include <stdlib.h>
#include "utilities.h"

std::string extract_sequence_name(const std::string &name) {
     char  cstr[1000000];
     strcpy(cstr, name.c_str());
     
     char * cptr = cstr;

     while( *cptr != '\t' && *cptr !=  ' ' && *cptr != '\0' )  cptr++; 
     (*cptr) ='\0';

     std::string sname(cstr);
     return sname;
}

/*
CML -- this function has been made robust to lines that are >1000 characters
by allocating more space for buf if required
*/
void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  if (strn.length() > 1000 )
      buf = (char *) malloc ((strn.length() + 1) * sizeof(char));
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.reserve(15);
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {

    unsigned long pos = str.find(stringtomatch);
    if(fromstart && pos ==0 ) return true;

    return !fromstart && pos >= 0;

}
