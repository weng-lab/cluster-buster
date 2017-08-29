// read a DNA sequence in FASTA format

#include <istream>
#include "MCFbio.hpp"

// Industrial strength fasta-format sequence reader
// appends the sequence to anything already in seq
std::istream & mcf::get_fasta(
  std::istream & strm,
  std::vector<unsigned> & seq,
  std::string & title,
  unsigned (*translator)(char))  // not allowed to repeat the default argument?
{
  std::string t;
  unsigned old_size = seq.size();
  char c;
  bool titflag = false;  // have we read a title line yet?

  while (strm >> c) {
    if (isalpha(c)) {
      seq.push_back(translator(c));
    } else if (c == '>') {
      if (titflag || seq.size() != old_size) {
	strm.unget();
	break;
      } else {
	std::getline(strm, t);
	titflag = true;
      }
    } else if (c == '#') {  // skip comments
      std::string junk;
      std::getline(strm, junk);
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof() && (titflag || seq.size() != old_size))
    strm.clear();

  if (strm)
    title = t;
  else
    seq.resize(old_size);
  return strm;
}
