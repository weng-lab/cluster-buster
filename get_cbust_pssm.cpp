// Read a position specific weight/count/score matrix in eFASTA format

#include <istream>
#include <sstream>
#include "get_cbust_pssm.hpp"

// Industrial strength PSSM reader
// (alphsize could be figured out from the number of entries per line...)
std::istream & cb::get_cbust_pssm(
  std::istream & strm,
  mcf::matrix<float> & mat,
  std::string & title,
  float & weight,
  float & gap,
  unsigned alphsize)  // not allowed to repeat the default argument?
{
  std::string t;
  float w = 1;  // default weight, if none is supplied
  float g = -1;  // dummy value for gap if none is supplied
  mcf::matrix<float> m(0, alphsize);
  char c = 0;
  bool titflag = false;  // have we read a title line yet?

  while (strm >> c) {
    if (c == '>') {
      if (titflag || m.rows() != 0) {
	strm.unget();
        break;
      } else {
	std::getline(strm, t);
	titflag = true;
	w = 1;  // reset weight to default
      }
    } else if (c == '#') {
      std::string temp;
      std::getline(strm, temp);
      std::istringstream is(temp);
      is >> temp;
      if (temp == "WEIGHT")
	is >> w;
      else if (temp == "GAP")
	is >> g;
    } else {
      strm.unget();
      std::vector<float> v;
      for (unsigned i = 0; i < alphsize; ++i) {
	float d;
	strm >> d;
	v.push_back(d);
      }
      if (!strm)
	return strm;  // failed to read alphsize doubles
      m.push_row(v.begin());
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof() && (titflag || m.rows() != 0))
    strm.clear();

  if (strm) {
    title = t;
    weight = w;
    gap = g;
    mat = m;
  }
  return strm;
}
