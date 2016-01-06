#include <iosfwd>  // "forward declaration" of istream
#include <string>
#include "matrix.hpp"

namespace cb {
  // Like get_weighted_pssm, but if a line beginning #GAP is present,
  // reads a numeric gap value from this line.
  std::istream & get_cbust_pssm(
    std::istream & strm,     // stream to read from
    mcf::matrix<float> & mat,  // store the matrix here
    std::string & title,     // store title line minus leading '>'
    float & weight,          // = 1 by default if not specified
    float & gap,             // = -1 by default if not specified
    unsigned alphsize = 4);  // perhaps it should guess alphsize from the input
}
