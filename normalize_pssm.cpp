#include "MCFbio.hpp"
#include "MCFgen.hpp" // normalize
#include <cassert>
#include <numeric> // accumulate

void mcf::normalize_pssm(matrix<float> &pssm,
                         const std::vector<float> &pseudos) {
  const unsigned r = pssm.rows();
  const unsigned c = pssm.cols();
  assert(pseudos.size() == c);

  for (std::vector<float>::iterator k = pssm[0]; k < pssm[r]; k += c) {
    for (unsigned i = 0; i < c; ++i)
      k[i] += pseudos[i];
    normalize(k, k + c);
  }
}
