#include <iosfwd>
#include <string>

namespace args {
  using namespace std;
  typedef unsigned uint;

  extern string seqfile;  // Name of fasta format sequence(s) file
  extern string matfile;  // Name of matrix file
  extern double score_thresh;  // Print clusters with scores >= this
  extern double motif_thresh;  // Print motifs inside clusters with scores >= this
  extern double e_gap;  // Expected distance between pairs of cis-elements
  extern bool gap_specified;  // Did the user specify a value for e_gap?
  enum format { BY_SEQUENCE, BY_SEQUENCE_CONCISE, BY_SCORE, BY_SCORE_CONCISE };
  extern format out_format;
  extern uint bg_range;  // Go up to this far either side of current base
                         // when counting local base abundances
  extern bool mask_lower;  // mask lowercase letters?
  extern double pseudo;  // Pseudocount to add to all matrix entries
  extern double tau;  // Transition probability to "end" state of HMM

  void parse(int argc, char ** argv);
  void print(ostream & strm, uint seq_num, uint mat_num);
}
