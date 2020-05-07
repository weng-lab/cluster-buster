#include <iosfwd>
#include <string>

namespace args {
using namespace std;
typedef unsigned uint;

extern string seqfile;      // Name of fasta format sequence(s) file
extern string matfile;      // Name of matrix file
extern double score_thresh; // Print clusters with scores >= this
extern double motif_thresh; // Print motifs inside clusters with scores >= this
extern double e_gap;        // Expected distance between pairs of cis-elements
extern bool gap_specified;  // Did the user specify a value for e_gap?
enum format { BY_SEQUENCE, BY_SEQUENCE_CONCISE, BY_SCORE, BY_SCORE_CONCISE, SEQUENCE_NAME_SORTED_BY_SCORE, BED };
extern format out_format;
extern uint bg_range;   // Go up to this far either side of current base
                        // when counting local base abundances
extern uint bg_padding; // Consider first X bp and last X bp of each sequence
                        // only for counting local base abundances but not for
                        // scoring to find motif clusters.
extern bool mask_lower; // mask lowercase letters?
extern double pseudo;   // Pseudocount to add to all matrix entries
extern uint keep_top_x_clusters_per_sequence; // Keep top X clusters per sequence
extern bool genomic_coordinates;              // Use genomic coordinates instead of relative coordinates
extern bool zero_based; // Is start coordinate in sequence name zero-based or one-based
extern double tau;      // Transition probability to "end" state of HMM
extern bool verbose;    // Be verbose or not

void parse(int argc, char **argv);
void print(ostream &strm, uint seq_num, uint mat_num);
}
