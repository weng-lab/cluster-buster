// Cluster-buster: find cis-element clusters in DNA sequences.
// Descendent of Cister and Comet.

// Method: employs a simple hidden Markov model of cis-element clusters.
// Any subsequence receives a score of:
// ln [ prob( subsequence | cluster model ) / prob( subsequence | background
// model ) ].
// Uses a heuristic to find optimal non-overlapping subsequences:
// First applies the forward algorithm to the entire sequence,
// and records all segments that exhibit maximal rises in the forward score.
// Then applies the backward algorithm starting from each segment endpoint,
// up to slightly before the segment start point.

#include <algorithm>
#include <cassert>
#include <cctype> // toupper
#include <cmath>  // pow, exp, log
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

#include "MCFbio.hpp"
#include "MCFgen.hpp" // die
#include "args.hpp"
#include "get_cbust_pssm.hpp"
#include "matrix.hpp"
#include "remove_overlapping_segments.hpp"

// reflect the matrix through a horizontal line:
template <class T> inline void reverse_matrix(mcf::matrix<T> &m) {
  const unsigned num_rows = m.rows();

  for (unsigned r = 0; r < num_rows / 2;
       ++r) // integer division: rounds fractions down
    std::swap_ranges(m[r], m[r + 1], m[num_rows - r - 1]);
}

// DNA-to-number translator that only recognizes uppercase
inline unsigned nolower_translator(char c) {
  switch (c) {
  case 'A':
    return 0u;
  case 'C':
    return 1u;
  case 'G':
    return 2u;
  case 'T':
    return 3u;
  default:
    return 4u;
  }
}

inline char number_to_revcomp_DNA(unsigned b) {
  static const char lookup[] = "tgcan";
  assert(b < 5);
  return lookup[b];
}

template <class T> struct byscore { // for sorting
  bool operator()(const T &a, const T &b) { return a.score > b.score; }
};

namespace cb {
using namespace std;
using mcf::matrix;

typedef unsigned uint;
typedef matrix<double> mat; // change to float??
typedef pair<uint, uint> segment;

struct s_segment {
  uint start;
  uint end;
  double score;
  s_segment(uint s, uint e, double sc) : start(s), end(e), score(sc) {}
  s_segment() {}
};

struct seq_info {
  string name;
  uint length;
  string chrom;
  uint genomic_pos;
  string extra_info;
  seq_info(const string &n, uint l) : name(n), length(l), chrom(n),
                                      genomic_pos(0), extra_info("") {}
  seq_info(const string &n, uint l, const string &c, uint p,
           const string &e) :
      name(n), length(l), chrom(c), genomic_pos(p), extra_info(e) {}
};

struct motif { // a predicted cis-element within a cluster
  string name;
  uint start, end; // zero-based coordinates, relative to whole sequence
  char strand;
  double score;
  string motif_seq;
  motif(uint mat_index, uint pos, uint clust_pos, double sc);
};

struct result { // a cis-element cluster
  uint seq_num;
  uint start, end; // zero-based coordinates
  double score;
  vector<double> motif_scores;
  string cluster_seq;
  vector<motif> hits;
  result(uint sn, uint s, uint e, double sc, const vector<double> &m,
         const vector<motif> &h);
};

const uint alphsize = 4;

vector<uint> seq;
vector<mat> raw_mats, mats;
vector<string> mat_names;
double gap_score;
uint max_motif_width;
string::size_type max_motif_name_len; // for pretty-printing
vector<result> results;

cb::seq_info get_chrom_and_pos(string &seq_name, uint length, bool zero_based);
void init_bg(uint *bg_counts, uint &bg_tot);
void get_bg(vector<double> &bg);
void get_hits(uint start, uint end, const vector<double> &bg,
              vector<motif> &hits);
void forward(uint start, uint end, const vector<double> &bg,
             vector<double> &scores, vector<pair<uint, uint> > &clusters);
uint backward(uint start, uint end, const vector<double> &bg,
              vector<double> &scores, uint ignore);
void scan_seq(uint seq_num);
void misc_init();
void fifth_column(const matrix<float> &m1, matrix<double> &m2);
void get_matrices();
void print_hits(ostream &strm, const seq_info &seq, const vector<motif> &hits);
void output_by_seq(ostream &strm, const seq_info &seq);
void output_by_seq_concise(ostream &strm, const seq_info &seq);
void output_by_seq_bed(ostream &strm, const seq_info &seq);
void output_by_score(ostream &strm, const vector<seq_info> &seqs);
void output_by_score_concise(ostream &strm, const vector<seq_info> &seqs);
void output_sequence_name_sorted_by_score(ostream &strm, const vector<seq_info> &seqs);
}

inline cb::seq_info cb::get_chrom_and_pos(string &seq_name, uint length, bool zero_based) {
  // Extract chromosome name, start position and extra info from
  // sequence name.
  //
  // Examples of sequence names from which chromosome names and
  // start positions can be extracted:
  //   - chr10:123456
  //   - chr10:123456-234567
  //   - chr10:123456@@gene_name
  //   - chr10:123456-234567@@gene_name

  string::size_type pos_extra_info = seq_name.find("@@");
  string chrom_and_pos;
  string chrom;
  uint genomic_pos;
  string extra_info;

  // Check if we have extra info (text after "@@").
  if (pos_extra_info != string::npos && pos_extra_info > 2) {
    chrom_and_pos = seq_name.substr(0, pos_extra_info);
    extra_info = seq_name.substr(pos_extra_info + 2);
  } else {
    chrom_and_pos = seq_name;
    extra_info = "";
  }

  string::size_type pos_chrom_pos = chrom_and_pos.find(":");

  // Check if we have a chromosome name before the ":".
  if (pos_chrom_pos != string::npos && pos_chrom_pos > 0) {
    chrom = chrom_and_pos.substr(0, pos_chrom_pos);

    // Extract the genomic start position from strings with:
    //   - only start position: "([0-9]+)"
    //   - start and end position: "([0-9]+)-[0-9]+"
    //   - no number ==> set to 0
    genomic_pos = atoi(chrom_and_pos.substr(pos_chrom_pos + 1).c_str());

    // Convert one-based to zero-based start position, if necessary.
    if (!zero_based) {
      genomic_pos -= 1;
    }
  } else {
    // No chromosomal position was found.
    chrom = seq_name;
    genomic_pos = 0;
  }

  return cb::seq_info(seq_name, length, chrom, genomic_pos, extra_info);
}

cb::motif::motif(uint mat_index, uint s, uint e, double sc)
    : name(mat_names[mat_index / 2]), start(s), end(e),
      strand(mat_index % 2 ? '-' : '+'),
      score(sc) { // Some old C++ systems don't think string has push_back():
  if (strand == '+')
    for (uint i = s; i <= e; ++i)
      motif_seq += mcf::number_to_DNA(seq[i]);
  else
    for (uint i = e; i != s - 1; --i) // deal with uint wrapping
      motif_seq += number_to_revcomp_DNA(seq[i]);
}

cb::result::result(uint sn, uint s, uint e, double sc, const vector<double> &m,
                   const vector<motif> &h)
    : seq_num(sn), start(s), end(e), score(sc), motif_scores(m), hits(h) {
  // don't need the body of this function for concise_output
  assert(e >= s);
  for (uint i = s; i <= e; ++i)
    cluster_seq += mcf::number_to_DNA(seq[i]);

  // put predicted cis-elements in uppercase:
  for (vector<motif>::const_iterator h = hits.begin(); h != hits.end(); ++h) {
    assert(h->start >= start && h->end >= h->start);
    uint x = h->start - start;
    uint y = h->end - start;
    for (string::iterator i = cluster_seq.begin() + x;
         i <= cluster_seq.begin() + y; ++i)
      *i = toupper(*i);
  }
}

// assumes bg_counts and bg_tot have been initialized
void cb::init_bg(uint *bg_counts, uint &bg_tot) {
  using args::bg_range;

  const vector<uint>::const_iterator end =
      seq.size() > bg_range ? seq.begin() + bg_range : seq.end();

  for (vector<uint>::const_iterator n = seq.begin(); n != end; ++n)
    if (*n < alphsize)
      ++bg_counts[*n], ++bg_tot;
}

// calculates background prob of each base by counting base freqs locally
// 'n' bases get garbage but finite values
// assumes bg is empty
void cb::get_bg(vector<double> &bg) {
  using args::bg_range;

  mcf::reserve_or_die(bg, seq.size());

  uint bg_counts[alphsize + 1] = {1, 1, 1, 1, 1}; // Laplace's rule
  uint bg_tot = alphsize;

  init_bg(bg_counts, bg_tot);

  for (vector<uint>::const_iterator n = seq.begin(); n != seq.end(); ++n) {
    uint i;

    if (uint(seq.end() - n) > bg_range && (i = *(n + bg_range)) < alphsize)
      ++bg_counts[i], ++bg_tot;

    bg.push_back(log(double(bg_tot) / bg_counts[*n]));
    //    cerr << mcf::number_to_DNA(*n) << " " << bg_counts[*n] /
    //    double(bg_tot) << endl;

    if (uint(n - seq.begin()) >= bg_range && (i = *(n - bg_range)) < alphsize)
      --bg_counts[i], --bg_tot;
  }
}

// Find matrix hits within a range of the sequence (standard PSSM search)
// Uses 'raw' matrices without absorbed transition probabilities
void cb::get_hits(uint start, uint end, const vector<double> &bg,
                  vector<motif> &hits) {
  for (uint n = start; n <= end; ++n) {
    for (vector<mat>::const_iterator m = raw_mats.begin(); m != raw_mats.end();
         ++m) {
      if (end - n + 1 < m->rows())
        continue;
      double s = 0;
      for (uint k = 0; k < m->rows(); ++k)
        s += (*m)[k][seq[n + k]] + bg[n + k];
      if (s > args::motif_thresh)
        hits.push_back(motif(m - raw_mats.begin(), n, n + m->rows() - 1, s));
    }
  }
}

// apply HMM forward algorithm to seq
// get all segments with maximal score increases: store them in "segs"
void cb::forward(uint start, uint end, const vector<double> &bg,
                 vector<double> &scores, vector<segment> &segs) {
  assert(start <= end && end < seq.size());
  assert(bg.size() == seq.size() && scores.size() == seq.size());

  uint lo = 0; // tracks position with lowest score so far

  for (uint n = 0; n < start; ++n) {
    scores.push_back(0.0);
  }

  for (uint n = start; n != end; ++n) {
    double score = gap_score;
    if (n > 0)
      score += scores[n - 1];

    for (vector<mat>::const_iterator m = mats.begin(); m != mats.end(); ++m) {
      double s;
      if (n >= m->rows())
        s = scores[n - m->rows()];
      else if (n + 1 == m->rows())
        s = 0;
      else
        continue;
      for (uint k = 0; k < m->rows(); ++k)
        s += (*m)[k][seq[n - k]] + bg[n - k]; // matrices are backwards
      score += log1p(exp(s - score)); // about 15% faster than log(1+...)
    }

    //    cerr << n << "  " << score << endl;
    scores.push_back(score);

    if (score <= scores[lo])
      lo = n;
    else if (score > scores[n - 1]) {
      uint seg_start = n - 1;
      for (uint i = segs.size() - 1; i != ~0u; --i) { // Never do ~0 without u
        const segment &c = segs[i];
        if (scores[c.second] >= score)
          break;
        if (scores[c.first] < scores[seg_start]) {
          seg_start = c.first;
          segs.resize(i);
          if (seg_start == lo)
            break;
        }
      }
      //      cerr << "#   " << seg_start << "   " << n << endl;
      segs.push_back(make_pair(seg_start, n));
    }
  }

  for (uint n = end; n < seq.size(); ++n) {
    scores.push_back(0.0);
  }
}

// apply backward algorithm from end to start, returning position of max score
unsigned cb::backward(uint start, uint end, const vector<double> &bg,
                      vector<double> &scores, uint ignore) {
  assert(start <= end && end < seq.size());
  assert(bg.size() == seq.size() && scores.size() == seq.size());

  uint hi = end;
  --start; // awful kludge to deal with uint wrapping

  for (uint n = end; n != start; --n) {
    double score = gap_score;
    if (n != end)
      score += scores[n + 1];

    for (vector<mat>::const_iterator m = mats.begin(); m != mats.end(); ++m) {
      if (uint(m - mats.begin()) / 2 == ignore)
        continue;
      double s;
      if (end - n >= m->rows())
        s = scores[n + m->rows()];
      else if (end - n + 1 == m->rows())
        s = 0;
      else
        continue;
      for (uint k = 0; k < m->rows(); ++k)
        s += (*m)[k][seq[n + k]] + bg[n + k];
      score += log1p(exp(s - score)); // about 15% faster than log(1+...)
    }

    //    cerr << n << "  " << score << endl;
    scores[n] = score;
    if (score > scores[hi])
      hi = n;
  }

  return hi;
}

// The main algorithm: find cis-element clusters in 1 sequence
void cb::scan_seq(uint seq_num) {
  vector<double> bg; // Holds background nucleotide probabilites
  get_bg(bg);

  vector<double> scores;
  mcf::reserve_or_die(scores, seq.size()); // could avoid this memory usage?

  vector<segment> segs;
  vector<s_segment> s_segs;

  if (args::bg_padding * 2 + 1 > seq.size()) {
    mcf::die("Sequence should be at least "
        + to_string(args::bg_padding * 2 + 1) + " bp long.");
  }

  for (uint i = 0; i < mats.size(); ++i)
    reverse_matrix(mats[i]);

  forward(args::bg_padding, seq.size() - args::bg_padding, bg, scores, segs);

  for (uint i = 0; i < mats.size(); ++i)
    reverse_matrix(mats[i]);

  for (vector<segment>::iterator s = segs.begin(); s != segs.end(); ++s) {
    //    cerr << s->first << "  " << s->second << endl;
    uint start =
        s->first + 2 > args::bg_padding + max_motif_width ? s->first + 2 - max_motif_width : args::bg_padding;
    start = backward(start, s->second, bg, scores, mat_names.size());
    if (scores[start] > args::score_thresh)
      s_segs.push_back(s_segment(start, s->second, scores[start]));
  }

  // Sort by cluster score.
  sort(s_segs.begin(), s_segs.end(), byscore<s_segment>());

  // Remove overlapping segements (returned segments are sorted by position).
  mcf::remove_overlapping_segments(s_segs, s_segs); // overkill

  if (args::keep_top_x_clusters_per_sequence > 0
        && args::keep_top_x_clusters_per_sequence < s_segs.size()) {
      // Sort by cluster score again.
      sort(s_segs.begin(), s_segs.end(), byscore<s_segment>());
      // Keep only the top X clusters per sequence.
      s_segs.resize(args::keep_top_x_clusters_per_sequence);
  }

  for (vector<s_segment>::const_iterator s = s_segs.begin(); s != s_segs.end();
       ++s) {
    vector<double> motif_scores; // each motif's score contribution
    for (uint m = 0; m != mat_names.size(); ++m) {
      backward(s->start, s->end, bg, scores, m);
      motif_scores.push_back(s->score - scores[s->start]);
      // motif_scores.push_back(s->score + log1p(-exp(scores[s->start] -
      // s->score)));
    }
    vector<motif> hits;
    if (args::out_format == args::BY_SEQUENCE ||
        args::out_format == args::BY_SCORE ||
        args::out_format == args::BED) {
      get_hits(s->start, s->end, bg, hits);
    }

    results.push_back(
        result(seq_num, s->start, s->end, s->score, motif_scores, hits));
  }
}

// Initialize miscellaneous variables: need to read matrices first!
void cb::misc_init() {
  gap_score = log(args::e_gap / (args::e_gap + 1.0));

  max_motif_width = 0;
  for (vector<cb::mat>::const_iterator m = cb::mats.begin();
       m != cb::mats.end(); ++m)
    max_motif_width = max(max_motif_width, m->rows());

  max_motif_name_len = 0;
  for (vector<string>::const_iterator i = mat_names.begin();
       i != mat_names.end(); ++i)
    max_motif_name_len = max(max_motif_name_len, i->size());
}

// add a column to a matrix, with all cells = minus_infinity
void cb::fifth_column(const matrix<float> &m1, matrix<double> &m2) {
  for (uint r = 0; r < m1.rows(); ++r) {
    copy(m1[r], m1[r + 1], m2[r]);
    m2[r][alphsize] = mcf::minus_infinity;
  }
}

// read matrices from a file and preprocess them
// Preprocessing steps:
// 1. normalize rows (counts -> probs)
// 2. take logs of all entries
// 3. add an extra column of all minus_infinities (for 'n' bases)
// 4. make reverse complements
// 5. absorb HMM transition probabilities
void cb::get_matrices() {
  ifstream file(args::matfile.c_str());
  if (!file)
    mcf::die("Sorry, couldn't open " + args::matfile);

  const vector<float> pseudos(alphsize, args::pseudo);
  matrix<float> matf;
  string title;
  vector<float> weights;
  float tot_weight = 0;
  float weight;
  float gap = -1;

  while (get_cbust_pssm(file, matf, title, weight, gap, alphsize)) {
    istringstream is(title);
    is >> title; // get first word (?)
    if (matf.rows() == 0)
      mcf::die("Empty matrix not allowed: " + title);
    mat_names.push_back(title);

    mcf::normalize_pssm(matf, pseudos);
    mcf::log_pssm(matf);
    matrix<double> matd(matf.rows(), alphsize + 1);
    fifth_column(matf, matd);
    raw_mats.push_back(matd);
    matf.rotate180(); // reverse complement
    fifth_column(matf, matd);
    raw_mats.push_back(matd);

    weights.insert(weights.end(), 2, weight);
    tot_weight += weight * 2;
    if (!args::gap_specified && gap >= 0)
      args::e_gap = gap;
  }

  if (!file.eof()) // catches some but not all errors
    mcf::die("Sorry, couldn't understand the matrix file " + args::matfile);
  if (raw_mats.size() == 0)
    mcf::die("No matrices read");
  if (tot_weight == 0)
    mcf::die("All matrix weights = 0: not allowed");

  mats = raw_mats;
  double continue_logprob = log(1.0 - args::tau);
  // HMM transition probability into cis-elements:
  double cis_start_logprob = continue_logprob - log(args::e_gap + 1.0);

  for (uint m = 0; m < mats.size(); ++m) {
    double w = weights[m] / tot_weight;
    w = (w == 0 ? mcf::minus_infinity : log(w));
    double term = cis_start_logprob + w - mats[m].rows() * continue_logprob;
    /*    cout << term << "  " << continue_logprob << " Matrix:\n";
    for (uint r = 0; r < mats[m].rows(); ++r) {
      for (uint c = 0; c < mats[m].cols(); ++c)
        cout << mats[m][r][c] << "\t";
      cout << "\n";
      }*/
    transform(mats[m][0], mats[m][1], mats[m][0],
              bind2nd(plus<double>(), term));
  }
}

void cb::print_hits(ostream &strm, const seq_info &seq,
                    const vector<motif> &hits) {
  for (vector<motif>::const_iterator h = hits.begin(); h != hits.end(); ++h)
    strm << h->name << "\t" << seq.genomic_pos + h->start + 1 << "\t"
         << seq.genomic_pos + h->end + 1 << "\t" << h->strand << "\t"
         << h->score << "\t" << h->motif_seq << '\n';
}

void cb::output_by_seq(ostream &strm, const seq_info &seq) {
  strm << '>' << seq.name << " (" << seq.length << " bp)\n\n";

  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    strm << "CLUSTER " << r - results.begin() + 1 << '\n'
         << "Location: " << seq.genomic_pos + r->start + 1 << " to "
                         << seq.genomic_pos + r->end + 1 << '\n'
         << "Score: " << r->score << '\n';
    vector<pair<double, string> > x;
    for (uint m = 0; m != mat_names.size(); ++m)
      x.push_back(make_pair(r->motif_scores[m], mat_names[m]));
    sort(x.begin(), x.end(), greater<pair<double, string> >());
    for (vector<pair<double, string> >::const_iterator i = x.begin();
         i != x.end(); ++i)
      strm << i->second << ": " << i->first << '\n';
    strm << r->cluster_seq << '\n';
    print_hits(strm, seq, r->hits);
    strm << '\n';
  }
}

void cb::output_by_seq_concise(ostream &strm, const seq_info &seq) {
  strm << '>' << seq.name << " (" << seq.length << " bp)\n";
  strm << "# Score\tStart\tEnd";
  for (vector<string>::const_iterator m = mat_names.begin();
       m != mat_names.end(); ++m)
    strm << "\t" << *m;
  strm << '\n';

  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    strm << r->score << "\t" << seq.genomic_pos + r->start + 1 << "\t"
         << seq.genomic_pos + r->end + 1;
    for (vector<double>::const_iterator m = r->motif_scores.begin();
         m != r->motif_scores.end(); ++m)
      strm << "\t" << *m;
    strm << '\n';
  }

  strm << '\n';
}

void cb::output_by_seq_bed(ostream &strm, const seq_info &seq) {
  static bool printed_header = false;

  if (!printed_header) {
    strm << "# chrom\t"
            "genomic_start__bed\t"
            "genomic_end__bed\t"
            "cluster_id_or_motif_name\t"
            "cluster_or_motif_score\t"
            "strand\t"
            "seq_name\t"
            "relative_start__bed\t"
            "relative_end__bed\t"
            "seq_number\t"
            "cluster_or_motif\t"
            "cluster_id\t"
            "motif_id\t"
            "motif_sequence\t"
            "motif_type_contribution_score\t"
            "extra_info\n";

    printed_header = true;
  }

  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    uint cluster_number = r - results.begin() + 1;
    string cluster_id = seq.name + "__cluster_" + to_string(cluster_number);

    strm << seq.chrom << "\t"
         << seq.genomic_pos + r->start << "\t"
         << seq.genomic_pos + r->end + 1 << "\t"
         << cluster_id << "\t"
         << r->score << "\t"
         << "+\t"
         << seq.name << "\t"
         << r->start << "\t"
         << r->end + 1 << "\t"
         << r->seq_num + 1 << "\t"
         << "cluster\t"
         << cluster_id << "\t"
         << "-\t"
         << "-\t"
         << "-\t"
         << ((seq.extra_info != "") ? seq.extra_info : "-") << "\n";

    for (vector<motif>::const_iterator h = r->hits.begin();
         h != r->hits.end(); ++h) {
      // Get motif type contribution score.
      vector<string>::iterator mat_names_iterator = find(mat_names.begin(),
                                                         mat_names.end(),
                                                         h->name);
      int matrix_name_index = distance(mat_names.begin(), mat_names_iterator);

      strm << seq.chrom << "\t"
           << seq.genomic_pos + h->start << "\t"
           << seq.genomic_pos + h->end + 1 << "\t"
           << h->name << "\t"
           << h->score << "\t"
           << h->strand << "\t"
           << seq.name << "\t"
           << h->start << "\t"
           << h->end + 1 << "\t"
           << r->seq_num + 1 << "\t"
           << "motif\t"
           << cluster_id << "\t"
           << cluster_id << "__motif_" << h->name << "\t"
           << h->motif_seq << "\t"
           << r->motif_scores[matrix_name_index] << "\t"
           << ((seq.extra_info != "") ? seq.extra_info : "-") << "\n";
    }
  }
}

void cb::output_by_score(ostream &strm, const vector<seq_info> &seqs) {
  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    strm << "CLUSTER " << r - results.begin() + 1 << '\n'
         << '>' << seqs[r->seq_num].name << " (" << seqs[r->seq_num].length
         << " bp)\n"
         << "Location: "
         << seqs[r->seq_num].genomic_pos + r->start + 1 << " to "
         << seqs[r->seq_num].genomic_pos + r->end + 1 << '\n'
         << "Score: " << r->score << '\n';
    vector<pair<double, string> > x;
    for (uint m = 0; m != mat_names.size(); ++m)
      x.push_back(make_pair(r->motif_scores[m], mat_names[m]));
    sort(x.begin(), x.end(), greater<pair<double, string> >());
    for (vector<pair<double, string> >::const_iterator i = x.begin();
         i != x.end(); ++i)
      strm << i->second << ": " << i->first << '\n';
    strm << r->cluster_seq << '\n';
    print_hits(strm, seqs[r->seq_num], r->hits);
    strm << endl;
  }
}

void cb::output_by_score_concise(ostream &strm, const vector<seq_info> &seqs) {
  strm << "# Score\tStart\tEnd\tSequence";
  for (vector<string>::const_iterator m = mat_names.begin();
       m != mat_names.end(); ++m)
    strm << "\t" << *m;
  strm << '\n';

  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    strm << r->score << "\t"
         << seqs[r->seq_num].genomic_pos + r->start + 1 << "\t"
         << seqs[r->seq_num].genomic_pos + r->end + 1 << "\t"
         << seqs[r->seq_num].name;
    for (vector<double>::const_iterator m = r->motif_scores.begin();
         m != r->motif_scores.end(); ++m)
      strm << "\t" << *m;
    strm << '\n';
  }

  strm << endl;
}

void cb::output_sequence_name_sorted_by_score(ostream &strm, const vector<seq_info> &seqs) {
  strm << "# Sequence name\tScore\tSequence number\tRank\n";

  uint rank_position = 1;

  for (vector<result>::const_iterator r = results.begin(); r != results.end();
       ++r) {
    strm << seqs[r->seq_num].name << "\t"
         << r->score << "\t"
         << r->seq_num << "\t"
         << rank_position << "\n";
    rank_position++;
  }

  strm << std::flush;
}

// I'll probably move this to a library
inline std::ifstream &open_or_die(const std::string &filename,
                                  std::ifstream &strm) {
  strm.open(filename.c_str());
  if (!strm)
    mcf::die("Sorry, couldn't open " + filename);
  return strm;
}

int main(int argc, char **argv) {
  using namespace std;

  args::parse(argc, argv);
  cb::get_matrices();
  cb::misc_init();

  vector<cb::seq_info> seqs;
  unsigned (*translator)(char) =
      (args::mask_lower ? nolower_translator : mcf::DNA_to_number);
  ifstream file;
  istream &in =
      args::seqfile.empty() ? (istream &)cin : open_or_die(args::seqfile, file);
  string seq_name;
  bool by_sequence = args::out_format == args::BY_SEQUENCE ||
                     args::out_format == args::BY_SEQUENCE_CONCISE ||
                     args::out_format == args::BED;
  cout.setf(ios::left, ios::adjustfield);
  cout.precision(3); // 3 sig figs

  // read and analyze the sequences 1 by 1:
  for (unsigned seq_num = 0; mcf::get_fasta(in, cb::seq, seq_name, translator);
       ++seq_num) {
    { // the "swap trick" to save memory:
      vector<unsigned> temp(cb::seq);
      temp.swap(cb::seq);
    } // the one-liner version of this trick confuses my compiler
    istringstream is(seq_name);
    is >> seq_name; // get first word (?)

    if (args::genomic_coordinates) {
        seqs.push_back(cb::get_chrom_and_pos(seq_name, cb::seq.size(), args::zero_based));
    } else {
        seqs.push_back(cb::seq_info(seq_name, cb::seq.size()));
    }

    if (args::verbose) {
      cout << std::flush;
      cerr << "Scanning " << seq_name << " (" << cb::seq.size() << " bp)..."
           << endl;
    }

    cb::scan_seq(seq_num);

    if (by_sequence && !cb::results.empty()) {
      sort(cb::results.begin(), cb::results.end(), byscore<cb::result>());

      if (args::out_format == args::BY_SEQUENCE) {
        cb::output_by_seq(cout, seqs.back());
      } else if (args::out_format == args::BY_SEQUENCE_CONCISE) {
        cb::output_by_seq_concise(cout, seqs.back());
      } else if (args::out_format == args::BED) {
        cb::output_by_seq_bed(cout, seqs.back());
      }

      cb::results.clear();
    }
    cb::seq.clear();
  }

  if (args::verbose) {
    cout << std::flush;
    cerr << endl;
  }

  if ((args::out_format == args::BY_SCORE ||
       args::out_format == args::BY_SCORE_CONCISE ||
       args::out_format == args::SEQUENCE_NAME_SORTED_BY_SCORE) &&
      !cb::results.empty()) {
    sort(cb::results.begin(), cb::results.end(), byscore<cb::result>());

    if (args::out_format == args::BY_SCORE) {
      cb::output_by_score(cout, seqs);
    } else if (args::out_format == args::BY_SCORE_CONCISE) {
      cb::output_by_score_concise(cout, seqs);
    } else if (args::out_format == args::SEQUENCE_NAME_SORTED_BY_SCORE) {
      cb::output_sequence_name_sorted_by_score(cout, seqs);
    }
  }

  if (args::out_format != args::SEQUENCE_NAME_SORTED_BY_SCORE &&
      args::out_format != args::BED) {
      cout.precision(6); // reset to default precision
      args::print(cout, seqs.size(),
                  cb::mat_names.size()); // print command line arguments
  }

  cout << std::flush;
}
