#include <ostream>
#include <unistd.h>
#include "MCFgen.hpp"  // tostring, die
#include "args.hpp"

namespace args {
  // definitions & defaults:
  string seqfile, matfile;
  double score_thresh = 5;  // ??
  double motif_thresh = 6;  // gives a reasonable number of motif predictions
  double e_gap = 35;
  bool gap_specified = false;
  format out_format = BY_SEQUENCE;
  uint bg_range = 100;
  bool mask_lower = false;
  double pseudo = 0.375;
  double tau = 0;
}

void args::parse(int argc, char ** argv)
{
  const string doc =
"Cluster-Buster is a program for finding interesting functional regions,\n"
"such as transcriptional enhancers, in DNA sequences. Its main premise is\n"
"that the functional region contains multiple smaller elements, e.g.\n"
"protein-binding sites, which thus occur in a cluster in the DNA sequence.\n"
"The user needs to specify a list of motifs (short sequence patterns): the\n"
"program will then look for dense clusters of these motifs in the sequence.\n"
"\n"
"Cluster-Buster requires two inputs: a file of DNA sequences in the standard\n"
"FASTA format, and a file of motifs. Any non-alphabetic characters in the\n"
"sequences are ignored, and any alphabetic characters except A, C, G, T\n"
"(uppercase or lowercase) are converted to 'n' and forbidden from matching\n"
"motifs.\n"
"\n"
"The motif file should contain matrices in the following format:\n"
"\n"
">TATA\n"
"1  0  0  19\n"
"20 0  0  0\n"
"2  1  0  17\n"
"19 0  1  0\n"
">CCAAT\n"
"(etc.)\n"
"\n"
"The rows of each matrix correspond to successive positions of the motif,\n"
"from 5' to 3', and the columns indicate the frequencies of A, C, G, and T,\n"
"respectively, in each position. These frequencies are usually obtained\n"
"from alignments of protein-binding sites.\n"
"\n"
"It is possible to assign 'weights' to the motifs indicating how important\n"
"they are for defining clusters. To do so, place a line like this within the\n"
"motif definition (below the line beginning '>'):\n"
"\n"
"# WEIGHT 2.7\n"
"\n"
"This motif will carry 2.7-fold more weight than a motif with weight 1 (the\n"
"default value). It is also possible to specify the gap parameter (see\n"
"below) in the motif file, by including a line like this:\n"
"\n"
"# GAP 22.4\n"
"\n"
"Specifying the gap parameter with the -g option (see below) overrides any\n"
"value given in the motif file. The program Cluster-Trainer can be used to\n"
"estimate good values for the weights and the gap parameter.\n"
"\n"
"Cluster-Buster compares each matrix to every location in the DNA sequence,\n"
"and calculates scores reflecting the goodness of the match. It then\n"
"identifies motif clusters as sequence regions with unusually strong\n"
"concentrations of high-scoring matches. Each motif cluster receives a score\n"
"which depends on the scores of its motifs, their weights, the tightness of\n"
"their clustering, and the manner in which they overlap one another.\n"
"\n"
"The location, score, and sequence of each motif cluster is printed. In\n"
"addition, the contribution of each motif type to each cluster is indicated\n"
"by a score. (This is the cluster score minus the score excluding that motif\n"
"type.) Finally, for high-scoring motif matches within clusters, the motif\n"
"type, location, strand, score and sequence is printed.\n"
"\n"
"The program's behavior may be modified with the following options. The\n"
"default values are designed to give sensible results.\n"
"\n"
"-h Help: print documentation. You already know this one.\n"
"-c Cluster score threshold (default = " + mcf::tostring(score_thresh) + ").\n"
"   Print details of all motif clusters with score >= this value.\n"
"-m Motif score threshold (default = " + mcf::tostring(motif_thresh) + ").\n"
"   Print details of all motif matches that have score >= this value and\n"
"   occur within printed clusters. This option has no effect on finding\n"
"   motif clusters: it just affects which motifs get printed.\n"
"-g Gap parameter (default = " + mcf::tostring(e_gap) + ").\n"
"   The expected distance in bp between neighboring motifs in a cluster. Low\n"
"   values will make the program more sensitive to finding very tight\n"
"   clusters even if the motif scores are not so high, and high values will\n"
"   make it more sensitive to finding clusters of high-scoring motifs even\n"
"   if the clustering is not very tight.\n"
"-r Range in bp for counting local nucleotide abundances (default = " + mcf::tostring(bg_range) + ").\n"
"   Abundances of A, C, G, and T vary significantly along natural DNA\n"
"   sequences. The program estimates the local nucleotide abundances at each\n"
"   position in the sequence by counting them up to this distance in both\n"
"   directions. These abundances affect motif and hence cluster scores: e.g.\n"
"   a GC-rich motif will receive a higher score if found in an AT-rich\n"
"   region, where it is more surprising, than if it is found in a GC-rich\n"
"   region.\n"
"-l Mask lowercase letters in the sequences (i.e. forbid motifs from\n"
"   matching them). Lowercase letters are often used to indicate repetitive\n"
"   regions. Repetitive regions, especially tandem repeats, can produce\n"
"   extremely high-scoring motif clusters, which may or may not be spurious.\n"
"-p Pseudocount (default = " + mcf::tostring(pseudo) + ").\n"
"   This value gets added to all entries in the motif matrices. Pseudocounts\n"
"   are a standard way of estimating underlying frequencies from a limited\n"
"   number of observations. If your matrices contain probabilities rather\n"
"   than counts, you should probably set this parameter to zero.\n"
"-f Output format (default = " + mcf::tostring(out_format) + ").\n"
"   0: Print the clusters in the first sequence sorted by score, then the\n"
"      clusters in the second sequence sorted by score, etc.\n"
"   1: Concise version of 0, omitting details of individual motif matches.\n"
"   2: Sort all clusters by score, regardless of which sequence they come\n"
"      from.\n"
"   3: Concise version of 2, omitting details of individual motif matches.\n"
"\n"
"Example usage: cbust -g20 -l mymotifs myseqs.fa\n"
"\n"
"For more information on the Cluster-Buster algorithm, see:\n"
"Cluster-Buster: Finding dense clusters of motifs in DNA sequences\n"
"by MC Frith, MC Li & Z Weng in Nucleic Acids Research 31:3666-8 2003\n"
"and references therein.\n"
"\n"
"Comments and questions to Martin Frith: mfrith@zlab.bu.edu.\n"
;

  const string usage =
"Usage: cbust [options] matrixfile myseqs.fa\n"
"\n"
"Options:\n"
"-h Help: print documentation\n"
"-c Cluster score threshold (" + mcf::tostring(score_thresh) + ")\n"
"-m Motif score threshold (" + mcf::tostring(motif_thresh) + ")\n"
"-g Expected gap (bp) between neighboring motifs in a cluster (" + mcf::tostring(e_gap) + ")\n"
"-r Range in bp for counting local nucleotide abundances (" + mcf::tostring(bg_range) + ")\n"
"-l Mask lowercase letters\n"
"-p Pseudocount (" + mcf::tostring(pseudo) + ")\n"
"-f Output format (" + mcf::tostring(out_format) + ")\n"
"   0: per sequence (default)\n"
"   1: per sequence, concise format\n"
"   2: sorted by cluster score\n"
"   3: sorted by cluster score, concise format\n"
    //"-e  transition probability to HMM end state (tau) (" + mcf::tostring(tau) + ")\n"
;

  int c;

  while ((c = getopt(argc, argv, "hc:m:g:f:r:lp:e:")) != -1)
    switch (c) {
    case 'h':
      cout << doc << endl;
      exit(0);
    case 'c':
      score_thresh = atof(optarg);
      if (score_thresh < 0) mcf::die("Cluster score threshold should be >= zero");
      break;
    case 'm':
      motif_thresh = atof(optarg);
      break;
    case 'g':
      e_gap = atof(optarg);
      gap_specified = true;
      if (e_gap < 0) mcf::die("Expected gap should be >= zero");
      break;
    case 'f':
      switch (atoi(optarg)) {
      case 0: out_format = BY_SEQUENCE; break;
      case 1: out_format = BY_SEQUENCE_CONCISE; break;
      case 2: out_format = BY_SCORE; break;
      case 3: out_format = BY_SCORE_CONCISE; break;
      default: mcf::die("Format should be 0, 1, 2, or 3");
      }
      break;
    case 'r':
      bg_range = atoi(optarg);  // what happens if user enters -ve number?
      break;
    case 'l':
      mask_lower = true;
      break;
    case 'p':
      pseudo = atof(optarg);
      if (pseudo < 0) mcf::die("Pseudocount should be >= zero");
      break;
    case 'e':
      tau = atof(optarg);
      if (tau < 0 || tau >= 1) mcf::die("Tau should be >= 0 and < 1");
      break;
    case '?':
      mcf::die(usage);
    }

  if (optind >= argc)  // there should be at least 1 non-option argument
    mcf::die(usage);

  matfile = argv[optind++];
  if (optind < argc)
    seqfile = argv[optind++];
}

void args::print(ostream & strm, uint seq_num, uint mat_num)
{
  strm
    << "Sequence file: " << seqfile << " (" << seq_num << " sequences)\n"
    << "Matrix file: " << matfile << " (" << mat_num << " matrices)\n"
    << "Expected gap: " << e_gap << '\n'
    << "Range for local abundances: " << bg_range << '\n'
    << "Lowercase filtering " << (mask_lower ? "ON\n" : "OFF\n")
    << "Pseudocount: " << pseudo << '\n'
    << "Cluster score threshold: " << score_thresh << '\n'
    << "Motif score threshold: " << motif_thresh << '\n'
    //    << "Tau: " << tau << '\n'
    << flush;
}
