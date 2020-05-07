# Cluster-Buster

The main application is detection of sequences that regulate gene transcription,
such as enhancers and silencers, but other types of biological regulation may be
mediated by motif clusters too.

## Original webiste

[Cluster-Buster overview page](http://cagt.bu.edu/page/ClusterBuster_about)

[Cluster-Buster download page](http://cagt.bu.edu/page/ClusterBuster_download)

[Cluster-Buster github page](https://github.com/weng-lab/cluster-buster)


## Publication

[Frith MC, Li MC, Weng Z (2003) Cluster-Buster: Finding dense clusters of motifs in DNA sequences. Nucleic Acids Res. 2003 Jul 1; 31(13):3666-8](http://dx.doi.org/10.1093/nar/gkg540)

## Contributors

* Martin Frith (mcfrith): m.frith ATaist.go.jp
* Gert Hulselmans (ghuls): hulselmansgert ATgmail.com

## Install

```bash
# Clone Cluster-Buster git repository
git clone https://github.com/weng-lab/cluster-buster

cd cluster-buster

# Compile Cluster-Buster with g++.
make cbust

# Compile Cluster-Buster with icpc (Intel compiler: log1p, exp and log might be faster).
make cbust_intel
```


## Help

```
$ cbust
Usage: cbust [options] matrixfile myseqs.fa

Options:
-h Help: print documentation
-V Show version
-v Verbose
-c Cluster score threshold (5)
-m Motif score threshold (6)
-g Expected gap (bp) between neighboring motifs in a cluster (35)
-r Range in bp for counting local nucleotide abundances (100)
-b Background padding in bp (0)
-l Mask lowercase letters
-p Pseudocount (0.375)
-t Keep top X clusters per sequence (0 (= all))
-G Use genomic coordinates (extracted from sequence name)
   0: zero-based start coordinate
   1: one-based start coordinate
-f Output format (0)
   0: per sequence (default)
   1: per sequence, concise format
   2: sorted by cluster score
   3: sorted by cluster score, concise format
   4: sorted by cluster score: seq name, score, seq number, rank
   5: BED file
```

```
$ cbust -h
Cluster-Buster is a program for finding interesting functional regions,
such as transcriptional enhancers, in DNA sequences. Its main premise is
that the functional region contains multiple smaller elements, e.g.
protein-binding sites, which thus occur in a cluster in the DNA sequence.
The user needs to specify a list of motifs (short sequence patterns): the
program will then look for dense clusters of these motifs in the sequence.

Cluster-Buster requires two inputs: a file of DNA sequences in the standard
FASTA format, and a file of motifs. Any non-alphabetic characters in the
sequences are ignored, and any alphabetic characters except A, C, G, T
(uppercase or lowercase) are converted to 'n' and forbidden from matching
motifs.

The motif file should contain matrices in the following format:

>TATA
1  0  0  19
20 0  0  0
2  1  0  17
19 0  1  0
>CCAAT
(etc.)

The rows of each matrix correspond to successive positions of the motif,
from 5' to 3', and the columns indicate the frequencies of A, C, G, and T,
respectively, in each position. These frequencies are usually obtained
from alignments of protein-binding sites.

It is possible to assign 'weights' to the motifs indicating how important
they are for defining clusters. To do so, place a line like this within the
motif definition (below the line beginning '>'):

# WEIGHT 2.7

This motif will carry 2.7-fold more weight than a motif with weight 1 (the
default value). It is also possible to specify the gap parameter (see
below) in the motif file, by including a line like this:

# GAP 22.4

Specifying the gap parameter with the -g option (see below) overrides any
value given in the motif file. The program Cluster-Trainer can be used to
estimate good values for the weights and the gap parameter.

Cluster-Buster compares each matrix to every location in the DNA sequence,
and calculates scores reflecting the goodness of the match. It then
identifies motif clusters as sequence regions with unusually strong
concentrations of high-scoring matches. Each motif cluster receives a score
which depends on the scores of its motifs, their weights, the tightness of
their clustering, and the manner in which they overlap one another.

The location, score, and sequence of each motif cluster is printed. In
addition, the contribution of each motif type to each cluster is indicated
by a score. (This is the cluster score minus the score excluding that motif
type.) Finally, for high-scoring motif matches within clusters, the motif
type, location, strand, score and sequence is printed.

The program's behavior may be modified with the following options. The
default values are designed to give sensible results.

-h Help: print documentation. You already know this one.
-V Show version.
-v Verbose: print sequence names which are being scanned.
-c Cluster score threshold (default = 5).
   Print details of all motif clusters with score >= this value.
-m Motif score threshold (default = 6).
   Print details of all motif matches that have score >= this value and
   occur within printed clusters. This option has no effect on finding
   motif clusters: it just affects which motifs get printed.
-g Gap parameter (default = 35).
   The expected distance in bp between neighboring motifs in a cluster. Low
   values will make the program more sensitive to finding very tight
   clusters even if the motif scores are not so high, and high values will
   make it more sensitive to finding clusters of high-scoring motifs even
   if the clustering is not very tight.
-r Range in bp for counting local nucleotide abundances (default = 100).
   Abundances of A, C, G, and T vary significantly along natural DNA
   sequences. The program estimates the local nucleotide abundances at each
   position in the sequence by counting them up to this distance in both
   directions. These abundances affect motif and hence cluster scores: e.g.
   a GC-rich motif will receive a higher score if found in an AT-rich
   region, where it is more surprising, than if it is found in a GC-rich
   region.
-b Background padding in bp (default = 0).
   Consider first X bp and last X bp of each sequence only for counting
   local base abundances but not for scoring to find motif clusters.
   This can be useful when scoring small regions.
-l Mask lowercase letters in the sequences (i.e. forbid motifs from
   matching them). Lowercase letters are often used to indicate repetitive
   regions. Repetitive regions, especially tandem repeats, can produce
   extremely high-scoring motif clusters, which may or may not be spurious.
-p Pseudocount (default = 0.375).
   This value gets added to all entries in the motif matrices. Pseudocounts
   are a standard way of estimating underlying frequencies from a limited
   number of observations. If your matrices contain probabilities rather
   than counts, you should probably set this parameter to zero.
-t Keep X top clusters per sequence (default = 0 (= all)).
   If set to 1 or higher, keep only this amount of best scoring clusters
   above the cluster threshold for each sequence. If set to 0, keep all
   clusters above the cluster threshold for each sequence.
-G Extract genomic coordinates from sequence name and output genomic
   coordinates instead of relative coordinates.
   Examples of sequence names from which chromosome names and start
   positions can be extracted:
     - chr10:123456
     - chr10:123456-234567
     - chr10:123456@@gene_name
     - chr10:123456-234567@@gene_name
   Specify if the start coordinate is zero- or one-based:
     0: zero-based (as generated by: bedtools getfasta)
     1: one-based
-f Output format (default = 0).
     0: Print the clusters in the first sequence sorted by score, then the
        clusters in the second sequence sorted by score, etc.
     1: Concise version of 0, omitting details of individual motif matches.
     2: Sort all clusters by score, regardless of which sequence they come
        from.
     3: Concise version of 2, omitting details of individual motif matches.
     4: Sort all clusters by score and output sequence name, score, number
        of the sequence in the FASTA file, ranking of the score.
     5: BED file with all info.

Example usage: cbust -g 20 -l mymotifs myseqs.fa

For more information on the Cluster-Buster algorithm, see:
Cluster-Buster: Finding dense clusters of motifs in DNA sequences
by MC Frith, MC Li & Z Weng in Nucleic Acids Research 31:3666-8 2003
and references therein.

Comments and questions to Martin Frith: m.frith ATaist.go.jp

```
