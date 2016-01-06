// bioinformatics routines of general use
/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#ifndef MCF_BIO_H
#define MCF_BIO_H

#include <iosfwd>  // "forward declaration" of istream
#include <string>
#include <vector>
#include <cassert>
#include <cmath>  // log

#include "matrix.hpp"  // my 2D matrix class

// Current policy: represent sequence residues as unsigned ints
// (consider unsigned chars?)

namespace mcf {
  char number_to_DNA(unsigned b);  // Translates 01234 to acgtn
  unsigned DNA_to_number(char c);  // aA->0, cC->1, gG->2, tT->3, other->4

  // Industrial strength fasta-format sequence reader
  // On error, seq and title are unchanged
  std::istream & get_fasta(
    std::istream & strm,          // stream to read from
    std::vector<unsigned> & seq,  // APPEND the sequence to this vector
    std::string & title,          // store title line minus leading '>'
    unsigned (*translator)(char)  // function translating chars to uints
    = DNA_to_number);             // default translator function

  // count residue types in a sequence
  // any residues with value >= alphsize will be ignored
  // if counts has < alphsize elements, it is resized
  void count_residues(
    const std::vector<unsigned> & seq,  // sequence who's residues to count
    std::vector<unsigned> & counts,     // place to store counts
    unsigned alphsize);                 // size of alphabet

  // Industrial strength reader for fasta-format PSSMs/PWMs/profiles
  // On error, mat and title are unchanged
  std::istream & get_simple_pssm(
    std::istream & strm,     // stream to read from
    matrix<float> & mat,     // store the matrix here
    std::string & title,     // store title line minus leading '>'
    unsigned alphsize = 4);  // perhaps it should guess alphsize from the input...

  // Like get_simple_pssm, but if a line beginning #WEIGHT is present,
  // reads a numeric weight from this line.
  std::istream & get_weighted_pssm(
    std::istream & strm,     // stream to read from
    matrix<float> & mat,     // store the matrix here
    std::string & title,     // store title line minus leading '>'
    float & weight,          // = 1 by default if not specified
    unsigned alphsize = 4);  // perhaps it should guess alphsize from the input...

  // normalize each row of the matrix (i.e. make each row sum to 1)
  // assumes no row has all zeros
  void normalize_pssm(
    matrix<float> & p,  // the matrix to normalize
    const std::vector<float> & pseudos);  // add these to matrix cells before normalizing
                                     // each residue type can have a different pseudocount

  // a very low number, but not so low that we risk overflow:
  const float minus_infinity = -1e25f;

  // replace each matrix cell with its natural log
  // log(zero) -> minus_infinity
  // assumes there are no negative numbers!
  void log_pssm(matrix<float> & p);

  // Count oligonucleotides of given length in a sequence
  // The sequence should consist of integers between 0 and alphsize-1, inclusive
  // Integers >= alphsize, e.g. representing masked repetitive sequence,
  // are treated as breaking the sequence into separate fragments
  // Oligos are encoded in a fairly obvious way:
  // E.g. the oligonucleotide (1, 3, 2) is encoded as
  // 2 + 3 * alph_size + 1 * alph_size^2
  // This code gives the index of the counts array that will contain each oligo's counts
  // counts is assumed to be pre-initialized
  void count_oligos(
    const std::vector<unsigned> & seq,  // sequence who's oligos to count
    std::vector<unsigned> & counts,     // place to store counts
    unsigned oli_len,                   // length of oligos to count
    unsigned alphsize);                 // size of alphabet
}

inline char mcf::number_to_DNA(unsigned b)
{
  static const char lookup[] = "acgtn";
  assert(b < 5);
  return lookup[b];
}

// should use static array instead of switch - ask Jianhua how
inline unsigned mcf::DNA_to_number(char c)
{
  switch(c) {
  case 'a':
  case 'A':
    return 0u;
  case 'c':
  case 'C':
    return 1u;
  case 'g':
  case 'G':
    return 2u;
  case 't':
  case 'T':
    return 3u;
  default:
    return 4u;
  }
}

inline void mcf::log_pssm(matrix<float> & p) 
{
  for (std::vector<float>::iterator i = p[0]; i < p[p.rows()]; ++i) {
    assert(*i >= 0);
    if (*i == 0)
      *i = minus_infinity;
    else
      *i = log(*i);
  }
}

#endif
