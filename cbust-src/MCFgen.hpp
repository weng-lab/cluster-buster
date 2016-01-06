/***** Useful, non-bio-specific C++ functions *****/
/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

// Basically, these are functions that I think should be
// in the standard library but aren't
// Better versions of many of these functions probably
// exist in libraries such as boost

#ifndef MCF_GEN_H
#define MCF_GEN_H

#include <iostream>  // die
#include <iterator>
#include <numeric>  // accumulate
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>

namespace mcf {
  // convert any type to a string, with default formatting
  template<class T> std::string tostring(T x);

  // perl-style die
  void die(const std::string & message);

  // return x to the power of y
  // only use this when S, T are integral types and y >= 0
  // doesn't check for overflow!!
  // pow() is dangerous for integers because it uses floats
  template<class S, class T> S int_pow(S x, T y);

  template<class T> void reserve_or_die(std::vector<T> & v, unsigned n);

  // return true if the range equals its reverse
  template<class It> bool is_reverse(It start, It end);

  // normalize a range (make it sum to 1)
  // no guarantee that the answer will be exactly 1
  // the range had better contain floating-point values!
  template<class It> void normalize(It start, It end);
}

template<class T> inline std::string mcf::tostring(T x)
{
  std::ostringstream temp;
  temp << x;
  return temp.str();
}

inline void mcf::die(const std::string & message)
{
  std::cerr << message << std::endl;
  exit(EXIT_FAILURE);
}

template<class S, class T> inline S mcf::int_pow(S x, T y)
{
  S ans = 1;
  for (; y > 0; --y)
    ans *= x;
  return ans;
}

template<class T>
inline void mcf::reserve_or_die(std::vector<T> & v, unsigned n)
{
  v.reserve(n);
  if (v.capacity() < n)
    die("Out of memory: couldn't reserve " + tostring(n) +
	" vector elements of size " + tostring(sizeof(T)));
}

template<class It> bool mcf::is_reverse(It start, It end)
{
  while (end > start)
    if (*(start++) != *(--end))
      return false;

  return true;
}

template<class It> void mcf::normalize(It start, It end)
{
  typename std::iterator_traits<It>::value_type tot =
    std::accumulate(start, end,  // added typename 16-2-2005:
		    typename std::iterator_traits<It>::value_type(0));

  assert(tot != 0);  // doesn't like being prefixed by std::

  tot = 1/tot;
  for (; start < end; ++start)
    *start *= tot;
}

/* Functions below here deprecated */

#include <cstdlib>  // malloc, rand, RAND_MAX

// use my matrix class instead!!
// allocate an x by y matrix:
// allocates just 1 chunk of memory, so can free it with free()
// (appalling low-level hackery!)
// is it possible to use new instead of malloc ???
template<class T> T ** new_matrix(unsigned x, unsigned y)
{
  T ** mat = (T **) malloc(x * (sizeof(void *) + y * sizeof(T)));

  T * p = (T *)(mat + x);

  for (unsigned i = 0u; i < x; ++i, p += y)
    mat[i] = p;

  return mat;
}

/*** functions for making random choices ***/
// random stuff is tricky!!! Use boost instead(?)

// return a random double between 0 (inclusive) and n (exclusive)
inline double rand(double n)
{
  return rand() / (RAND_MAX + 1.0) * n;
}

// return a random float between 0 (inclusive) and n (exclusive)
// the loop seems to be needed for float, but not double
inline float rand(float n)
{
  float r;

  do {
    r = float(rand(double(n)));
  } while (r == n);

  return r;
}

// return a random unsigned int between 0 (inclusive) and n (exclusive)
// from C FAQ, 13.16
inline unsigned rand(unsigned n)
{
  return unsigned(rand(double(n)));
}

// randomly choose an element in the range [start end)
// weighted by the values of the elements
// the total of the element values is passed in & assumed to be correct
// the "end" is theoretically redundant, but if T is a floating type,
// the "total" might be slightly imprecise
template<class It, class T> It random_choice(It start, It end, T total)
{
  // choose a random number between 0 (inclusive) and total (exclusive):
  T choice = rand(total);

  T x = T(0);
  --end;  // if we get to the last element, we will choose it for sure

  while (start < end) {
    x += *start;
    if (x > choice)
      break;
    ++start;
  }

  return start;
}

#endif
