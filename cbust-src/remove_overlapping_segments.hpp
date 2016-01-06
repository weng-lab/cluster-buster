#include <set>
#include <vector>
#include <cassert>

// This function takes a vector, 'in', of "segments".
// A segment can be any type that contains data members 'start' and 'end'.
// It returns the segments, in 'out', with overlapping ones removed.
// The 'in' vector is assumed to be sorted in order of preference,
// so that when 2 segments overlap, the one later in 'in' will be removed.
// The returned segments are sorted by position.

namespace mcf {
  template<class T>
  void remove_overlapping_segments(const std::vector<T> & in, std::vector<T> & out);
}

namespace {  // anonymous namespace: accessible from this file only
  template<class T>
  struct lesspos {  // sort criterion for segments
    bool operator()(const T & a, const T & b) const
    { return a.end < b.start; }
  };
}

template<class T>
void mcf::remove_overlapping_segments(const std::vector<T> & in, std::vector<T> & out)
{
  std::set<T, lesspos<T> > s;
  // set is better than vector + binary_search, because vector insertion
  // is slow.

  for (typename std::vector<T>::const_iterator i = in.begin(); i != in.end(); ++i) {
    assert(i->start <= i->end);
    typename std::set<T, lesspos<T> >::iterator j = s.lower_bound(*i);
    if (j == s.end() || lesspos<T>()(*i, *j))
      s.insert(j, *i);
  }

  out.assign(s.begin(), s.end());
}
