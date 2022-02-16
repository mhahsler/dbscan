#include <Rcpp.h>
using namespace Rcpp;

// Computes the mutual reachability distance defined for HDBSCAN
//
// The mutual reachability distance is a summary at what level two points together
// will connect. The mutual reachability distance is defined as:
// mrd(a, b) = max[core_distance(a), core_distance(b), distance(a, b)]
//
// Input:
// * dm: distances as a dist object (vector) of size (n*(n-1))/2 where n
//       is the number of points.
//       Note: we divide by 2 early to stay within the numb10er range of int.
// * cd: the core distances as a vector of length n
//
// Returns:
// a vector (dist object) in the same order as dm
// [[Rcpp::export]]
NumericVector mrd(NumericVector dm, NumericVector cd) {
  R_xlen_t n = cd.length();
  if (dm.length() != (n / 2 * (n-1)))
    stop("number of mutual reachability distance values and size of the distances do not agree.");

  NumericVector res = NumericVector(dm.length());
  for (R_xlen_t i = 0, idx = 0; i < n; ++i) {
    for (R_xlen_t j = i+1; j < n; ++j, ++idx) {
      res[idx] = std::max(dm[idx], std::max(cd[i], cd[j]));
    }
  }
  return res;
}
