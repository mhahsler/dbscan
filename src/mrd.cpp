//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

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
//       Note: we divide by 2 early to stay within the number range of int.
// * cd: the core distances as a vector of length n
//
// Returns:
// a vector (dist object) in the same order as dm
// [[Rcpp::export]]
NumericVector mrd(NumericVector dm, NumericVector cd) {
  R_xlen_t n = cd.length();
  if (dm.length() != (n * (n-1) / 2))
    stop("number of mutual reachability distance values and size of the distance matrix do not agree.");

  NumericVector res = NumericVector(dm.length());
  for (R_xlen_t i = 0, idx = 0; i < n; ++i) {
//    Rprintf("i = %ill of %ill, idx = %ill\n", i, n, idx);
    for (R_xlen_t j = i+1; j < n; ++j, ++idx) {
      res[idx] = std::max(dm[idx], std::max(cd[i], cd[j]));
    }
  }
  return res;
}
