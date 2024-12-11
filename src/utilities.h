//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#ifndef UTILITIES_H
#define UTILITIES_H

#include <Rcpp.h>
using namespace Rcpp;

// std::to_string is apparently a c++11 only thing that crashes appveyor, so using ostringstream it is!
namespace patch {
template <typename T> std::string to_string(const T& n) {
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
} // namespace patch

template <typename T, typename C>
bool contains(const T& container, const C& key) {
  if (std::find(container.begin(), container.end(), key) != container.end()) {
    return true;
  }
  return false;
}

// [[Rcpp::export]]
IntegerVector lowerTri(IntegerMatrix m);

// [[Rcpp::export]]
NumericVector combine(const NumericVector& t1, const NumericVector& t2);

IntegerVector combine(const IntegerVector& t1, const IntegerVector& t2);

// Faster version of above combine function, assuming you can precompute and store
// the containers needing to be concatenated
// [[Rcpp::export]]
IntegerVector concat_int(List const& container);

#endif
