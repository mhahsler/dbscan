#include <Rcpp.h>
using namespace Rcpp;

// Computes the mutual reachability distance from the the distance vector 'dm' and the core distance 
// [[Rcpp::export]]
NumericVector mrd(NumericVector dm, NumericVector cd) {
  NumericVector res = NumericVector(dm.length());
  int n = cd.length();
  for (int i = 0, c_ = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      res[c_] = std::max(dm[c_], std::max(cd[i], cd[j]));
      c_++;
    }
  }
  return res;
}

// Full Matrix version
// [[Rcpp::export]]
NumericMatrix mrd_m(NumericMatrix dm, NumericVector cd) {
  int nrow = dm.nrow(), ncol = dm.ncol();
  NumericMatrix out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out(i, j) = i == j ? cd[i] : std::max(dm(i, j), std::max(cd[i], cd[j]));
    }
  }
  return out;
}