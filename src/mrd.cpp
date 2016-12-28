#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix mrd(NumericMatrix dm, NumericVector cd) {
  int nrow = dm.nrow(), ncol = dm.ncol();
  NumericMatrix out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      out(i, j) = i == j ? cd[i] : std::max(dm(i, j), std::max(cd[i], cd[j]));
    }
  }
  return out;
}