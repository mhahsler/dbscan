#ifndef MST_H
#define MST_H

#include <Rcpp.h>
#include "lt.h"

using namespace Rcpp;

// Functions to compute MST and build hclust object out of the resulting tree
NumericMatrix mst_dist_cpp(const NumericVector x_dist, const R_xlen_t n);

NumericMatrix mst_mrd(const NumericMatrix x, const NumericVector coredist);

List hclustMergeOrder(NumericMatrix mst, IntegerVector o);

#endif
