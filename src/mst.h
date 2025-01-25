#ifndef MST_H
#define MST_H

#include <Rcpp.h>
#include "lt.h"

using namespace Rcpp;

// Functions to compute MST and build hclust object out of the resulting tree
NumericMatrix mst(const NumericVector x_dist, const R_xlen_t n);

List hclustMergeOrder(NumericMatrix mst, IntegerVector o);

#endif
