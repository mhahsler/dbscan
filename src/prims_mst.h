#ifndef PRIMS_MST_H
#define PRIMS_MST_H

#include <Rcpp.h>
using namespace Rcpp;

#include <queue> // priority_queue

// Functions to compute MST and build hclust object out of resulting tree
NumericMatrix prims(const NumericVector x_dist, const R_xlen_t n);

List hclustMergeOrder(NumericMatrix mst, IntegerVector o);

#endif
