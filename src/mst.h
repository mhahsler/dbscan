#ifndef MST_H
#define MST_H

#include <Rcpp.h>
using namespace Rcpp;

#include <queue>
#include <bits/stdc++.h>
#include "lt.h"

// Functions to compute MST and build hclust object out of resulting tree
NumericMatrix mst(const NumericVector x_dist, const R_xlen_t n);

List hclustMergeOrder(NumericMatrix mst, IntegerVector o);

#endif
