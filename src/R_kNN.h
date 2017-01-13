#ifndef R_KNN_H
#define R_KNN_H

#include <Rcpp.h>
#include "ANN/ANN.h"

using namespace Rcpp;

// returns knn + dist
// [[Rcpp::export]]
List kNN_int(NumericMatrix data, int k,
             int type, int bucketSize, int splitRule, double approx);

#endif