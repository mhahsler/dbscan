//----------------------------------------------------------------------
//                                DBSCAN density
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3 or later
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
#include "ann_raii.h"
#include "ANN/ANN.h"
#include "regionQuery.h"

using namespace Rcpp;

// faster implementation of counting point densities from a matrix
// using a kd-tree
// [[Rcpp::export]]
IntegerVector dbscan_density_int(
    NumericMatrix data, double eps,
    int type, int bucketSize, int splitRule, double approx) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;

  ANNPoints dataPts;
  int nrow = NA_INTEGER;
  int ncol= NA_INTEGER;

  // copy data for kd-tree
  nrow = data.nrow();
  ncol = data.ncol();
  dataPts = ANNPoints(nrow, ncol);
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  auto kdTree = make_ann_tree(dataPts.get(), nrow, ncol, type,
    bucketSize, splitRule);
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  std::vector<int> N;
  IntegerVector count(nrow);

  for (int i=0; i<nrow; ++i) {
    //Rprintf("processing point %d\n", i+1);
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    N = regionQuery(i, dataPts.get(), kdTree.get(), eps2, approx);
    count[i] = N.size();
  }

  return count;
}
