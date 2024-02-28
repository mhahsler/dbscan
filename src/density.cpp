//----------------------------------------------------------------------
//                                DBSCAN density
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
#include "ANN/ANN.h"
#include "regionQuery.h"

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector dbscan_density_int(
    NumericMatrix data, double eps,
    int type, int bucketSize, int splitRule, double approx) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;

  ANNpointSet* kdTree = NULL;
  ANNpointArray dataPts = NULL;
  int nrow = NA_INTEGER;
  int ncol= NA_INTEGER;

  // copy data for kd-tree
  nrow = data.nrow();
  ncol = data.ncol();
  dataPts = annAllocPts(nrow, ncol);
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  if (type==1) kdTree = new ANNkd_tree(dataPts, nrow, ncol, bucketSize,
    (ANNsplitRule) splitRule);
  else kdTree = new ANNbruteForce(dataPts, nrow, ncol);
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  std::vector<int> N;
  IntegerVector count(nrow);

  for (int i=0; i<nrow; ++i) {
    //Rprintf("processing point %d\n", i+1);
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    N = regionQuery(i, dataPts, kdTree, eps2, approx);
    count[i] = N.size();
  }

  // cleanup
  if (kdTree != NULL) delete kdTree;
  if (dataPts != NULL)  annDeallocPts(dataPts);
  // annClose(); is now done globally in the package

  return count;
}

