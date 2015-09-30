//----------------------------------------------------------------------
//                   Fixed Radius Nearest Neighbors
// File:                        frNN.cpp
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
#include "ANN/ANN.h"
#include "R_regionQuery.h"

using namespace Rcpp;

// [[Rcpp::export]]
List frNN_int(NumericMatrix data, double eps, int type, int bucketSize, int splitRule, double approx) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;

  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  ANNpointSet* kdTree = NULL;
  if (type==1){
    kdTree = new ANNkd_tree(dataPts, nrow, ncol, bucketSize,
      (ANNsplitRule)  splitRule);
  } else{
    kdTree = new ANNbruteForce(dataPts, nrow, ncol);
  }
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  // frNN
  std::vector< IntegerVector > id; id.resize(nrow);
  std::vector< NumericVector > dist; dist.resize(nrow);

  for (int p=0; p<nrow; p++) {
    if (!(p % 100)) Rcpp::checkUserInterrupt();

    //Rprintf("processing point %d\n", p+1);
    nn N = regionQueryDist(p, dataPts, kdTree, eps2, approx);

    // fix index
    //std::transform(N.first.begin(), N.first.end(),
    //  N.first.begin(), std::bind2nd( std::plus<int>(), 1 ) );

    // take sqrt of distance since the tree stores d^2
    //std::transform(N.second.begin(), N.second.end(),
    //  N.second.begin(), static_cast<double (*)(double)>(std::sqrt));

    // remove self matches
    IntegerVector ids = IntegerVector(N.first.begin(), N.first.end());
    LogicalVector take = ids != p;
    ids = ids[take];
    NumericVector dists = NumericVector(N.second.begin(), N.second.end())[take];

    id[p] = ids+1;
    dist[p] = sqrt(dists);




  }

  // cleanup
  delete kdTree;
  annDeallocPts(dataPts);
  annClose();

  // prepare results
  List ret;
  ret["id"] = id;
  ret["dist"] = dist;
  ret["eps"] = eps;
  return ret;
}
