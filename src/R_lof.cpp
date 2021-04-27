//----------------------------------------------------------------------
//                  Find the Neighbourhood for LOF
// File:                    R_lof.cpp
//----------------------------------------------------------------------
// Copyright (c) 2021 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// LOF needs to find the k-NN distance and then how many points are within this
// neighborhood.

#include <Rcpp.h>
#include "R_regionQuery.h"
//#include "ANN/ANN.h"

using namespace Rcpp;

// returns knn-dist and the neighborhood size as a matrix
// [[Rcpp::export]]
List lof_kNN(NumericMatrix data, int minPts,
  int type, int bucketSize, int splitRule, double approx) {

  // minPts includes the point itself; k does not!
  int k = minPts - 1;

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

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];
  nn N;

  // results
  List id(nrow);
  List dist(nrow);
  NumericVector k_dist(nrow);

  for (int i=0; i<nrow; i++) {
    //Rprintf("processing point %d\n", p+1);
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = dataPts[i];

    // find k-NN distance
    kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    k_dist[i] = ANN_ROOT(dists[k]); // this is a squared distance!

    // find k-NN neighborhood which can be larger than k with tied distances
    // This works under Linux and Windows, but not under Solaris: The points at the
    // k_distance may not be included.
    //nn N = regionQueryDist_point(queryPt, dataPts, kdTree, dists[k], approx);

    // Make the comparison robust.
    // Compare doubles: http://c-faq.com/fp/fpequal.html
    double minPts_dist = dists[k] + DBL_EPSILON * dists[k];
    nn N = regionQueryDist_point(queryPt, dataPts, kdTree, minPts_dist, approx);

    IntegerVector ids = IntegerVector(N.first.begin(), N.first.end());
    NumericVector dists = NumericVector(N.second.begin(), N.second.end());

    // remove self matches -- not an issue with query points
    LogicalVector take = ids != i;
    ids = ids[take];
    dists = dists[take];

    id[i] = ids+1;
    dist[i] = sqrt(dists);
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  // annClose(); is now done globally in the package

  // all k_dists are squared
  //k_dist = sqrt(k_dist);

  // prepare results
  List ret;
  ret["k_dist"] = k_dist;
  ret["ids"] = id;
  ret["dist"] = dist;
  return ret;
}
