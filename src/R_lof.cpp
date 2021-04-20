//----------------------------------------------------------------------
//                  Find the Neighbourhood for LOF
// File:                    R_lof.cpp
//----------------------------------------------------------------------
// Copyright (c) 2021 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// LOF needs to find the k-NN distance and then how many points are withing this
// neighbourhood.

#include <Rcpp.h>
#include "R_regionQuery.h"
//#include "ANN/ANN.h"

using namespace Rcpp;

// returns knn-dist and the neighborhood size as a matrix
// [[Rcpp::export]]
List lof_kNN(NumericMatrix data, int k,
  int type, int bucketSize, int splitRule, double approx) {

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
  ANNpoint queryPt;

  // results
  List id(nrow);
  List dist(nrow);
  NumericVector k_dist(nrow);

  for (int i=0; i<nrow; i++) {
    //Rprintf("processing point %d\n", p+1);
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    queryPt = dataPts[i];

    //if(type==1) kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    //else kdTree->annkSearch(queryPt, k+1, nnIdx, dists);

    // find k-NN distance
    kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    k_dist[i] = dists[k]; // this is squared

    // find k-NN neighborhood which can be larger than k with tied distances
    nn N = regionQueryDist_point(queryPt, dataPts, kdTree, k_dist[i] , approx);

    IntegerVector ids = IntegerVector(N.first.begin(), N.first.end());
    NumericVector dists = NumericVector(N.second.begin(), N.second.end());

    // remove self matches -- not an issue with query points
    LogicalVector take = ids != i;
    ids = ids[take];
    dists = dists[take];

    id[i] = ids+1;
    dist[i] = sqrt(dists);
  }

  k_dist = sqrt(k_dist);

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  // annClose(); is now done globally in the package


  // prepare results
  List ret;
  ret["dist"] = dist;
  ret["ids"] = id;
  ret["k_dist"] = k_dist;
  return ret;
}
