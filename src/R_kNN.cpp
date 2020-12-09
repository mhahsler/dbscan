//----------------------------------------------------------------------
//                  Find the k Nearest Neighbors
// File:                    kNNdist.cpp
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// Note: does not return self-matches!

#include <Rcpp.h>
#include "ANN/ANN.h"

using namespace Rcpp;

// returns knn + dist
// [[Rcpp::export]]
List kNN_int(NumericMatrix data, int k,
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

  NumericMatrix d(nrow, k);
  IntegerMatrix id(nrow, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i<nrow; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = dataPts[i];

    if(type==1) kdTree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    else kdTree->annkSearch(queryPt, k+1, nnIdx, dists);

    // remove self match
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;

    NumericVector ndists = NumericVector(dists, dists+k+1)[take];
    d(i, _) = sqrt(ndists);
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  // annClose(); is now done globaly in the package


  // prepare results
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  ret["sort"] = true;
  return ret;
}

// returns knn + dist using data and query
// [[Rcpp::export]]
List kNN_query_int(NumericMatrix data, NumericMatrix query, int k,
  int type, int bucketSize, int splitRule, double approx) {

  // FIXME: check ncol for data and query

  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }

  // copy query
  int nrow_q = query.nrow();
  int ncol_q = query.ncol();
  ANNpointArray queryPts = annAllocPts(nrow_q, ncol_q);
  for(int i = 0; i < nrow_q; i++){
    for(int j = 0; j < ncol_q; j++){
      (queryPts[i])[j] = query(i, j);
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

  NumericMatrix d(nrow_q, k);
  IntegerMatrix id(nrow_q, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  for (int i=0; i<nrow_q; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = queryPts[i];

    if(type==1) kdTree->annkSearch(queryPt, k, nnIdx, dists, approx);
    else kdTree->annkSearch(queryPt, k, nnIdx, dists);

    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k);
    id(i, _) = ids + 1;

    NumericVector ndists = NumericVector(dists, dists+k);
    d(i, _) = sqrt(ndists);
  }

  // cleanup
  delete kdTree;
  delete [] dists;
  delete [] nnIdx;
  annDeallocPts(dataPts);
  annDeallocPts(queryPts);
  // annClose(); is now done globaly in the package

  // prepare results (ANN returns points sorted by distance)
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  ret["sort"] = true;
  return ret;
}
