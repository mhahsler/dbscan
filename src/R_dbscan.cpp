//----------------------------------------------------------------------
//                                DBSCAN
// File:                        dbscan.cpp
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
IntegerVector dbscan_int(
    NumericMatrix data, double eps, int minPts, NumericVector weights,
    int borderPoints, int type, int bucketSize, int splitRule, double approx) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;
  bool weighted = FALSE;
  double Nweight = 0.0;

  if (weights.size() != 0) {
    if (weights.size() != data.nrow())
      stop("length of weights vector is incompatible with data.");
    weighted = TRUE;
  }

  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  ANNpointSet* kdTree = NULL;
  if (type==1)kdTree = new ANNkd_tree(dataPts, nrow, ncol, bucketSize,
    (ANNsplitRule) splitRule);
  else kdTree = new ANNbruteForce(dataPts, nrow, ncol);
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  // DBSCAN
  std::vector<bool> visited(nrow, false);
  std::vector< std::vector<int> > clusters; // vector of vectors == list

  for (int i=0; i<nrow; i++) {
    //Rprintf("processing point %d\n", i+1);
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    if (visited[i]) continue;


    std::vector<int> N = regionQuery(i, dataPts, kdTree, eps2, approx);

    // noise points stay unassigned for now
    //if (weighted) Nweight = sum(weights[IntegerVector(N.begin(), N.end())]) +
    if (weighted) {
      // This should work, but Rcpp has a problem with the sugar expression!
      // Assigning the subselection forces it to be materialized.
      // Nweight = sum(weights[IntegerVector(N.begin(), N.end())]) +
      // weights[i];
      NumericVector w = weights[IntegerVector(N.begin(), N.end())];
      Nweight = sum(w);
    } else Nweight = N.size();

    if (Nweight < minPts) continue;

    // start new cluster and expand
    std::vector<int> cluster;
    cluster.push_back(i);
    visited[i] = true;

    while (!N.empty()) {
      int j = N.back();
      N.pop_back();

      if (visited[j]) continue; // point already processed
      visited[j] = true;

      std::vector<int> N2 = regionQuery(j, dataPts, kdTree, eps2, approx);

      if (weighted) {
        // Nweight = sum(weights(NumericVector(N2.begin(), N2.end())) +
        // weights[j]
        NumericVector w = weights[IntegerVector(N2.begin(), N2.end())];
        Nweight = sum(w);
      } else Nweight = N2.size();

      if (Nweight >= minPts) { // expand neighborhood
        // this is faster than set_union and does not need sort! visited takes
        // care of duplicates.
        std::copy(N2.begin(), N2.end(),
          std::back_inserter(N));
      }

      // for DBSCAN* (borderPoints==FALSE) border points are considered noise
      if(Nweight >= minPts || borderPoints) cluster.push_back(j);
    }

    // add cluster to list
    clusters.push_back(cluster);
  }

  // prepare cluster vector
  // unassigned points are noise (cluster 0)
  IntegerVector id(nrow, 0);
  for (std::size_t i=0; i<clusters.size(); i++) {
    for (std::size_t j=0; j<clusters[i].size(); j++) {
      id[clusters[i][j]] = i+1;
    }
  }

  // cleanup
  delete kdTree;
  annDeallocPts(dataPts);
  annClose();

  return wrap(id);
}

