//----------------------------------------------------------------------
//                                OPTICS
// File:                        optics.cpp
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

void update(
    std::pair< std::vector<int>, std::vector<double> > &N,
    int p,
    std::vector<int> &seeds,
    double eps2,
    int minPts,
    std::vector <bool> &visited,
    std::vector<int> &orderedPoints,
    std::vector<double> &reachdist,
    std::vector<double> &coredist
    ){

  std::vector<int>::iterator pos_seeds;

  // find core distance
  // Note: we hopefully cannot get here if ds is not at least minPts long!
  std::vector<double> ds = N.second;
  std::sort(ds.begin(), ds.end());
  coredist[p] = ds[minPts-1];

  while(!N.first.empty()) {
    int o = N.first.back();
    double o_d = N.second.back();
    N.first.pop_back();
    N.second.pop_back();

    if(visited[o]) continue;

    double newreachdist = std::max(coredist[p], o_d);

    pos_seeds = std::find(seeds.begin(), seeds.end(), o);
    if(pos_seeds == seeds.end()) {
      reachdist[o] = newreachdist;
      seeds.push_back(o);
    } else if(newreachdist < reachdist[o]) reachdist[o] = newreachdist;
  }
}


// [[Rcpp::export]]
List optics_int(NumericMatrix data, double eps, int minPts,
  int type, int bucketSize, int splitRule, double approx) {

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

  // OPTICS
  std::vector<bool> visited(nrow, false);
  std::vector<int> orderedPoints; orderedPoints.reserve(nrow);
  std::vector<double> reachdist(nrow, INFINITY); // we used Inf as undefined
  std::vector<double> coredist(nrow, INFINITY);
  nn N, N2;
  std::vector<int> seeds;

  for (int p=0; p<nrow; p++) {
    if (!(p % 100)) Rcpp::checkUserInterrupt();
    //Rprintf("processing point %d\n", p+1);

    if (visited[p]) continue;

    // ExpandClusterOrder
    // Note:
    N = regionQueryDist(p, dataPts, kdTree, eps2, approx);

    // find core distance
    if(N.second.size() >= (size_t) minPts) {
      std::vector<double> ds = N.second;
      std::sort(ds.begin(), ds.end()); // sort inceasing
      coredist[p] = ds[minPts-1];
    }

    // mark visited and output p
    visited[p] = true;
    orderedPoints.push_back(p);

    if (coredist[p] == INFINITY) continue; // core-dist is undefined

    // updateable priority queue does not exist in C++ STL!
    seeds.clear();

    // update
    update(N, p, seeds, eps2, minPts, visited, orderedPoints,
      reachdist, coredist);

    while (!seeds.empty()) {
      // find smallest dist
      std::vector<int>::iterator q_it = seeds.begin();
      for (std::vector<int>::iterator it = seeds.begin();
        it!=seeds.end(); ++it) {
        if (reachdist[*it] < reachdist[*q_it]) q_it = it;
      }

      int q = *q_it;
      seeds.erase(q_it);

      visited[q] = true;
      // FIXME: set core dist?

      orderedPoints.push_back(q);

      N2 = regionQueryDist(q, dataPts, kdTree, eps2, approx);

      // contains q?
      if(N2.first.size() < (size_t) minPts-1) continue; // q has no core dist.
      update(N2, q, seeds, eps2, minPts, visited, orderedPoints,
        reachdist, coredist);
    }
  }

  // cleanup
  delete kdTree;
  annDeallocPts(dataPts);
  annClose();

  // prepare results (R index starts with 1)
  List ret;
  ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end())+1;
  ret["reachdist"] = sqrt(NumericVector(reachdist.begin(), reachdist.end()));
  ret["coredist"] = sqrt(NumericVector(coredist.begin(), coredist.end()));
  return ret;
}

