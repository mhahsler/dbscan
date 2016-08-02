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
    int minPts,
    std::vector <bool> &visited,
    std::vector<int> &orderedPoints,
    std::vector<double> &reachdist,
    std::vector<double> &coredist, 
    std::vector<int> &pre){

  std::vector<int>::iterator pos_seeds;
  double newreachdist;
  int o;
  double o_d;

  while(!N.first.empty()) {
    o = N.first.back();
    o_d = N.second.back();
    N.first.pop_back();
    N.second.pop_back();

    if(visited[o]) continue;

    newreachdist = std::max(coredist[p], o_d);

    if(reachdist[o] == INFINITY) {
      reachdist[o] = newreachdist;
      seeds.push_back(o);
    } else {
      // o was not visited and has a reachability distance must be
      // already in seeds!
      if(newreachdist < reachdist[o]) {
        reachdist[o] = newreachdist;
        pre[o] = p;
      }
    }
  }
}


// [[Rcpp::export]]
List optics_int(NumericMatrix data, double eps, int minPts,
  int type, int bucketSize, int splitRule, double approx, List frNN) {

  // kd-tree uses squared distances
  double eps2 = eps*eps;

  ANNpointSet* kdTree = NULL;
  ANNpointArray dataPts = NULL;
  int nrow = NA_INTEGER;
  int ncol= NA_INTEGER;

  if(frNN.size()) {
    // no kd-tree
    nrow = (as<List>(frNN["id"])).size();
  }else{

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
    //Rprintf("kd-tree ready. starting OPTICS.\n");

  }


  // OPTICS
  std::vector<bool> visited(nrow, false);
  std::vector<int> orderedPoints; orderedPoints.reserve(nrow);
  std::vector<int> pre(nrow, NA_INTEGER);
  std::vector<double> reachdist(nrow, INFINITY); // we used Inf as undefined
  std::vector<double> coredist(nrow, INFINITY);
  nn N;
  std::vector<int> seeds;
  std::vector<double> ds;

  for (int p=0; p<nrow; p++) {
    if (!(p % 100)) Rcpp::checkUserInterrupt();
    //Rprintf("processing point %d\n", p+1);

    if (visited[p]) continue;

    // ExpandClusterOrder
    //N = regionQueryDist(p, dataPts, kdTree, eps2, approx);
    if(frNN.size())   N = std::make_pair(
      as<std::vector<int> >(as<List>(frNN["id"])[p]),
      as<std::vector<double> >(as<List>(frNN["dist"])[p]));
    else              N = regionQueryDist(p, dataPts, kdTree, eps2, approx);

    visited[p] = true;

    // find core distance
    if(N.second.size() >= (size_t) minPts) {
      ds = N.second;
      std::sort(ds.begin(), ds.end()); // sort inceasing
      coredist[p] = ds[minPts-1];
    }
    int tmp_p = NA_INTEGER; 
    if (pre[p] == NA_INTEGER) { tmp_p = p; }
    orderedPoints.push_back(p);

    if (coredist[p] == INFINITY) continue; // core-dist is undefined

    // updateable priority queue does not exist in C++ STL so we use a vector!
    //seeds.clear();

    // update
    update(N, p, seeds, minPts, visited, orderedPoints,
      reachdist, coredist, pre);

    int q;
    while (!seeds.empty()) {
      // get smallest dist (to emulate priority queue). All should have already
      // a reachability distance <Inf from update().
      std::vector<int>::iterator q_it = seeds.begin();
      for (std::vector<int>::iterator it = seeds.begin();
        it!=seeds.end(); ++it) {
        // Note: The second part of the if statement ensures that ties are
        // always broken consistenty (higher ID wins to produce the same
        // results as the elki implementation)!
        if (reachdist[*it] < reachdist[*q_it] ||
          (reachdist[*it] == reachdist[*q_it] && *q_it < *it)) q_it = it;
      }
      q = *q_it;
      seeds.erase(q_it);

      //N2 = regionQueryDist(q, dataPts, kdTree, eps2, approx);
      if(frNN.size())   N = std::make_pair(
        as<std::vector<int> >(as<List>(frNN["id"])[q]),
        as<std::vector<double> >(as<List>(frNN["dist"])[q]));
      else              N = regionQueryDist(q, dataPts, kdTree, eps2, approx);

      visited[q] = true;

      // update core distance
      if(N.second.size() >= (size_t) minPts) {
        ds = N.second;
        std::sort(ds.begin(), ds.end());
        coredist[q] = ds[minPts - 1];
      }
      if (pre[q] == NA_INTEGER) { pre[q] = tmp_p; }
      orderedPoints.push_back(q);

      if(N.first.size() < (size_t) minPts) continue; //  == q has no core dist. 

      // update seeds
      update(N, q, seeds, minPts, visited, orderedPoints,
        reachdist, coredist, pre);
    }
  }

  // cleanup
  if(kdTree != NULL) {
    delete kdTree;
    annDeallocPts(dataPts);
    annClose();
  }

  // prepare results (R index starts with 1)
  List ret;
  ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) + 1;
  ret["reachdist"] = sqrt(NumericVector(reachdist.begin(), reachdist.end()));
  ret["coredist"] = sqrt(NumericVector(coredist.begin(), coredist.end()));
  ret["predecessor"] = IntegerVector(pre.begin(), pre.end()) + 1;
  return ret;
}

