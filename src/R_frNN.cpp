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
List frNN_int(NumericMatrix data, double eps, int type,
  int bucketSize, int splitRule, double approx) {

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
  //std::vector< IntegerVector > id; id.resize(nrow);
  //std::vector< NumericVector > dist; dist.resize(nrow);
  List id(nrow);
  List dist(nrow);

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

    IntegerVector ids = IntegerVector(N.first.begin(), N.first.end());
    NumericVector dists = NumericVector(N.second.begin(), N.second.end());

    // remove self matches
    LogicalVector take = ids != p;
    ids = ids[take];
    dists = dists[take];

    //Rprintf("Found neighborhood size %d\n", ids.size());
    id[p] = ids+1;
    dist[p] = sqrt(dists);
  }

  // cleanup
  delete kdTree;
  annDeallocPts(dataPts);
  // annClose(); is now done globally in the package

  // prepare results
  List ret;
  ret["dist"] = dist;
  ret["id"] = id;
  ret["eps"] = eps;
  return ret;
}

// [[Rcpp::export]]
List frNN_query_int(NumericMatrix data, NumericMatrix query, double eps, int type,
  int bucketSize, int splitRule, double approx) {

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

  // frNN
  //std::vector< IntegerVector > id; id.resize(nrow);
  //std::vector< NumericVector > dist; dist.resize(nrow);
  List id(nrow_q);
  List dist(nrow_q);

  for (int p=0; p<nrow_q; p++) {
    if (!(p % 100)) Rcpp::checkUserInterrupt();

    //Rprintf("processing point %d\n", p+1);
    ANNpoint queryPt = queryPts[p];
    nn N = regionQueryDist_point(queryPt, dataPts, kdTree, eps2, approx);

    // fix index
    //std::transform(N.first.begin(), N.first.end(),
    //  N.first.begin(), std::bind2nd( std::plus<int>(), 1 ) );

    // take sqrt of distance since the tree stores d^2
    //std::transform(N.second.begin(), N.second.end(),
    //  N.second.begin(), static_cast<double (*)(double)>(std::sqrt));

    IntegerVector ids = IntegerVector(N.first.begin(), N.first.end());
    NumericVector dists = NumericVector(N.second.begin(), N.second.end());

    // remove self matches -- not an issue with query points
    //LogicalVector take = ids != p;
    //ids = ids[take];
    //dists = dists[take];

    //Rprintf("Found neighborhood size %d\n", ids.size());
    id[p] = ids+1;
    dist[p] = sqrt(dists);
  }

  // cleanup
  delete kdTree;
  annDeallocPts(dataPts);
  annDeallocPts(queryPts);
  // annClose(); is now done globally in the package

  // prepare results
  List ret;
  ret["dist"] = dist;
  ret["id"] = id;
  ret["eps"] = eps;
  ret["sort"] = false;
  return ret;
}
