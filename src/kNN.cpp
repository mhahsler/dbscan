//----------------------------------------------------------------------
//                  Find the k Nearest Neighbors
// File:                    R_kNNdist.cpp
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3 or later
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// Note: does not return self-matches!

#include "kNN.h"
#include "ann_raii.h"

// returns knn + dist
List kNN_int(NumericMatrix data, int k,
  int type, int bucketSize, int splitRule, double approx) {

  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNPoints dataPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  auto kdTree = make_ann_tree(dataPts.get(), nrow, ncol, type,
    bucketSize, splitRule);
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  NumericMatrix d(nrow, k);
  IntegerMatrix id(nrow, k);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  std::vector<ANNdist> dists(k+1);
  std::vector<ANNidx> nnIdx(k+1);

  for (int i=0; i<nrow; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = dataPts[i];

    kdTree->annkSearch(queryPt, k+1, nnIdx.data(), dists.data(), approx);

    // remove self match; with duplicate points the query point may not be
    // among the k+1 hits, in which case we keep the k nearest
    IntegerVector ids = IntegerVector(nnIdx.begin(), nnIdx.end());
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = head(ids, k) + 1;

    NumericVector ndists = NumericVector(dists.begin(), dists.end())[take];
    d(i, _) = sqrt(head(ndists, k));
  }


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
  ANNPoints dataPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }

  // copy query
  int nrow_q = query.nrow();
  int ncol_q = query.ncol();
  ANNPoints queryPts(nrow_q, ncol_q);
  for(int i = 0; i < nrow_q; i++){
    for(int j = 0; j < ncol_q; j++){
      (queryPts[i])[j] = query(i, j);
    }
  }
  //Rprintf("Points copied.\n");

  // create kd-tree (1) or linear search structure (2)
  auto kdTree = make_ann_tree(dataPts.get(), nrow, ncol, type,
    bucketSize, splitRule);
  //Rprintf("kd-tree ready. starting DBSCAN.\n");

  NumericMatrix d(nrow_q, k);
  IntegerMatrix id(nrow_q, k);

  // Note: does not return itself with query
  std::vector<ANNdist> dists(k);
  std::vector<ANNidx> nnIdx(k);

  for (int i=0; i<nrow_q; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    ANNpoint queryPt = queryPts[i];
    kdTree->annkSearch(queryPt, k, nnIdx.data(), dists.data(), approx);

    IntegerVector ids = IntegerVector(nnIdx.begin(), nnIdx.end());
    id(i, _) = ids + 1;

    NumericVector ndists = NumericVector(dists.begin(), dists.end());
    d(i, _) = sqrt(ndists);
  }

  // prepare results (ANN returns points sorted by distance)
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  ret["sort"] = true;
  return ret;
}
