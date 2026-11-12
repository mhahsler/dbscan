//----------------------------------------------------------------------
//                                OPTICS
// File:                        R_optics.cpp
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3 or later
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
#include "ann_raii.h"
#include "ANN/ANN.h"
#include "regionQuery.h"
#include "lt.h"

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
      pre[o] = p;
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

  // k-d tree uses squared distances
  double eps2 = eps*eps;

  ANNPoints dataPts;
  std::unique_ptr<ANNpointSet> kdTree;
  int nrow = NA_INTEGER;
  int ncol= NA_INTEGER;

  if(frNN.size()) {
    // no k-d tree
    nrow = (as<List>(frNN["id"])).size();
  }else{

    // copy data for k-d tree
    nrow = data.nrow();
    ncol = data.ncol();
    dataPts = ANNPoints(nrow, ncol);
    for (int i = 0; i < nrow; i++){
      for (int j = 0; j < ncol; j++){
        (dataPts[i])[j] = data(i, j);
      }
    }
    //Rprintf("Points copied.\n");

    // create k-d tree (1) or linear search structure (2)
    kdTree = make_ann_tree(dataPts.get(), nrow, ncol, type,
      bucketSize, splitRule);
    //Rprintf("k-d tree ready. starting OPTICS.\n");

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
    if (!(p % 10)) Rcpp::checkUserInterrupt();
    //Rprintf("processing point %d\n", p+1);

    if (visited[p]) continue;

    // ExpandClusterOrder
    //N = regionQueryDist(p, dataPts.get(), kdTree.get(), eps2, approx);
    if(frNN.size())   N = std::make_pair(
      as<std::vector<int> >(as<List>(frNN["id"])[p]),
      as<std::vector<double> >(as<List>(frNN["dist"])[p]));
    else              N = regionQueryDist(p, dataPts.get(), kdTree.get(), eps2, approx);

    visited[p] = true;

    // find core distance
    if(N.second.size() >= (size_t) minPts) {
      ds = N.second;
      std::sort(ds.begin(), ds.end()); // sort increasing
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
        // always broken consistently (lower ID wins to produce the same
        // results as the current elki implementation)!
        if (reachdist[*it] < reachdist[*q_it] ||
          (reachdist[*it] == reachdist[*q_it] && *it < *q_it)) q_it = it;
      }
      q = *q_it;
      seeds.erase(q_it);

      //N2 = regionQueryDist(q, dataPts.get(), kdTree.get(), eps2, approx);
      if(frNN.size())   N = std::make_pair(
        as<std::vector<int> >(as<List>(frNN["id"])[q]),
        as<std::vector<double> >(as<List>(frNN["dist"])[q]));
      else              N = regionQueryDist(q, dataPts.get(), kdTree.get(), eps2, approx);

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

  // prepare results (R index starts with 1)
  List ret;
  ret["order"] = IntegerVector(orderedPoints.begin(), orderedPoints.end()) + 1;
  ret["reachdist"] = sqrt(NumericVector(reachdist.begin(), reachdist.end()));
  ret["coredist"] = sqrt(NumericVector(coredist.begin(), coredist.end()));
  ret["predecessor"] = IntegerVector(pre.begin(), pre.end()) + 1;
  return ret;
}


// ---------------------------------------------------------------------------
// OPTICS for eps = Inf
//
// For an unrestricted neighborhood, OPTICS does not need neighborhood queries.
// Every point is a core point and the reachability distance of an unvisited
// point o from the visited core point p is max(coredist(p), d(p, o)). OPTICS
// repeatedly visits the unvisited point with the smallest reachability
// distance (ties: lowest id), which is Prim's algorithm with these weights
// (starting with the first point). The predecessor is the first visited point
// that gave the smallest reachability distance.
//
// The distances are computed on the fly, so memory use is O(n) (+ O(n*d)
// for a data matrix) and time is O(n^2 * d). The original implementation
// optics_int() is kept for finite eps.
//
// coredist: core distances for all points (finite).
// dist_fn(i, j): distance between points i and j (0-based, not squared).
template <typename DistFun>
static List optics_prim_impl(const R_xlen_t n, const NumericVector &coredist,
                             DistFun dist_fn) {
  if (coredist.length() != n)
    stop("coredist needs to have one value for each point.");
  for (R_xlen_t i = 0; i < n; ++i)
    if (!R_FINITE(coredist[i]))
      stop("optics_prim needs finite core distances.");

  std::vector<int> order; order.reserve(n);
  std::vector<int> pre(n, NA_INTEGER);
  std::vector<double> reach(n, INFINITY);   // reachability distance (Inf = undefined)
  std::vector<double> weight(n, INFINITY);  // best reachability so far
  std::vector<bool> visited(n, false);

  int cur = 0;
  for (R_xlen_t k = 0; k < n; ++k) {
    if (!(k % 256)) Rcpp::checkUserInterrupt();

    order.push_back(cur);
    visited[cur] = true;
    if (k == n - 1) break;

    const double cd = coredist[cur];
    int next = -1;
    double next_weight = INFINITY;

    for (int i = 0; i < n; ++i) {
      if (visited[i]) continue;

      double w = std::max(cd, dist_fn(cur, i));
      if (w < weight[i]) {
        weight[i] = w;
        pre[i] = cur;
      }
      // strict '<' in ascending order of i: the lowest id wins ties
      if (weight[i] < next_weight) {
        next_weight = weight[i];
        next = i;
      }
    }

    if (next < 0) stop("optics_prim could not find the next point (non-finite distances?).");
    reach[next] = next_weight;
    cur = next;
  }

  IntegerVector order_r(order.begin(), order.end());
  IntegerVector pre_r(n);
  for (R_xlen_t i = 0; i < n; ++i)
    pre_r[i] = (pre[i] == NA_INTEGER) ? NA_INTEGER : pre[i] + 1;

  List ret;
  ret["order"] = order_r + 1;  // R index starts with 1
  ret["reachdist"] = NumericVector(reach.begin(), reach.end());
  ret["coredist"] = clone(coredist);
  ret["predecessor"] = pre_r;
  return ret;
}

// Data matrix: Euclidean distances are computed on the fly.
// [[Rcpp::export]]
List optics_prim(const NumericMatrix x, const NumericVector coredist) {
  const R_xlen_t n = x.nrow();
  const R_xlen_t d = x.ncol();
  if (n < 1) stop("x needs to contain at least 1 object.");

  // row-major copy of the data for cache-friendly access to a point
  std::vector<double> pts(n * d);
  for (R_xlen_t j = 0; j < d; ++j)
    for (R_xlen_t i = 0; i < n; ++i)
      pts[i * d + j] = x(i, j);

  return optics_prim_impl(n, coredist, [&](int a, int b) {
    const double *pa = &pts[(R_xlen_t) a * d];
    const double *pb = &pts[(R_xlen_t) b * d];
    double sum = 0.0;
    for (R_xlen_t j = 0; j < d; ++j) {
      double diff = pa[j] - pb[j];
      sum += diff * diff;
    }
    return std::sqrt(sum);
  });
}

// dist object (lower triangle vector) with n points.
// [[Rcpp::export]]
List optics_prim_dist(const NumericVector x_dist, const R_xlen_t n,
                      const NumericVector coredist) {
  if (n < 1) stop("x needs to contain at least 1 object.");
  if (x_dist.length() != n * (n - 1) / 2)
    stop("the length of the dist object does not match n.");

  return optics_prim_impl(n, coredist, [&](int a, int b) {
    return (double) x_dist[LT_POS0((R_xlen_t) n, (R_xlen_t) a, (R_xlen_t) b)];
  });
}
