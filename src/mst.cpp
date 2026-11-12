//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3 or later
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#include "mst.h"

// coreFromDist indexes through the a dist vector to retrieve the core distance;
// this might be useful in some situations. For example, you can get the core distance
// from only a dist object, without needing the original data. In experimentation, the
// kNNdist ended up being faster than this.
//
// // [[Rcpp::export]]
// NumericVector coreFromDist(const NumericVector dist, const int n, const int minPts){
//   NumericVector core_dist = NumericVector(n);
//   NumericVector row_dist = NumericVector(n - 1);
//   for (R_xlen_t i = 0; i < n; ++i){
//     for (R_xlen_t j = 0; j < n; ++j){
//       if (i == j) continue;
//       R_xlen_t index = LT_POS0(n, j, i)
//       row_dist.at(j > i ? j  - 1 : j) = dist.at(index);
//     }
//     std::sort(row_dist.begin(), row_dist.end());
//     core_dist[i] = row_dist.at(minPts-2); // one for 0-based indexes, one for inclusive minPts condition
//   }
//   return(core_dist);
// }


// Prim's Algorithm
// this implementation for dense dist objects avoids the use of a min-heap.
// [[Rcpp::export]]
Rcpp::NumericMatrix mst_dist_cpp(const NumericVector x_dist, const R_xlen_t n) {
  Rcpp::NumericMatrix mst = NumericMatrix(n - 1, 3);
  colnames(mst) = CharacterVector::create("from", "to", "weight");

  // vector to store the parent of vertex
  std::vector<int> parent(n);
  std::vector<double> weight(n, INFINITY);
  std::vector<bool> visited(n, false);

  // first node is always the root of MST.
  parent[0] = -1;
  weight[0] = 0;

  int next_node = 0;
  double next_weight;
  int node;

  while (next_node >= 0) {
    node = next_node;
    next_node = -1;
    next_weight = INFINITY;

    visited[node] = true;
    if (node > 0) {
      mst(node-1, 1) = parent[node] +1;
      mst(node-1, 0) = node + 1;
      mst(node-1, 2) = weight[node];
    }

    for (int i = 1; i < n; i++) { // 0 is always the first node
      if (visited[i] || node == i) continue;

      double the_weight = x_dist[LT_POS0(n, node, i)];
      if (the_weight < weight[i]) {
        weight[i] = the_weight;
        parent[i] = node;
      }

      // find minimum weight node
      if (weight[i] < next_weight) {
        next_weight = weight[i];
        next_node = i;
      }

    }
  }

  return(mst);
}

// Prim's Algorithm on the mutual reachability graph computed from the data.
// Same algorithm (and tie-breaking) as mst_dist_cpp(), but the Euclidean distances are
// computed on the fly and mrd(a, b) = max(coredist(a), coredist(b), d(a, b)) is
// evaluated when needed. No distance matrix is stored, so memory use is O(n * d)
// instead of O(n^2).
//
// Input:
// * x: numeric data matrix (n x d) with finite values
// * coredist: core distances (length n)
//
// Returns: the MST as a matrix with columns from, to, weight (as in mst_dist_cpp()).
// [[Rcpp::export]]
Rcpp::NumericMatrix mst_mrd(const NumericMatrix x, const NumericVector coredist) {
  const R_xlen_t n = x.nrow();
  const R_xlen_t d = x.ncol();

  if (n < 2)
    stop("x needs to contain at least 2 objects.");
  if (coredist.length() != n)
    stop("coredist needs to have one value for each row of x.");

  // row-major copy of the data (O(n*d)) for cache-friendly access to a point
  std::vector<double> pts(n * d);
  for (R_xlen_t j = 0; j < d; ++j)
    for (R_xlen_t i = 0; i < n; ++i)
      pts[i * d + j] = x(i, j);

  Rcpp::NumericMatrix mst = NumericMatrix(n - 1, 3);
  colnames(mst) = CharacterVector::create("from", "to", "weight");

  std::vector<int> parent(n);
  std::vector<double> weight(n, INFINITY);
  std::vector<bool> visited(n, false);

  parent[0] = -1;
  weight[0] = 0;

  int next_node = 0;
  double next_weight;
  int node;

  while (next_node >= 0) {
    node = next_node;
    next_node = -1;
    next_weight = INFINITY;

    visited[node] = true;
    if (node > 0) {
      mst(node - 1, 1) = parent[node] + 1;
      mst(node - 1, 0) = node + 1;
      mst(node - 1, 2) = weight[node];
    }

    const double *pn = &pts[node * d];
    const double cn = coredist[node];

    for (int i = 1; i < n; i++) { // 0 is always the first node
      if (visited[i] || node == i) continue;

      // Euclidean distance (same accumulation order as stats::dist)
      const double *pi = &pts[i * d];
      double sum = 0.0;
      for (R_xlen_t j = 0; j < d; ++j) {
        double diff = pn[j] - pi[j];
        sum += diff * diff;
      }
      double the_weight = std::max(std::sqrt(sum), std::max(cn, (double)coredist[i]));

      if (the_weight < weight[i]) {
        weight[i] = the_weight;
        parent[i] = node;
      }

      if (weight[i] < next_weight) {
        next_weight = weight[i];
        next_node = i;
      }
    }

    if ((node & 0xFFF) == 0) checkUserInterrupt();
  }

  return(mst);
}

//
// // [[Rcpp::export]]
// IntegerVector order_(NumericVector x) {
//   if (is_true(any(duplicated(x)))) {
//     Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
//   }
//   NumericVector sorted = clone(x).sort();
//   return match(sorted, x);
// }


// Single link hierarchical clustering
// used by GLOSH.R and hdbscan.R

void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind) {
  // base case
  if (merge(i, j) < 0) {
    order.at(ind++) = -merge(i, j);
  }
  else {
    visit(merge, order, merge(i, j) - 1, 0, ind);
    visit(merge, order, merge(i, j) - 1, 1, ind);
  }
}

IntegerVector extractOrder(IntegerMatrix merge){
  IntegerVector order = IntegerVector(merge.nrow()+1);
  int ind = 0;
  visit(merge, order, merge.nrow() - 1, 0, ind);
  visit(merge, order, merge.nrow() - 1, 1, ind);
  return(order);
}

// [[Rcpp::export]]
List hclustMergeOrder(NumericMatrix mst, IntegerVector o){
  int npoints = mst.nrow() + 1;
  NumericVector dist = mst(_, 2);

  // Extract order, reorder indices
  NumericVector left = mst(_, 0), right = mst(_, 1);
  IntegerVector left_int = as<IntegerVector>(left[o-1]), right_int = as<IntegerVector>(right[o-1]);

  // Labels and resulting merge matrix
  IntegerVector labs = -seq_len(npoints);
  IntegerMatrix merge = IntegerMatrix(npoints - 1, 2);

  // Replace singletons as negative and record merge of non-singletons as positive
  for (int i = 0; i < npoints - 1; ++i) {
    int lab_left = labs.at(left_int.at(i)-1), lab_right = labs.at(right_int.at(i)-1);
    merge(i, _) = IntegerVector::create(lab_left, lab_right);
    for (int c = 0; c < npoints; ++c){
      if (labs.at(c) == lab_left || labs.at(c) == lab_right){
        labs.at(c) = i+1;
      }
    }
  }
  //IntegerVector int_labels = seq_len(npoints);
  List res = List::create(
    _["merge"] = merge,
    _["height"] = dist[o-1],
    _["order"] = extractOrder(merge),
    _["labels"] = R_NilValue, //as<StringVector>(int_labels)
    _["method"] = "robust single",
    _["dist.method"] = "unknown"
  );
  res.attr("class") = "hclust";
  return res;
}
