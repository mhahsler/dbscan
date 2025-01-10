//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// Header
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
// this implementation is based on
// https://www.geeksforgeeks.org/prims-algorithm-in-cpp/
// [[Rcpp::export]]
Rcpp::NumericMatrix mst(const NumericVector x_dist, const R_xlen_t n) {
  Rcpp::NumericMatrix mst = NumericMatrix(n - 1, 3);
  colnames(mst) = CharacterVector::create("from", "to", "weight");

  // vector to store the parent of vertex
  std::vector<int> parent(n);
  std::vector<double> weight(n);
  std::vector<bool> visited(n);

  std::priority_queue<std::pair<double, int>,
                 std::vector<std::pair<double, int>>,
                 std::greater<std::pair<double, int>>> pq;

  // Initialize weights and visited
  for (int i = 0; i < n; i++) {
    weight[i] = INT_MAX;
    visited[i] = false;
  }

  // First node is always the root of MST. pick it first.
  parent[0] = -1;
  weight[0] = -INFINITY;

  // Push the source vertex to the min-heap
  pq.push({0, 0});

  while (!pq.empty()) {
    int node = pq.top().second;
    pq.pop();
    visited[node] = true;
    for (int i = 0; i < n; i++) {
      double the_weight = x_dist[LT_POS0(n, node, i)];
      if (!visited[i] && node != i
            && the_weight < weight[i]) {
        pq.push({the_weight, i});
        weight[i] = the_weight;
        parent[i] = node;
      }
    }
  }

  for (int i = 1; i < n; i++) {
    mst(i-1, 1) = parent[i] +1;
    mst(i-1, 0) = i + 1;
    mst(i-1, 2) = x_dist[LT_POS0(n, i, parent[i])];
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
    _["dist.method"] = "mutual reachability"
  );
  res.attr("class") = "hclust";
  return res;
}
