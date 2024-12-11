//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// Header 
#include "prims_mst.h"

// inline indexing of distance matrix 
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

// Structures to do priority queue
struct edge
{
  int to;
  double weight;
  edge(int to_id, double cost) : to(to_id), weight(cost) { }
};

struct compare_edge
{
  bool operator()(const edge& e1, const edge& e2) const
  { return e1.weight > e2.weight; }
};

// For some reason, STL nor standards < C++11 give reasonable means of reserving container
// memory for the priority queue, so make own!
template <class T, class Container = std::vector<T>,
          class Compare = std::less<typename Container::value_type> >
class r_priority_queue : public std::priority_queue<T, Container , Compare>
{
public:
  r_priority_queue(size_t size)
  {
    this->c.reserve(size);
  }
};

// coreFromDist indexes through the a dist vector to retrieve the core distance;
// this might be useful in some situations. For example, you can get the core distance
// from only a dist object, without needing the original data. In experimentation, the
// kNNdist ended up being faster than this.

// [[Rcpp::export]]
NumericVector coreFromDist(const NumericVector dist, const int n, const int minPts){
  NumericVector core_dist = NumericVector(n);
  NumericVector row_dist = NumericVector(n - 1);
  for (R_xlen_t i = 0; i < n; ++i){
    for (R_xlen_t j = 0; j < n; ++j){
      if (i == j) continue;
      R_xlen_t index = i > j ? INDEX_TF(n, j, i) : INDEX_TF(n, i, j);
      row_dist.at(j > i ? j  - 1 : j) = dist.at(index);
    }
    std::sort(row_dist.begin(), row_dist.end());
    core_dist[i] = row_dist.at(minPts-2); // one for 0-based indexes, one for inclusive minPts condition
  }
  return(core_dist);
}

// [[Rcpp::export]]
NumericMatrix prims(const NumericVector x_dist, const R_xlen_t n) {

  // Resulting MST
  NumericMatrix mst = NumericMatrix(n-1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));
  // r_priority_queue<edge, std::vector<edge>, compare_edge> _heap(n);
  // std::priority_queue<edge, std::vector<edge>, compare_edge>* heap = &_heap;

  NumericVector row_entry;
  double min = std::numeric_limits<double>::infinity(), priority = 0.0;
  R_xlen_t c_node = 0;
  int min_id = n-1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen
    // Compare all the new edge weight w/ the "current best" edge weights
    for (R_xlen_t i = 0; i < n; ++i) {
      if (i == c_node) continue;
      if (v_selected[i] < 0) {
        R_xlen_t index = i > c_node ? INDEX_TF(n, c_node, i) : INDEX_TF(n, i, c_node); // bigger index always on the right
        priority = x_dist[index]; // get distance
        if (priority < fringe[i].weight) {
          fringe[i].weight = priority;
          fringe[i].to = c_node; // i indexes the 'from' node
        }
        if (fringe[i].weight < min) { // an edge 'on the fringe' might be less than any of the current nodes weights
          min = fringe[i].weight;
          min_id = i;
        }
      }
    }

    // Extract and insert the minimum edge
    mst(n_edges, _) = Rcpp::NumericVector::create(min_id+1, c_node+1, min);
    v_selected[c_node] = 1;
    c_node = min_id;
  }

  return(mst);
}


// [[Rcpp::export]]
IntegerVector order_(NumericVector x) {
  if (is_true(any(duplicated(x)))) {
    Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}


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
