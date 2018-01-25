#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include "pr_queue_k.h"

#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)

// Given a 'dist' object and a integer k, return the k-th nearest neighbor for each point
// [[Rcpp::export]]
List knn_dist(NumericVector dist_x, const int k, int all) {
  const int n = as<int>(dist_x.attr("Size"));
  NumericMatrix knn_dist = NumericMatrix(n, all ? k : 1);
  IntegerMatrix knn_ids = IntegerMatrix(n, all ? k : 1);
  ANNmin_k* pr_queue = new ANNmin_k(k);
  double cdist = 0.0;
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      if (i != j){
        cdist = dist_x[INDEX_TF(n, min(i, j), max(i, j))];
        if (cdist < pr_queue->max_key()){
          pr_queue->insert(cdist, j);
        }
      }
    }
    if (all){
      // Transfer knn distances 
      NumericVector dists = NumericVector(k);
      for (int ki = 0; ki < k; ++ki){ dists[ki] = pr_queue->ith_smallest_key(ki); }
      knn_dist(i, _) = dists;
      
      // And knn ids 
      IntegerVector ids = IntegerVector(k);
      for (int ki = 0; ki < k; ++ki){ ids[ki] = pr_queue->ith_smallest_info(ki); }
      knn_ids(i, _) = ids + 1; 
    } else { 
      knn_dist[i] = pr_queue->max_key(); 
      knn_ids[i] = pr_queue->ith_smallest_info(k-1)+1; // conversion to 1-based
    }
    pr_queue->reset();
  }
  
  List res = List::create(_["dist"] = knn_dist, _["id"] = knn_ids, _["k"] = k, _["sort"] = true);
  StringVector class_names = StringVector(2);
  class_names.at(0) = "kNN";
  class_names.at(1) = "NN";
  res.attr("class") = class_names;
  return(res);
}


/*** R
# x <- cbind(rnorm(40000), rnorm(40000))
# truth <- dbscan::kNN(x, k = 150L)
# dist_x <- dist(x)
# test <- dbscan:::knn_dist(dist_x = dist_x, k = 150L, all = 1L)
# all(truth$dist == test$dist)
# all(truth$id == test$id)
  */
