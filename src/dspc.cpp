#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"
#include <string>
#include <unordered_map>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
StringVector intToStr(IntegerVector iv){
  StringVector res = StringVector(iv.length());
  int ci = 0; 
  for (IntegerVector::iterator i = iv.begin(); i != iv.end(); ++i){
    res[ci++] = patch::to_string(*i);
  }
  return(res);
}

std::unordered_map<std::string, double> toMap(List map){
  std::vector<std::string> keys = map.names();
  std::unordered_map<std::string, double> hash_map = std::unordered_map<std::string, double>();
  const int n = map.size();
  for (int i = 0; i < n; ++i){
    hash_map.emplace((std::string) keys.at(i), (double) map.at(i));
  }
  return(hash_map);
}

NumericVector retrieve(StringVector keys, std::unordered_map<std::string, double> map){
  int n = keys.size(), i = 0;
  NumericVector res = NumericVector(n);
  for (StringVector::iterator it = keys.begin(); it != keys.end(); ++it){ res[i++] = map[as< std::string >(*it)]; }
  return(res);
}


// Density Separation code
// [[Rcpp::export]]
NumericMatrix dspc(List config) {
  
  // Load configuration from list
  const int n = config["n"];
  const int ncl = config["ncl"]; 
  const int n_pairs = config["n_pairs"]; 
  List node_ids = config["node_ids"];
  List acp = config["acp"];
  NumericVector xdist = config["xdist"];
  
  // Conversions and basic setup 
  std::unordered_map<std::string, double> acp_map = toMap(acp);
  double min_mrd = std::numeric_limits<double>::infinity(); 
  NumericMatrix min_mrd_dist = NumericMatrix(n_pairs, 3);
  
  // Loop through cluster combinations, and for each combination 
  int c = 0; 
  for (int ci = 0; ci < ncl; ++ci) {
    for (int cj = (ci+1); cj < ncl; ++cj){
      Rcpp::checkUserInterrupt();
      IntegerVector i_idx = node_ids[ci], j_idx = node_ids[cj]; // i and j cluster point indices
      for (IntegerVector::iterator i = i_idx.begin(); i != i_idx.end(); ++i){
        for (IntegerVector::iterator j = j_idx.begin(); j != j_idx.end(); ++j){
          const int lhs = *i < *j ? *i : *j, rhs = *i < *j ? *j : *i;
          double dist_ij = xdist[INDEX_TF(n, lhs - 1, rhs - 1)]; // dist(p_i, p_j)
          double acd_i = acp_map[patch::to_string(*i)]; // all core distance for p_i
          double acd_j = acp_map[patch::to_string(*j)]; // all core distance for p_i
          double mrd_ij = std::max(std::max(acd_i, acd_j), dist_ij); // mutual reachability distance of the pair
          if (mrd_ij < min_mrd){
            min_mrd = mrd_ij; 
          }
        }
      }
      min_mrd_dist(c++, _) = NumericVector::create(ci+1, cj+1, min_mrd);
      min_mrd = std::numeric_limits<double>::infinity(); 
    }
  }
  return(min_mrd_dist);
}


/*** R
# intToStr(1:4)
# 
# # acp_dist_map[as.character(1:3)]
# # test_func(1:3, acp_dist_map)
# km <- structure(as.list(4:6), names = as.character(1:3))
# test_int <- 1:3
# 
# config <- list(n = 3, ncl = 3L, node_ids = as.list(1:3), acp = km, xdist=dist(1:3))
# dspc(config)
# # toMap(km)
*/
