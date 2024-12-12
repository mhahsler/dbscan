#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Includes
#include "utilities.h"
#include "mst.h"
#include "ANN/ANN.h"
#include "kNN.h"
#include <string>
#include <unordered_map>


// inline indexing of distance matrix
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

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


NumericVector dist_subset_arma(const NumericVector& dist, IntegerVector idx){
  // vec v1 = as<vec>(v1in);
  // uvec idx = as<uvec>(idxin) - 1;
  // vec subset = v1.elem(idx);
  // return(wrap(subset));
  return(NumericVector::create());
}



// Provides a fast of extracting subsets of a dist object. Expects as input the full dist
// object to subset 'dist', and a (1-based!) integer vector 'idx' of the points to keep in the subset
// [[Rcpp::export]]
NumericVector dist_subset(const NumericVector& dist, IntegerVector idx){
  const int n = dist.attr("Size");
  const int cl_n = idx.length();
  NumericVector new_dist = Rcpp::no_init((cl_n * (cl_n - 1))/2);
  int ii = 0;
  for (IntegerVector::iterator i = idx.begin(); i != idx.end(); ++i){
    for (IntegerVector::iterator j = i; j != idx.end(); ++j){
      if (*i == *j) { continue; }
      const int ij_idx = INDEX_TF(n, (*i < *j ? *i : *j) - 1, (*i < *j ? *j : *i) - 1);
      new_dist[ii++] = dist[ij_idx];
    }
  }
  new_dist.attr("Size") = cl_n;
  new_dist.attr("class") = "dist";
  return(new_dist);
}

// Returns true if a given distance is less than 32-bit floating point precision
bool remove_zero(ANNdist cdist){
  return(cdist <= std::numeric_limits<float>::epsilon());
}

ANNdist inv_density(ANNdist cdist){
  return(1.0/cdist);
}

// // [[Rcpp::export]]
// List all_pts_core_sorted_dist(const NumericMatrix& sorted_dist, const List& cl, const int d, const bool squared){
//   // The all core dists to return
//   List all_core_res = List(cl.size());
//
//   // Do the kNN searches per cluster; note that k varies with the cluster
//   int i = 0;
//   for (List::const_iterator it = cl.begin(); it < cl.end(); ++it, ++i){
//     const IntegerVector& cl_pts = (*it);
//     const int k = cl_pts.length();
//
//     // Initial vector to record the per-point all core dists
//     NumericVector all_core_cl = Rcpp::no_init_vector(k);
//
//     // For each point in the cluster, get the all core points dist
//     int j = 0;
//     for (IntegerVector::const_iterator pt_id = cl_pts.begin(); pt_id != cl_pts.end(); ++pt_id, ++j){
//       const NumericMatrix::ConstColumn& knn_dist = sorted_dist.column((*pt_id) - 1);
//
//       // Calculate the all core points distance for this point
//       std::vector<ANNdist> ndists = std::vector<ANNdist>(knn_dist.begin(), knn_dist.begin()+k);
//       std::remove_if(ndists.begin(), ndists.end(), remove_zero);
//       std::transform(ndists.begin(), ndists.end(), ndists.begin(), [=](ANNdist cdist){ return std::pow(1.0/cdist, d); });
//       ANNdist sum_inv_density = std::accumulate(ndists.begin(), ndists.end(), (ANNdist) 0.0);
//       double acdist = std::pow(sum_inv_density/(k - 1.0), -(1.0 / double(d))); // Apply all core points equation
//       all_core_cl[j] = acdist;
//       // return(List::create(_["ndists"] = acdist, _["denom"] = sum_inv_density/(k - 1.0), _["k"] = k));
//     }
//     all_core_res[i] = all_core_cl;
//   }
//   return(all_core_res);
// }

// // [[Rcpp::export]]
// List all_pts_core(const NumericMatrix& data, const List& cl, const bool squared){
//   // copy data
//   int nrow = data.nrow();
//   int ncol = data.ncol();
//   ANNpointArray dataPts = annAllocPts(nrow, ncol);
//   for(int i = 0; i < nrow; i++){
//     for(int j = 0; j < ncol; j++){
//       (dataPts[i])[j] = data(i, j);
//     }
//   }
//
//   // create kd-tree (1) or linear search structure (2)
//   ANNpointSet* kdTree = new ANNkd_tree(dataPts, nrow, ncol, 30, (ANNsplitRule)  5);
//
//   // The all core dists to
//   List all_core_res = List(cl.size());
//
//   // Do the kNN searches per cluster; note that k varies with the cluster
//   int i = 0;
//   for (List::const_iterator it = cl.begin(); it < cl.end(); ++it, ++i){
//     const IntegerVector& cl_pts = (*it);
//     const int k = cl_pts.length();
//
//     // Initial vector to record the per-point all core dists
//     NumericVector all_core_cl = Rcpp::no_init_vector(k);
//
//     // For each point in the cluster, get the all core points dist
//     int j = 0;
//     ANNdistArray dists = new ANNdist[k];
//     ANNidxArray nnIdx = new ANNidx[k];
//     for (IntegerVector::const_iterator pt_id = cl_pts.begin(); pt_id != cl_pts.end(); ++pt_id, ++j){
//       // Do the search
//       ANNpoint queryPt = dataPts[(*pt_id) - 1]; // use original data points
//       kdTree->annkSearch(queryPt, k, nnIdx, dists);
//
//       // V2.
//       std::vector<ANNdist> ndists = std::vector<ANNdist>(dists, dists+k);
//       std::remove_if(ndists.begin(), ndists.end(), remove_zero);
//       std::transform(ndists.begin(), ndists.end(), ndists.begin(), [=](ANNdist cdist){ return std::pow(1.0/cdist, ncol); });
//       ANNdist sum_inv_density = std::accumulate(ndists.begin(), ndists.end(), (ANNdist) 0.0);
//       double acdist = std::pow(sum_inv_density/(k - 1.0), -(1.0 / double(ncol))); // Apply all core points equation
//       all_core_cl[j] = acdist;
//       // return(List::create(_["ndists"] = acdist, _["denom"] = sum_inv_density/(k - 1.0), _["k"] = k));
//     }
//     delete [] dists;
//     delete [] nnIdx;
//     all_core_res[i] = all_core_cl;
//   }
//
//   // cleanup
//   delete kdTree;
//   annDeallocPts(dataPts);
//   annClose();
//
//   // Return the all point core distance
//   if(!squared){ for (int i = 0; i < cl.size(); ++i){ all_core_res[i] = Rcpp::sqrt(all_core_res[i]); } }
//   return(all_core_res);
// }



// NumericVector all_pts_core(const NumericVector& dist, IntegerVector cl, const int d){
//   const int n = dist.attr("Size");
//   const int cl_n = cl.length();
//   NumericVector all_pts_cd = NumericVector(cl_n);
//   NumericVector tmp = NumericVector(cl_n);
//   int knn_i = 0, ii = 0;
//   for (IntegerVector::iterator i = cl.begin(); i != cl.end(); ++i){
//     for (IntegerVector::iterator j = cl.begin(); j != cl.end(); ++j){
//       if (*i == *j) { continue; }
//       const int idx = INDEX_TF(n, (*i < *j ? *i : *j) - 1, (*i < *j ? *j : *i) - 1);
//       double dist_ij = dist[idx];
//       tmp[knn_i++] = 1.0 / (dist_ij == 0.0 ? std::numeric_limits<double>::epsilon() : dist_ij);
//     }
//     all_pts_cd[ii++] = pow(sum(pow(tmp, d))/(cl_n - 1.0), -(1.0 / d));
//     knn_i = 0;
//   }
//   return(all_pts_cd);
// }


// [[Rcpp::export]]
Rcpp::LogicalVector XOR(Rcpp::LogicalVector lhs, Rcpp::LogicalVector rhs) {
  R_xlen_t i = 0, n = lhs.size();
  Rcpp::LogicalVector result(n);
  for ( ; i < n; i++) {  result[i] = (lhs[i] ^ rhs[i]); }
  return result;
}

// [[Rcpp::export]]
NumericMatrix dspc(const List& cl_idx, const List& internal_nodes, const IntegerVector& all_cl_ids, const NumericVector& mrd_dist) {

  // Setup variables
  const int ncl = cl_idx.length(); // number of clusters
  NumericMatrix res = Rcpp::no_init_matrix((ncl * (ncl - 1))/2, 3); // resulting separation measures

  // Loop through cluster combinations, and for each combination
  int c = 0;
  double min_edge = std::numeric_limits<double>::infinity();
  for (int ci = 0; ci < ncl; ++ci) {
    for (int cj = (ci+1); cj < ncl; ++cj){
      Rcpp::checkUserInterrupt();

      // Do lots of indexing to get the relative indexes corresponding to internal nodes
      const IntegerVector i_idx = internal_nodes[ci], j_idx = internal_nodes[cj]; // i and j cluster point indices
      const IntegerVector rel_i_idx = match(as<IntegerVector>(cl_idx[ci]), all_cl_ids)[i_idx - 1];
      const IntegerVector rel_j_idx = match(as<IntegerVector>(cl_idx[cj]), all_cl_ids)[j_idx - 1];
      IntegerVector int_idx = combine(rel_i_idx, rel_j_idx);

      // Get the pairwise MST
      NumericMatrix pairwise_mst = mst_prims(dist_subset(mrd_dist, int_idx), int_idx.length());

      // Do lots of indexing / casting
      const IntegerVector from_int = seq_len(rel_i_idx.length());
      const NumericVector from_idx = as<NumericVector>(from_int);
      const NumericVector from = pairwise_mst.column(0), to = pairwise_mst.column(1), height = pairwise_mst.column(2);

      // Find which distances in the MST cross to both clusters
      LogicalVector cross_edges = XOR(Rcpp::in(from, from_idx), Rcpp::in(to, from_idx));

      // The minimum weighted edge of these cross edges is the density separation between the two clusters
      min_edge = min(as<NumericVector>(height[cross_edges]));

      // Save the minimum edge
      res(c++, _) = NumericVector::create(ci+1, cj+1, min_edge);
      min_edge = std::numeric_limits<double>::infinity();
    }
  }
  return(res);
}


// Density Separation code
// NumericMatrix dspc(List config, const NumericVector& xdist) {
//
//   // Load configuration from list
//   const int n = config["n"];
//   const int ncl = config["ncl"];
//   const int n_pairs = config["n_pairs"];
//   List node_ids = config["node_ids"];
//   List acp = config["acp"];
//
//   // Conversions and basic setup
//   std::unordered_map<std::string, double> acp_map = toMap(acp);
//   double min_mrd = std::numeric_limits<double>::infinity();
//   NumericMatrix min_mrd_dist = NumericMatrix(n_pairs, 3);
//
//   // Loop through cluster combinations, and for each combination
//   int c = 0;
//   for (int ci = 0; ci < ncl; ++ci) {
//     for (int cj = (ci+1); cj < ncl; ++cj){
//       Rcpp::checkUserInterrupt();
//       IntegerVector i_idx = node_ids[ci], j_idx = node_ids[cj]; // i and j cluster point indices
//       for (IntegerVector::iterator i = i_idx.begin(); i != i_idx.end(); ++i){
//         for (IntegerVector::iterator j = j_idx.begin(); j != j_idx.end(); ++j){
//           const int lhs = *i < *j ? *i : *j, rhs = *i < *j ? *j : *i;
//           double dist_ij = xdist[INDEX_TF(n, lhs - 1, rhs - 1)]; // dist(p_i, p_j)
//           double acd_i = acp_map[patch::to_string(*i)]; // all core distance for p_i
//           double acd_j = acp_map[patch::to_string(*j)]; // all core distance for p_i
//           double mrd_ij = std::max(std::max(acd_i, acd_j), dist_ij); // mutual reachability distance of the pair
//           if (mrd_ij < min_mrd){
//             min_mrd = mrd_ij;
//           }
//         }
//       }
//       min_mrd_dist(c++, _) = NumericVector::create(ci+1, cj+1, min_mrd);
//       min_mrd = std::numeric_limits<double>::infinity();
//     }
//   }
//   return(min_mrd_dist);
// }


