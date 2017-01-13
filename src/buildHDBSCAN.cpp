#include <Rcpp.h>
using namespace Rcpp;

#include <stack>

// std::to_string is apparently a c++11 only thing that crashes appveyor, so using ostringstream it is!
namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


template <typename T> bool contains (const T& container, std::string key)
{
  if (std::find(container.begin(), container.end(), key) != container.end()){
    return true; 
  } else {
    return false; 
  }
}

IntegerVector which_cpp( NumericVector x, double value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

IntegerVector which_cpp( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] == value) y.push_back(i); }
  return wrap(y);
}

IntegerVector which_geq( IntegerVector x, int value) {
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  for(int i = 0; i < nx; ++i) { if (x[i] >= value) y.push_back(i); }
  return wrap(y);
}


// [[Rcpp::export]]
NumericVector combine(const NumericVector& t1, const NumericVector& t2){
  std::size_t n = t1.size() + t2.size();
  NumericVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}

IntegerVector combine(const IntegerVector& t1, const IntegerVector& t2){
  std::size_t n = t1.size() + t2.size();
  IntegerVector output = Rcpp::no_init(n);
  std::copy(t1.begin(), t1.end(), output.begin());
  std::copy(t2.begin(), t2.end(), output.begin()+t1.size());
  return output;
}




// [[Rcpp::export]]
List buildDendrogram(List hcl, int minPts = 5) {
  
  // Extract relevent info
  NumericMatrix merge = hcl["merge"]; 
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"]; 
  int n = merge.nrow() + 1; 
  
  List new_br; 
  List z = List::create();
  StringVector labels = as<StringVector>(pt_order);
  
  List bcontains = List(), bheight = List();
  for (int k = 0; k < n-1; ++k){
    //String pt_label = labels.at(k);
    bcontains[patch::to_string(k)] = IntegerVector();
    bheight[patch::to_string(k)] = IntegerVector(); 
  }
  List hdbscan = List(); 
  
  
  // Global variables to keep track of left+right merge points and distance values 
  IntegerVector l_contains, r_contains; 
  NumericVector l_eps, r_eps; 
  
  
  int k, cid = 0, ccid=0, bid = 0; 
  for (k = 0; k < n-1; k++){ 
    
    std::cout << "bid:" << bid << std::endl; 
    int lm = merge(k, 0), rm = merge(k, 1);
    Rcout << "lm: " << lm << ", rm: " << rm << std::endl; 
    std::string l_str = patch::to_string(lm), r_str = patch::to_string(rm);
    IntegerVector m(2);
    m.at(0) = merge(k, 0), m.at(1) = merge(k, 1);
    
    // First Case: Both are singletons, so need to create leaves
    if (all(m < 0).is_true()){
      // Left 
      IntegerVector left = IntegerVector::create(-lm);
      left.attr("label") = patch::to_string(-lm); //labels.at((-lm) - 1);
      left.attr("members") = 1;
      left.attr("height") = 0;
      left.attr("leaf") = true;

      // Right 
      IntegerVector right = IntegerVector::create(-rm);
      right.attr("label") = patch::to_string(-rm);
      right.attr("members") = 1;
      right.attr("height") = 0;      
      right.attr("leaf") = true;

      // Merge 
      new_br = List::create(left, right);
      new_br.attr("midpoint") = 0.5;
      new_br.attr("members") = 2;
      new_br.attr("class") = "dendrogram";
      new_br.attr(".bid") = ++bid;
      
      // Store info for later to prevent recursion
      l_contains = left, r_contains = right; 

    } 
    // Second case: 1 is a singleton, the other is a branch
    else if (any(m < 0).is_true()){
      bool isL = lm < 0;

      // Create the leaf from the negative entry
      List leaf = List::create(isL ? -lm : -rm);
      leaf.attr("members") = 1;
      leaf.attr("leaf") = true;
      leaf.attr("height") = 0;
      leaf.attr("label") = patch::to_string(isL ? -lm : -rm);

      // Merge the leaf with the other existing branch
      std::string branch_key = isL ? patch::to_string(rm - 1) : patch::to_string(lm - 1);
      List sub_branch = z[branch_key];
      new_br = isL ? List::create(leaf, sub_branch) : List::create(sub_branch, leaf);
      z[branch_key] = R_NilValue;

      // Set attributes of new branch
      int sub_members = sub_branch.attr("members");
      double mid_pt = sub_branch.attr("midpoint");
      new_br.attr("members") = int(sub_members) + 1;
      new_br.attr("midpoint") = (int(isL ? 0 : sub_members) + mid_pt) / 2;
      new_br.attr("class") = "dendrogram";
      new_br.attr(".bid") = ++bid;

      int cbid = sub_branch.attr(".bid");
      IntegerVector cb_contains = bcontains[patch::to_string(cbid)];
      NumericVector cb_height = bheight[patch::to_string(cbid)];
      

      // If the sub branch was part of a cluster, record the cluster info
      if (!Rf_isNull(sub_branch.attr(".last_cid"))){
        ccid = sub_branch.attr(".last_cid");
        List check = hdbscan[patch::to_string(ccid)];

        NumericVector eps = check["eps"];
        IntegerVector contains = check["contains"];


        // Push the noise info
        eps.push_back(eps_dist.at(k));
        contains.push_back(isL ? -lm : -rm);
      }
      
      l_contains = IntegerVector::create(isL ? -lm : -rm), r_contains = cb_contains;
      l_eps = NumericVector::create(eps_dist.at(k)), r_eps = cb_height;
    } else {
      // Create the new branch
      List l_branch = z[patch::to_string(lm - 1)];
      List r_branch = z[patch::to_string(rm - 1)];
      new_br = List::create(l_branch, r_branch);
      //new_br.attr(".bid") = ++bid;

      // Store attribute valeus in locals to get around proxy
      int left_members = l_branch.attr("members");
      int right_members = r_branch.attr("members");
      double l_mid = l_branch.attr("midpoint");
      double r_mid = r_branch.attr("midpoint");
      int l_bid = l_branch.attr(".bid");
      int r_bid = r_branch.attr(".bid");

      l_contains = bcontains[patch::to_string(l_bid)], r_contains = bcontains[patch::to_string(r_bid)]; 
      l_eps = bheight[patch::to_string(l_bid)], r_eps = bheight[patch::to_string(r_bid)]; 
      
      // Agglomerative HDBSCAN rule: if both members have >= minPts, then it was a true split,
      // so the next (going up) cluster needs to be formed and initialized with the current eps
      if (left_members >= minPts && right_members >= minPts){
        
        // Mark the left branch as a cluster 
        List l_cluster = List::create(_["contains"] = l_contains, _["eps"] = l_eps);
        hdbscan[patch::to_string(++cid)] = l_cluster;
        l_branch.attr(".last_cid") = cid;
        
        // Mark the right branch as a cluster 
        List r_cluster = List::create(_["contains"] = r_contains, _["eps"] = r_eps);
        hdbscan[patch::to_string(++cid)] = r_cluster;
        l_branch.attr(".last_cid") = cid;
        
        // Create new cluster going up 
        List check = List::create(
          _["contains"] = combine(l_contains, r_contains),
          _["eps"] = NumericVector(left_members + right_members, eps_dist.at(k))
        );
        hdbscan[patch::to_string(++cid)] = check; // note the CID increment
        new_br.attr(".last_cid") = cid;
      } else {
        Rcout << "here4" << std::endl; 
        // One branch is noise, the other is an existing cluster
        bool isL = left_members < minPts;
        Rcout << "here5: " << isL << std::endl; 
        
        
        if ((isL ? r_branch.attr(".last_cid") : l_branch.attr(".last_cid")) == R_NilValue){
          
        } else {
          ccid = isL ? r_branch.attr(".last_cid") : l_branch.attr(".last_cid");
          Rcout << "ccid: " << ccid << std::endl; 
          List check = hdbscan[patch::to_string(ccid)];
          NumericVector eps = check["eps"];
          IntegerVector contains = check["contains"];
          
          //  get branch info
          List sub_branch = isL ? r_branch : l_branch;
          int cbid = sub_branch.attr(".bid");
          IntegerVector cb_contains = bcontains[patch::to_string(cbid)];
          NumericVector cb_height = bheight[patch::to_string(cbid)];
          check["contains"] = combine(contains, cb_contains);
          check["eps"] = combine(eps, cb_height);
          
          Rcout << "here5" << std::endl; 
          
          // Deallocate old branch info
          bcontains[patch::to_string(cbid)] = R_NilValue;
          bheight[patch::to_string(cbid)] = R_NilValue;
          
          // Retain original cluster ID
          Rcout << "here6" << std::endl; 
          new_br.attr(".last_cid") = sub_branch.attr(".last_cid");
        }
        
      }

      // Set up new branch attributes
      new_br.attr("members") = left_members + right_members;
      new_br.attr("midpoint") = (left_members + l_mid + r_mid) / 2;
      new_br.attr("class") = "dendrogram";
      new_br.attr(".bid") = ++bid;
      
      // record the branch info 

      // *Experimental* deallocate unneeded memory
      z[patch::to_string(lm - 1)] = R_NilValue;
      z[patch::to_string(rm - 1)] = R_NilValue;
    }
    //new_br.attr(".bid") = ++bid; 
    
    new_br.attr("height") = eps_dist.at(k);
    z[patch::to_string(k)] = new_br;
  }
  Rcout << k << std::endl; 
  List res = z[patch::to_string(k - 1)];
  res.attr("class") = "dendrogram";
  return(List::create(_["dend"] = res, _["hdbscan"]=hdbscan, 
                      _["bcontains"] = bcontains, _["bheight"] = bheight));
  //return(labels);
  // "The minimum epsilon that x_i was still considered as part of the cluster 
}

// [[Rcpp::export]]
List n_members(List hcl, int minPts = 5){
  // Extract relevent info
  NumericMatrix merge = hcl["merge"]; 
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"]; 
  int n = merge.nrow() + 1; 
  
  // First pass - agglomerative step: record information like member sizes, point indices for 
  // agglomeration of non-singletons, etc.
  IntegerVector member_sizes = IntegerVector(n-1, 0); 
  List contains = List();
  List contains_k = List();
  List hdbscan = List(); 
  List eps_index = List(); 
  List eps_death = List(); 
  List eps_info; 
  int cid = 0; 
  
  IntegerVector cl_tracker = IntegerVector(n-1 , 0);
  List cl_hierarchy = List(), cl_merge = List();
  
  IntegerVector c_contains, ck_contains; 
  
  int cl_cid = 0; 
  for (int k = 0; k < n-1; ++k){
    
    // Current Merge
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm);
    
    //Rcout << "lm: " << lm << ", rm: " << rm << std::endl; 
    //std::string l_str = patch::to_string(lm), r_str = patch::to_string(rm);
  
    // Trivial case: merge of singletons, create a temporary *assumed* cluster to be resolved on merger
    if (all(m < 0).is_true()){
      member_sizes[k] = 2;  
      cl_tracker.at(k) = ++cl_cid;
      contains[patch::to_string(cl_cid)] = IntegerVector::create(-lm, -rm);
      contains_k[patch::to_string(cl_cid)] = IntegerVector::create(k, k);
      
      // If minPts == 2, then the epsilon death distance has been reached
      if (minPts == 2){
        eps_death[patch::to_string(cl_cid)] = eps_dist.at(k);
      } 
      // Rcout << "Creating: " << cl_cid << std::endl;  
    } else if (any(m < 0).is_true()) {
      int pos_merge = (lm < 0 ? rm : lm), merge_size = member_sizes[pos_merge - 1];
      member_sizes[k] = merge_size + 1; 
      
      std::string contains_index = patch::to_string(cl_tracker.at(pos_merge - 1)); 
      c_contains = contains[contains_index]; 
      c_contains.push_front(-(lm < 0 ? lm : rm));
      contains[contains_index] = c_contains;//combine(IntegerVector::create(), );
      
      ck_contains = contains_k[contains_index];
      ck_contains.push_front(k);
      contains_k[contains_index] = ck_contains;
      
      cl_tracker.at(k) = cl_tracker.at((lm < 0 ? rm : lm) - 1);

      // Rcout << "Trying: " << contains_index << std::endl; 
      // eps_info = eps_index[contains_index];
      // if (eps_info["eps_death"] == R_NilValue && merge_size + 1 >= minPts){
      //   eps_info["death"] = eps_dist.at(k);
      // }
      if (!eps_death.containsElementNamed(contains_index.c_str()) && merge_size + 1 >= minPts){
        eps_death[contains_index] = eps_dist.at(k); 
      }
    } else {
      // Record Member Sizes 
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];
      member_sizes[k] = merge_size1 + merge_size2; 
      
      // Cluster ID's of branches
      std::string l_index = patch::to_string(cl_tracker.at(lm - 1)); 
      std::string r_index = patch::to_string(cl_tracker.at(rm - 1)); 
      
      // Get contains and indices 
      IntegerVector l_contains = contains[l_index], r_contains = contains[r_index]; 
      IntegerVector l_ci = contains_k[l_index], r_ci = contains_k[r_index]; 
      
      // The HDBSCAN step 
      if (merge_size1 >= minPts && merge_size2 >= minPts){

        // Merge was a true split point, so increment the cluster ID and record 
        cl_tracker.at(k) = ++cl_cid;
        std::string cl_cid_str = patch::to_string(cl_cid); 
        cl_hierarchy[cl_cid_str] = List::create(
          _["contains"] = IntegerVector::create(cl_tracker.at(lm - 1), cl_tracker.at(rm - 1))
        );
        
        // Record integer IDs of containing points & merge steps
        contains[cl_cid_str] = combine(l_contains, r_contains);
        contains_k[cl_cid_str] = combine(l_ci, r_ci);
        eps_death[cl_cid_str] = eps_dist.at(k);
        
        // Record eps births 
        // eps_info = eps_index[l_index]; 
        // eps_info["eps_birth"] = eps_dist.at(lm - 1);
        // eps_info = eps_index[r_index]; 
        // eps_info["eps_birth"] = eps_dist.at(rm - 1);
      //   // Left Cluster 
      //   NumericVector l_eps = NumericVector(l_contains.size());
      //   for (int i = 0; i < l_contains.size(); ++i){ l_eps.at(i) = eps_dist.at(l_ci.at(i)); }
      //   hdbscan[patch::to_string(cl_tracker.at(lm - 1))] = List::create(
      //     _["eps"] = l_eps, 
      //     _["contains"] = l_contains
      //   );

        
        // Right Cluster 
        // NumericVector r_eps = NumericVector(r_contains.size());
        // for (int i = 0, r_size = r_contains.size(); i < r_size; ++i){ r_eps.at(i) = eps_dist.at(r_ci.at(i)); }
        // hdbscan[patch::to_string(cl_tracker.at(rm - 1))] = List::create(
        //   _["eps"] = r_eps,
        //   _["contains"] = r_contains
        // );
      } else {
        // Small correction: Two cluster IDs assigned to noise points that were merged to form a single 
        // cluster. Correct by choosing one of the IDs and replacing the other.
        for (int i = 0, right = cl_tracker.at(rm - 1), n=cl_tracker.length(); i < n; ++i){
          if (cl_tracker[i] == right) { cl_tracker[i] = cl_tracker.at(lm - 1); }
        }
        cl_tracker.at(k) = cl_tracker.at(lm - 1);
        
        // Keep left index
        contains[l_index] = combine(l_contains, r_contains);
        contains_k[l_index] = combine(l_ci, r_ci);
      
        
        // eps_info = eps_index[l_index];
        // if (eps_info["eps_death"] == R_NilValue && l_contains.length() + r_contains.length() >= minPts){
        //   eps_info["eps_death"] = eps_dist.at(k);
        // }
        
        // Remove right index to cut down unnecessary memory usage
        contains[r_index] = R_NilValue;
        contains_k[r_index] = R_NilValue;
        // eps_index[r_index] = R_NilValue;
        
        // Adjust cl_cid 
      }
    }
    // 
    // List eps_info = eps_index[patch::to_string(cl_cid)]; 
    // if (eps_info["eps_death"] == R_NilValue){
    //   eps_info["eps_death"] = eps_dist.at(k);
    //   
    // }
  }
  
  // Post processing 1: Gather up relevent information, calculate stability scores, etc. 
  NumericVector cl_eps; 
  IntegerVector cl_contains, cl_k;
  StringVector cl_ids = as<StringVector>(contains.names());
  for (int i = 0; i < cl_ids.size(); ++i){
    std::string cluster_id = Rcpp::as< std::string >(cl_ids.at(i));
    if (contains[cluster_id] != R_NilValue){
      // Get cluster point indices and cluster labels
      cl_contains = contains[cluster_id], cl_k = contains_k[cluster_id];
      IntegerVector cl_merge_steps = which_cpp(cl_tracker, std::stoi(cluster_id));
      
      // Retrieve the distance the cluster formed (going down the hierarchy) 
      double eps_birth = eps_dist.at(max(cl_merge_steps));
      //IntegerVector cluster_size = member_sizes[cl_merge_steps];
      
      // Fill in eps values 
      cl_eps = NumericVector(cl_contains.size());
      for (int i = 0, cl_size = cl_contains.size(); i < cl_size; ++i){
        cl_eps.at(i) = std::max((double) eps_death[cluster_id], (double) eps_dist.at(int(cl_k[i]))); // never record distances below eps_death
      }
      
      
      // Save for return 
      NumericVector scores = (1/cl_eps) - (1/eps_birth); 
      hdbscan[cluster_id] = List::create(
        _["eps"] = cl_eps,
        _["eps_birth"] = eps_birth, 
        _["eps_death"] = eps_death[cluster_id],
        _["contains"] = cl_contains, 
        _["score"] = (double) sum(scores)
      );
    }
  }
  
  return(List::create(
      _["member_sizes"] = member_sizes, 
      _["contains"] = contains, 
      _["contains_k"] = contains_k,
      _["clusters"] = hdbscan, 
      _["cl_tracker"] = cl_tracker, 
      _["cl_hierarchy"] = cl_hierarchy, 
      _["merge_steps"] = contains_k, 
      _["names"] = cl_ids, 
      _["eps_info"] = eps_index, 
      _["eps_death"] = eps_death));
}

// [[Rcpp::export]]
IntegerVector all_children(List hier, int key, IntegerVector sub_tree){
  IntegerVector res = IntegerVector(); 
  
  List elem = hier[patch::to_string(key).c_str()];
  IntegerVector children = elem["contains"];
  int n_children; 
  bool finished;
  
  while (!finished){
    finished = true; 
    for (n_children = 0; n_children < children.length(); ++n_children){
      int child_id = children.at(n_children);
      res.push_back(child_id);
      if (hier.containsElementNamed(patch::to_string(child_id).c_str())){
        finished = false; 
      }
    }
  }
 
  
  // if (!hier.containsElementNamed(patch::to_string(key).c_str())){
  //   sub_tree.push_back(key);
  // } else {
  //   
  //  
  //   for (int i = 0; i < children.size(); ++i){
  //     all_children(hier, children.at(i), sub_tree);
  //   }
  // }
  // return(sub_tree);
}//dbscan:::all_children()

// Post processing 2: Compute stability scores for cluster objects in the hierarchy
double computeSalientScores(List info, std::string cid, IntegerVector res){
  // List clusters = info["clusters"];
  // List cl_hierarchy = info["cl_hierarchy"];
  // IntegerVector cl_tracker = info["cl_tracker"];
  // 
  // // Base case 
  // if (!cl_hierarchy.containsElementNamed(cid.c_str())){
  //   List cl = clusters[cid];
  //   res.push_back(cid); // assume the leaf will be a salient cluster until proven otherwise 
  //   return((double) cl["score"]);
  // } else {
  // // Non-base case: at a merge of clusters, determine which to keep  
  //   List cl = clusters[cid];
  //   NumericVector child_scores = NumericVector();
  //   IntegerVector child_ids = cl_hierarchy[patch::to_string(cid)];
  //   for (int i = 0; i < child_ids.length(); ++i){
  //     int child_id = child_ids.at(i);
  //     child_scores.push_back(computeSalientScores(info, patch::to_string(child_id), res));
  //   }
  //   double old_score = (double) cl["score"]; 
  //   double new_score = std::max(old_score, (double) sum(child_scores));
  //   cl["new_score"] = new_score;
  //   
  //   // If the score is unchanged, this cluster is better, so remove the children
  //   if (new_score == old_score) {
  //     for (int i = 0; i < child_ids.length(); ++i){
  //       res = res[res != child_ids.at(i)];
  //     }
  //     //res.push_back(cid);
  //   }
  //   return(new_score);
  // }
  // 
  // 
  // std::string c_node = patch::to_string(cl_tracker[cl_tracker.length() - 1]);
  // std::stack<int> depth = std::stack<int>(); 
  // StringVector labels = cl_hierarchy.attr("names");
  // bool is_done = false;
  // while(!is_done){
  //   if (clusters.containsElementNamed(c_node.c_str())){
  //     List cluster_info = clusters[c_node];
  //     for (int i = 0; i < cluster_info.length(); ++i){
  //       int child_id = cluster_info.at(i);
  //       depth.push(child_id);
  //     }
  //   } else {
  //     //depth.push(child_id);
  //   }
  // }
} 


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
dbscan:::buildDendrogram(cl)
dbscan:::n_members(hcl)

post_process <- function(res){
  
  pos_merges <- hcl$merge[which(hcl$merge[, 1] > 0 & hcl$merge[, 2] > 0),]
  true_splits <- which(apply(pos_merges, 1, function(row) res$member_sizes[row[[1]]] > 5 && res$member_sizes[row[[2]]] > 5))
  sapply(unique(res$cl_tracker), function(id) res$contains[max(which(res$cl_tracker == id))])
  cluster_hierarchy <- pos_merges[true_splits,]
  cluster_hierarchy <- cluster_hierarchy[nrow(cluster_hierarchy):1,]
  
  max(res$hdbscan$`3`$eps)
  res$member_sizes[which(res$cl_tracker == 1)] >= 5
  
}

*/
