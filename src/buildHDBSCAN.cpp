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
List buildDendrogram(List hcl) {
  
  // Extract hclust info
  NumericMatrix merge = hcl["merge"]; 
  NumericVector height = hcl["height"];
  IntegerVector order = hcl["order"]; 
  int n = merge.nrow() + 1, k; 
  
  List new_br, z = List::create();
  for (k = 0; k < n-1; k++){ 
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm); 
    
    // First Case: Both are singletons, so need to create leaves
    if (all(m < 0).is_true()){
      // Left 
      IntegerVector left = IntegerVector::create(-lm);
      left.attr("label") = patch::to_string(-lm);
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
      new_br.attr("midpoint") = (int(isL ? 1 : sub_members) + mid_pt) / 2;
      new_br.attr("class") = "dendrogram";
    } else {
      // Create the new branch
      List l_branch = z[patch::to_string(lm - 1)], r_branch = z[patch::to_string(rm - 1)];
      new_br = List::create(l_branch, r_branch);

      // Store attribute valeus in locals to get around proxy
      int left_members = l_branch.attr("members"), right_members = r_branch.attr("members");
      double l_mid = l_branch.attr("midpoint"), r_mid = r_branch.attr("midpoint");

      // Set up new branch attributes
      new_br.attr("members") = left_members + right_members;
      new_br.attr("midpoint") = (left_members + l_mid + r_mid) / 2;
      new_br.attr("class") = "dendrogram";

      // Deallocate unneeded memory along the way
      z[patch::to_string(lm - 1)] = R_NilValue;
      z[patch::to_string(rm - 1)] = R_NilValue;
    }
    new_br.attr("height") = height.at(k);
    z[patch::to_string(k)] = new_br;
  }
  List res = z[patch::to_string(k - 1)];
  res.attr("class") = "dendrogram";
  return(res);
}

// List buildCondensedTree(List info){
//   List clusters = info["clusters"];
//   List cl_hierarchy = info["cl_hierarchy"];
//   IntegerVector cl_tracker = info["cl_tracker"];
//   int root_node = info[".root_node"];
//   
//   IntegerVector children = cl_hierarchy[patch::to_string(root_node).c_str()];
//   List branches = List(clusters.length()); 
//   std::stack<int> to_do = std::stack<int>(); 
//   to_do.push(root_node);
//   
//   List members = List(), midpoints = List(); 
//   std::vector<std::string> cl_labels = (clusters.names());
//   for (int i = 0, cid = 0; i < cl_labels.size(); ++i){
//     if (!cl_hierarchy.containsElementNamed(cl_labels.at(i).c_str())){
//       IntegerVector leaf = IntegerVector::create(++cid); //clusters[patch::to_string(parent)]
//       leaf.attr("label") = cl_labels.at(i);
//       leaf.attr("members") = 1;
//       leaf.attr("height") = 0;
//       leaf.attr("midpoint") = 0.5; 
//       leaf.attr("leaf") = true;
//       branches[cl_labels.at(i)] = leaf;
//       members[cl_labels.at(i)] = 1; 
//       midpoints[cl_labels.at(i)] = 0.5;
//     }
//   }
//   
//   while (to_do.size() != 0){
//     int parent = to_do.top();
//     if (!cl_hierarchy.containsElementNamed(patch::to_string(parent).c_str())){
//       to_do.pop();
//     } else {
//       children = cl_hierarchy[patch::to_string(parent).c_str()];
//       List child_branches = List(children.length());
//       bool process = true, accumulate = true;  
//       for (int n_children = 0, n_members = 0; n_children < children.length(); ++n_children){
//         int child_id = children.at(n_children);
//         if (branches.containsElementNamed(patch::to_string(child_id).c_str())){
//           child_branches.at(n_children) = branches[patch::to_string(child_id)];;
//         } else { 
//           process = false; 
//           to_do.push(child_id);
//         }
//       }
//       if (process){
//         to_do.pop();
//         List parent_cl = clusters[patch::to_string(parent)]; 
//         std::string l_str = patch::to_string(children.at(0)); 
//         std::string r_str = patch::to_string(children.at(1)); 
//         List l_branch = branches[l_str], r_branch = branches[r_str];
//         
//         // Unfortunately this is necessary
//         int l_members = members[l_str]; 
//         int r_members = members[r_str]; 
//         float l_mid = midpoints[l_str]; 
//         float r_mid = midpoints[r_str]; 
//         
//         // Make the new branch
//         child_branches.attr("label") = patch::to_string(parent);
//         child_branches.attr("members") = l_members + r_members;
//         
//         //      attr(zk, "midpoint") <- (.memberDend(zk[[1L]]) + 
//         // attr(z[[X[1 + isL]]], "midpoint"))/2
//         bool isL = (bool) !cl_hierarchy.containsElementNamed(l_str.c_str()); // is left a leaf 
//         if (!isL && cl_hierarchy.containsElementNamed(r_str.c_str())){ // is non-singleton merge
//           child_branches.attr("midpoint") = (l_members + l_mid + r_mid) / 2;
//         } else { // contains a leaf 
//           Rcout << "HERE";
//           int sub_members = isL ? r_members : l_members; 
//           float mid_pt = isL ? r_mid : l_mid; 
//           child_branches.attr("midpoint") = ((isL ? 1 : sub_members) + mid_pt) / 2;
//         }
//         
// 
//         child_branches.attr("height") = (float) parent_cl["eps_birth"];
//         child_branches.attr("class") = "dendrogram";
//         
//         // Need to store this for later 
//         branches[patch::to_string(parent)] = child_branches; 
//         midpoints[patch::to_string(parent)] = (float) child_branches.attr("midpoint"); 
//         members[patch::to_string(parent)] = (float) child_branches.attr("members"); 
//         
//         // Don't need these anymore
//         branches[l_str] = R_NilValue;
//         branches[r_str] = R_NilValue;
//       }
//     }
//   }
//   return(branches[patch::to_string(root_node)]);
// }

// [[Rcpp::export]]
IntegerVector all_children(List hier, int key){
  IntegerVector res = IntegerVector(); 
  
  // If the key doesn't exist return an empty vector
  if (!hier.containsElementNamed(patch::to_string(key).c_str())){
    return(res);
  }
  
  // Else, do iterative 'recursive' type function to extract all the IDs of 
  // all sub trees
  IntegerVector children = hier[patch::to_string(key).c_str()];
  std::stack<int> to_do = std::stack<int>(); 
  to_do.push(key);
  while (to_do.size() != 0){
    int parent = to_do.top();
    if (!hier.containsElementNamed(patch::to_string(parent).c_str())){
      to_do.pop();
    } else {
      children = hier[patch::to_string(parent).c_str()];
      to_do.pop();
      for (int n_children = 0; n_children < children.length(); ++n_children){
        int child_id = children.at(n_children);
        res.push_back(child_id);
        to_do.push(child_id);
      }
    }
  }
  return(res);
}

// Post processing 2: Compute stability scores for cluster objects in the hierarchy
// [[Rcpp::export]]
double computeSalientScores(const List info, std::string cid, std::list<int>& sc){
  List clusters = info["clusters"];
  List cl_hierarchy = info["cl_hierarchy"];
  IntegerVector cl_tracker = info["cl_tracker"];
  
  // Base case
  if (!cl_hierarchy.containsElementNamed(cid.c_str())){
    List cl = clusters[cid];
    sc.push_back(stoi(cid)); // assume the leaf will be a salient cluster until proven otherwise
    return((double) cl["score"]); // Use precomputed scores for leaves 
  } else {
  // Non-base case: at a merge of clusters, determine which to keep
    List cl = clusters[cid];
    
    // Get child stability scores
    NumericVector child_scores = NumericVector();
    IntegerVector child_ids = cl_hierarchy[cid];
    for (int i = 0, clen = child_ids.length(); i < clen; ++i){
      int child_id = child_ids.at(i);
      child_scores.push_back(computeSalientScores(info, patch::to_string(child_id), sc));
    }
    
    // Compare and update stability scores 
    double old_score = (double) cl["score"];
    double new_score = std::max(old_score, (double) sum(child_scores));
    cl["new_score"] = new_score;
    
    // If the score is unchanged, remove the children and add parent
    if (new_score == old_score && stoi(cid) != ((int) info[".root_node"])) {
      IntegerVector children = all_children(cl_hierarchy, stoi(cid));
      for (int i = 0, clen = children.length(); i < clen; ++i){
        sc.remove(children.at(i)); // use list for slightly better random deletion performance
      }
      sc.push_back(stoi(cid));
    }
    
    // Return this sub trees score
    return(new_score);
  }
} 


/* HDBSCAN*
  Processing step to compute all the relevent information needed for HDBSCAN.
  The cluster stability scores are computed via the tree traversal rely on a separate function
  Requires information associated with hclust elements. See ?hclust in R for more info. 
  merge := an (n-1) x 2 matrix representing the MST from the mutual reachability graph
  height := the (epsilon) distance each new set of clusters formed from the MST 
  order := the point indices of the original data the negative entries in merge refer to */ 
// [[Rcpp::export]]
List hdbscan_fast(List hcl, int minPts = 5){
  // Extract hclust info
  NumericMatrix merge = hcl["merge"]; 
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"]; 
  int n = merge.nrow() + 1, cl_cid = 0, k; 
  
  IntegerVector cl_tracker = IntegerVector(n-1 , 0), //  Which cluster does each merge step represent
                member_sizes = IntegerVector(n-1, 0); // Size each step -- may not be needed 
  List contains = List(), // Point indices per cluster 
       contains_k = List(), // Merge step (k) each point participated in, per cluster
       clusters = List(), // Final cluster information
       eps_death = List(), // epsilon deaths, per cluster
       eps_birth = List(), // epsilon births, per cluster
       cl_hierarchy = List(); // Keeps track of hierarchy, which cluster contains who 
  
  // Temporaries 
  IntegerVector c_contains, ck_contains; 
  for (k = 0; k < n-1; ++k){
    // Current Merge
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm);

    // Trivial case: merge of singletons, create a temporary *assumed* cluster to be resolved on merger
    if (all(m < 0).is_true()){
      member_sizes[k] = 2;  
      cl_tracker.at(k) = ++cl_cid;
      contains[patch::to_string(cl_cid)] = IntegerVector::create(-lm, -rm);
      contains_k[patch::to_string(cl_cid)] = IntegerVector::create(k, k);
      
      // If minPts == 2, then the epsilon death distance has been reached
      eps_death[patch::to_string(cl_cid)] = eps_dist.at(k); 
      // if (minPts == 2){
      //   eps_death[] = eps_dist.at(k);
      // } 
    } else if (any(m < 0).is_true()) {
      int pos_merge = (lm < 0 ? rm : lm), merge_size = member_sizes[pos_merge - 1];
      
      // Record member sizes
      member_sizes[k] = merge_size + 1; 
      
      // Which cluster is being agglomerate
      std::string contains_index = patch::to_string(cl_tracker.at(pos_merge - 1)); 
      c_contains = contains[contains_index]; 
      c_contains.push_front(-(lm < 0 ? lm : rm));
      contains[contains_index] = c_contains;
      
      // Record the merge step 
      ck_contains = contains_k[contains_index];
      ck_contains.push_front(k);
      contains_k[contains_index] = ck_contains; // necessary: do not remove 
      cl_tracker.at(k) = cl_tracker.at((lm < 0 ? rm : lm) - 1);
      
      // If the merge causes more than minPts to be fused, record density-level as eps_death 
      // if (!eps_death.containsElementNamed(contains_index.c_str()) && merge_size + 1 >= minPts){
      //   eps_death[contains_index] = eps_dist.at(k); 
      // }
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
        cl_hierarchy[cl_cid_str] = IntegerVector::create(cl_tracker.at(lm - 1), cl_tracker.at(rm - 1)); 
        
        // Record the distance the clusters appeared before the merge
        eps_birth[l_index] = eps_dist.at(lm - 1);
        eps_birth[r_index] = eps_dist.at(rm - 1);

        // Record integer IDs of containing points & merge steps
        contains[cl_cid_str] = combine(l_contains, r_contains);
        contains_k[cl_cid_str] = combine(l_ci, r_ci);; 
        eps_death[cl_cid_str] = eps_dist.at(k);
      } else {
        // Small correction: Two cluster IDs assigned to noise points that were merged to form a single 
        // cluster. Correct by choosing one of the IDs and replacing the other. This is necessary because of
        // information not known from computing the cluster agglomeratively, rather than divisively
        for (int i = 0, right = cl_tracker.at(rm - 1), n=cl_tracker.length(); i < n; ++i){
          if (cl_tracker[i] == right) { cl_tracker[i] = cl_tracker.at(lm - 1); }
        }
        cl_tracker.at(k) = cl_tracker.at(lm - 1);
        
        // Keep left index
        contains[l_index] = combine(l_contains, r_contains);
        contains_k[l_index] = combine(l_ci, r_ci);
        eps_death[l_index] = std::min((double) eps_death[l_index], (double) eps_death[r_index]);
        
        // Remove right index to cut down unnecessary memory usage
        contains[r_index] = R_NilValue;
        contains_k[r_index] = R_NilValue;
      }
    }
  }
  eps_birth[patch::to_string(cl_cid)] = eps_dist.at(k-1); // root node
  
  // Post processing 1: Compute epsilon values from recorded data, and from there initial stability scores, etc.
  NumericVector cl_eps; 
  IntegerVector cl_contains, cl_k;
  StringVector cl_ids = as<StringVector>(contains.names());
  for (int i = 0; i < cl_ids.size(); ++i){
    std::string cluster_id = Rcpp::as< std::string >(cl_ids.at(i));
    if (contains[cluster_id] != R_NilValue){
      // Get cluster point indices and merge steps used
      cl_contains = contains[cluster_id], cl_k = contains_k[cluster_id];
      
      // Retrieve the largest distance the cluster formed (going down the hierarchy) 
      double cl_birth = (double) eps_birth[cluster_id], cl_death = (double) eps_death[cluster_id];
    
      // Fill in eps values - never record distances below eps_death
      cl_eps = NumericVector(cl_contains.size());
      for (int i = 0, cl_size = cl_contains.size(); i < cl_size; ++i){
        cl_eps.at(i) = std::max(cl_death, (double) eps_dist.at(int(cl_k[i])));
      }
      
      // Per cluster information 
      NumericVector stability = (1/cl_eps) - (1/cl_birth); 
      clusters[cluster_id] = List::create(
        _["eps"] = cl_eps,
        _["eps_birth"] = cl_birth, 
        _["eps_death"] = cl_death,
        _["contains"] = cl_contains, 
        _["score"] = (double) sum(stability)
      );
    }
  }
  
  List res = List::create(
    _["member_sizes"] = member_sizes, 
    _["contains"] = contains, 
    _["clusters"] = clusters, 
    _["cl_tracker"] = cl_tracker, 
    _["cl_hierarchy"] = cl_hierarchy, 
    _["merge_steps"] = contains_k, 
    _["names"] = cl_ids, 
    _["eps_death"] = eps_death);
  
  // Post processing 2: Compute Salient Clusters
  std::string root_node = patch::to_string(cl_tracker[cl_tracker.length() - 1]);
  std::list<int>* sc = new std::list<int>();
  res[".root_node"] = stoi(root_node);
  computeSalientScores(res, root_node, *sc);
  res["salient_clusters"] = wrap(*sc);
  
  // Post processing 3: Unfold cluster assignments
  IntegerVector cluster = IntegerVector(n, 0);
  for (std::list<int>::iterator it = sc->begin(); it != sc->end(); it++){
    List cl = clusters[patch::to_string(*it)];
    IntegerVector cl_contains = cl["contains"];
    cluster[(cl_contains - 1)] = *it;
  }
  
  res["cluster"] = cluster; 
  return(res);
}
