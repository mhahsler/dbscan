#include <Rcpp.h>
using namespace Rcpp;


#include <unordered_map>
#include <stack>
#include <queue>

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
  List labels = List(); // allows to avoid type inference 
  if (!hcl.containsElementNamed("labels") || hcl["labels"] == R_NilValue){
    labels = seq_along(order); 
  } else { 
    labels = as<StringVector>(hcl["labels"]); 
  }
  

  int n = merge.nrow() + 1, k; 
  List new_br, z = List(n);
  for (k = 0; k < n-1; k++){ 
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm); 
    
    // First Case: Both are singletons, so need to create leaves
    if (all(m < 0).is_true()){
      // Left 
      IntegerVector left = IntegerVector::create(-lm);
      left.attr("label") = labels.at(-(lm + 1));
      left.attr("members") = 1;
      left.attr("height") = 0;
      left.attr("leaf") = true;

      // Right 
      IntegerVector right = IntegerVector::create(-rm);
      right.attr("label") = labels.at(-(rm + 1));
      right.attr("members") = 1;
      right.attr("height") = 0;      
      right.attr("leaf") = true;

      // Merge 
      new_br = List::create(left, right);
      new_br.attr("midpoint") = 0.5;
      new_br.attr("members") = 2;
      //new_br.attr("class") = "dendrogram";
    } 
    // Second case: 1 is a singleton, the other is a branch
    else if (any(m < 0).is_true()){
      bool isL = lm < 0;

      // Create the leaf from the negative entry
      IntegerVector leaf = IntegerVector::create(isL ? -lm : -rm);
      leaf.attr("members") = 1;
      leaf.attr("leaf") = true;
      leaf.attr("height") = 0;
      leaf.attr("label") = labels.at(isL ? -(lm + 1) : -(rm + 1));

      // Merge the leaf with the other existing branch
      int branch_key = isL ? rm - 1 : lm - 1;
      List sub_branch = z[branch_key];
      new_br = isL ? List::create(leaf, sub_branch) : List::create(sub_branch, leaf);
      z.at(branch_key) = R_NilValue;

      // Set attributes of new branch
      int sub_members = sub_branch.attr("members");
      double mid_pt = sub_branch.attr("midpoint");
      new_br.attr("members") = int(sub_members) + 1;
      new_br.attr("midpoint") = (int(isL ? 1 : sub_members) + mid_pt) / 2;
      //new_br.attr("class") = "dendrogram";
    } else {
      // Create the new branch
      List l_branch = z.at(lm - 1), r_branch = z.at(rm - 1);
      new_br = List::create(l_branch, r_branch);

      // Store attribute valeus in locals to get around proxy
      int left_members = l_branch.attr("members"), right_members = r_branch.attr("members");
      double l_mid = l_branch.attr("midpoint"), r_mid = r_branch.attr("midpoint");

      // Set up new branch attributes
      new_br.attr("members") = left_members + right_members;
      new_br.attr("midpoint") = (left_members + l_mid + r_mid) / 2;
      new_br.attr("class") = "dendrogram";

      // Deallocate unneeded memory along the way
      z.at(lm - 1) = R_NilValue;
      z.at(rm - 1) = R_NilValue;
    }
    new_br.attr("height") = height.at(k);
    z.at(k) = new_br;
  }
  List res = z.at(k - 1);
  res.attr("class") = "dendrogram";
  return(res);
}

// [[Rcpp::export]]
IntegerVector all_children(List hier, int key, bool leaves_only = false){
  IntegerVector res = IntegerVector(); 
  
  // If the key doesn't exist return an empty vector
  if (!hier.containsElementNamed(patch::to_string(key).c_str())){
    return(res);
  }
  
  // Else, do iterative 'recursive' type function to extract all the IDs of 
  // all sub trees
  IntegerVector children = hier[patch::to_string(key).c_str()];
  std::queue<int> to_do = std::queue<int>(); 
  to_do.push(key);
  while (to_do.size() != 0){
    int parent = to_do.front();
    if (!hier.containsElementNamed(patch::to_string(parent).c_str())){
      to_do.pop();
    } else {
      children = hier[patch::to_string(parent).c_str()];
      to_do.pop();
      for (int n_children = 0; n_children < children.length(); ++n_children){
        int child_id = children.at(n_children);
        if (leaves_only){
          if (!hier.containsElementNamed(patch::to_string(child_id).c_str())) {
            res.push_back(child_id);
          }
        } else { res.push_back(child_id); }
        to_do.push(child_id);
      }
    }
  }
  return(res);
}

// [[Rcpp::export]]
NumericMatrix node_xy(List hier){
  IntegerVector children = all_children(hier, 0);
  NumericMatrix res = NumericMatrix(children.length(), 2);

  Rcout << "here" << std::endl; 
  // Do iterative 'recursive' type function to extract all the IDs of sub trees
  std::queue<int> to_do = std::queue<int>();
  to_do.push(0);
  int leaf_counter = 0, row_index = 0;
  while (to_do.size() != 0){
    int parent = to_do.front();
    if (!hier.containsElementNamed(patch::to_string(parent).c_str())){
      leaf_counter++;
      row_index++;
      if (row_index < children.length()){
        res(row_index, _) = NumericVector::create(leaf_counter, (double) to_do.front());
      }
      to_do.pop();
    } else {
      IntegerVector local_ch = hier[patch::to_string(parent).c_str()];
      if (row_index < children.length())
        res(row_index, _) = NumericVector::create(leaf_counter, (double) to_do.front());
      to_do.pop();
      for (int n_children = 0; n_children < local_ch.length(); ++n_children){
        int child_id = local_ch.at(n_children);
        if (!hier.containsElementNamed(patch::to_string(child_id).c_str())) { // is leaf
          leaf_counter++;
        }
        row_index++;
        to_do.push(child_id);
      }
    }
    if (row_index < children.length()){
      Rcout << "todo: " << to_do.front() <<  std::endl;
      //res(row_index, _) = NumericVector::create(leaf_counter, (double) to_do.front());
    }
  }
  return(res);
}



// [[Rcpp::export]]
List buildCondensedTree(List hdbscan){
  //List clusters = hdbscan.attr("cluster");
  List cl_hierarchy = hdbscan.attr("cl_hierarchy");
  IntegerVector all_childs = all_children(cl_hierarchy, 0);
  NumericMatrix node_xy = NumericMatrix(all_childs.length()+1, 2); 
  //IntegerVector cl_tracker = info["cl_tracker"];
  //int root_node = info[".root_node"];

  IntegerVector children = cl_hierarchy["0"];
  
  List dendrogram = List(hdbscan.length());
  // to_do.push(root_node);

  List members = List(), midpoints = List();
  std::vector<std::string> cl_labels = (hdbscan.names());
  for (int i = 0, cid = 0; i < cl_labels.size(); ++i) {
    if (!cl_hierarchy.containsElementNamed(cl_labels.at(i).c_str())){
      IntegerVector leaf = IntegerVector::create(++cid); //clusters[patch::to_string(parent)]
      List cl = hdbscan[cl_labels.at(i).c_str()]; 
      leaf.attr("label") = cl_labels.at(i);
      leaf.attr("members") = 1;
      leaf.attr("height") = cl["eps_death"];
      leaf.attr("midpoint") = 0;
      leaf.attr("leaf") = true;
      leaf.attr(".info") = cl; 
      dendrogram[cl_labels.at(i)] = leaf;
      members[cl_labels.at(i)] = 1;
      midpoints[cl_labels.at(i)] = 0;
    }
  }

  // Build the hierarchy top-down 
  std::stack<int> to_do = std::stack<int>(); //(size_t) hdbscan.length()
  to_do.push(0);
  int row_index = 0, leaf_index = 0;

  while (to_do.size() != 0) {
    int parent = to_do.top();
    //Rcout << "Processing: " << parent << std::endl; 
    // if (!cl_hierarchy.containsElementNamed(patch::to_string(parent).c_str())){
    //   to_do.pop();
    //   Rcout << "here" << std::endl; 
    // } else {
      children = cl_hierarchy[patch::to_string(parent).c_str()];
      List child_branches = List(children.length());
      bool process = true;
      for (int c_i = 0, n_members = 0; c_i < children.length(); ++c_i){
        int child_id = children.at(c_i);
        //Rcout << "Pushing: " << child_id << " to slot " << c_i << std::endl; 
        if (dendrogram.containsElementNamed(patch::to_string(child_id).c_str())) {
          Rcout << "Pushing singleton: " << child_id << " i: " << leaf_index << std::endl; 
          child_branches.at(c_i) = dendrogram[patch::to_string(child_id)];
          leaf_index++;
        } else {
          process = false;
          Rcout << "Pushing non-singleton: " << child_id << " at ri: " << row_index << std::endl;
          to_do.push(child_id);
          row_index++; 
        }
      }
      if (process){
        to_do.pop();
        List parent_cl = hdbscan[patch::to_string(parent)];
        node_xy(row_index--, _) = NumericVector::create(row_index, parent_cl["eps_birth"]);
        //Rcout << "RI: " << row_index << std::endl; 
        
        std::string l_str = patch::to_string(children.at(0));
        std::string r_str = patch::to_string(children.at(1));
        List l_branch = dendrogram[l_str], r_branch = dendrogram[r_str];

        // Unfortunately this is necessary
        int l_members = members[l_str], r_members = members[r_str];
        float l_mid = midpoints[l_str], r_mid = midpoints[r_str];

        // Make the new branch
        child_branches.attr("label") = patch::to_string(parent);
        child_branches.attr("members") = l_members + r_members;

        bool isL = (bool) !cl_hierarchy.containsElementNamed(l_str.c_str()); // is left a leaf
        if (!isL && cl_hierarchy.containsElementNamed(r_str.c_str())){ // is non-singleton merge
          child_branches.attr("midpoint") = (l_members + l_mid + r_mid) / 2;
        } else { // contains a leaf
          //Rcout << "HERE";
          int sub_members = isL ? r_members : l_members;
          float mid_pt = isL ? r_mid : l_mid;
          child_branches.attr("midpoint") = ((isL ? 1 : sub_members) + mid_pt) / 2;
        }
        child_branches.attr("height") = (float) parent_cl["eps_birth"];
        child_branches.attr(".info") = parent_cl;
        child_branches.attr("class") = "dendrogram";

        // Need to store this for later
        dendrogram[patch::to_string(parent)] = child_branches;
        midpoints[patch::to_string(parent)] = (float) child_branches.attr("midpoint");
        members[patch::to_string(parent)] = (float) child_branches.attr("members");

        // Don't need these anymore
        dendrogram[l_str] = R_NilValue, dendrogram[r_str] = R_NilValue;
      }
    }
  // }
  List root = dendrogram["0"]; 
  root.attr("node_xy") = node_xy;
  return(root);
}

// Compute stability scores for cluster objects in the hierarchy
// [[Rcpp::export]]
double computeSalientScores(const List hdbscan, std::string cid, std::list<int>& sc, List cl_hierarchy){
  // Base case: at a leaf
  if (!cl_hierarchy.containsElementNamed(cid.c_str())){
    List cl = hdbscan[cid];
    sc.push_back(stoi(cid)); // assume the leaf will be a salient cluster until proven otherwise
    return((double) cl["score"]); // Use precomputed scores for leaves 
  } else {
  // Non-base case: at a merge of clusters, determine which to keep
    List cl = hdbscan[cid];
    
    // Get child stability scores
    NumericVector child_scores = NumericVector();
    IntegerVector child_ids = cl_hierarchy[cid];
    for (int i = 0, clen = child_ids.length(); i < clen; ++i){
      int child_id = child_ids.at(i);
      child_scores.push_back(computeSalientScores(hdbscan, patch::to_string(child_id), sc, cl_hierarchy));
    }
    
    // Compare and update stability scores 
    double old_score = (double) cl["score"];
    double new_score = std::max(old_score, (double) sum(child_scores));
    cl["new_score"] = new_score;
    
    // If the score is unchanged, remove the children and add parent
    if (new_score == old_score && cid != "0") {
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
List hdbscan_fast(const List hcl, const int minPts){
  // Extract hclust info
  NumericMatrix merge = hcl["merge"]; 
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"]; 
  int n = merge.nrow() + 1, k; 
  
  IntegerVector cl_tracker = IntegerVector(n-1 , 0), //  Which cluster does each merge step represent
                member_sizes = IntegerVector(n-1, 0); // Size each step
  
  List clusters = List(), // Final cluster information
       cl_hierarchy = List(); // Keeps track of hierarchy, which cluster contains who 

  // The primary information needed  
  std::unordered_map<std::string, IntegerVector> contains = std::unordered_map<std::string, IntegerVector>(); 
  std::unordered_map<std::string, NumericVector> eps = std::unordered_map<std::string, NumericVector>(); 
  
  // Supplemental information for either conveniance or to reduce memory
  std::unordered_map<std::string, int> n_children = std::unordered_map<std::string, int>(); 
  std::unordered_map<std::string, double> eps_death = std::unordered_map<std::string, double>(); 
  std::unordered_map<std::string, double> eps_birth = std::unordered_map<std::string, double>(); 
  
  // First pass: Agglomerate up the hierarchy, recording member sizes
  for (k = 0; k < n-1; ++k){
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm);
    if (all(m < 0).is_true()){
      member_sizes[k] = 2;
    } else if (any(m < 0).is_true()) {
      int pos_merge = (lm < 0 ? rm : lm), merge_size = member_sizes[pos_merge - 1];
      member_sizes[k] = merge_size + 1;
    } else {
      // Record Member Sizes
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];
      member_sizes[k] = merge_size1 + merge_size2;
    }
  }  
  
  // Initialize root 
  contains["0"] = pt_order; 
  eps["0"] = NumericVector(n, eps_dist.at(eps_dist.length()-1));   
  eps_birth["0"] = eps_dist.at(eps_dist.length()-1); 
  
  int global_cid = 0; 
  // Second pass: Divisively split the hierarchy, recording the epsilon and point index values as needed
  for (k = n-2; k >= 0; --k){
    // Current Merge
    int lm = merge(k, 0), rm = merge(k, 1), cid = cl_tracker.at(k);;
    IntegerVector m = IntegerVector::create(lm, rm);
    std::string cl_cid = patch::to_string(cid);
    // Rcout << "lm: " << lm << ", rm: " << rm  << ", cid:" << cid << std::endl; 
    
    // Trivial case: merge of singletons, create a temporary *assumed* cluster to be resolved on merger
    if (all(m < 0).is_true()){
      contains[cl_cid].push_back(-lm), contains[cl_cid].push_back(-rm);
      eps[cl_cid].push_back(eps_dist.at(k)), eps[cl_cid].push_back(eps_dist.at(k)); 
      eps_death[cl_cid] = std::min((double) eps_dist.at(k), (double) eps_death[cl_cid]); 
    } else if (any(m < 0).is_true()) {
      // Record new point info and mark the non-singleton with the cluster id
      contains[cl_cid].push_back(-(lm < 0 ? lm : rm)), eps[cl_cid].push_back(eps_dist.at(k));
      cl_tracker.at((lm < 0 ? rm : lm) - 1) = cid;
    } else {
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];

      // The HDBSCAN step 
      if (merge_size1 >= minPts && merge_size2 >= minPts){
        // Record death of current cluster
        eps_death[cl_cid] = eps_dist.at(k);
        
        // Mark the lower merge steps as new clusters 
        cl_hierarchy[cl_cid] = IntegerVector::create(global_cid+1, global_cid+2);
        std::string l_index = patch::to_string(global_cid+1), r_index = patch::to_string(global_cid+2);
        cl_tracker.at(lm - 1) = ++global_cid, cl_tracker.at(rm - 1) = ++global_cid; 
        
        // Record the distance the new clusters appeared and initialize containers
        contains[l_index] = IntegerVector(), contains[r_index] = IntegerVector(); ; 
        eps[l_index] = NumericVector(), eps[r_index] = NumericVector(); 
        eps_birth[l_index] = eps_dist.at(lm - 1),  eps_birth[r_index] = eps_dist.at(rm - 1);
        eps_death[l_index] = eps_dist.at(lm - 1),  eps_death[r_index] = eps_dist.at(lm - 1); 
        n_children[cl_cid] = merge_size1 + merge_size2; 
      } else {
        // Inherit cluster identity 
        cl_tracker.at(lm - 1) = cid,  cl_tracker.at(rm - 1) = cid; 
      }
    }
  }
  
  // Aggregate data into a returnable list 
  List res = List(); 
  NumericVector outlier_scores = NumericVector((int) n, -1.0);
  for (std::unordered_map<std::string, IntegerVector>::iterator key = contains.begin(); key != contains.end(); ++key){
    int nc = n_children[key->first]; 
    res[key->first] = List::create(
      _["contains"] = key->second, 
      _["eps"] = eps[key->first],
      _["eps_birth"] = eps_birth[key->first], 
      _["eps_death"] = eps_death[key->first], 
      _["score"] = sum(1/eps[key->first] - 1/eps_birth[key->first]) + (nc * 1/eps_death[key->first] - nc * 1/eps_birth[key->first]),
      _["n_children"] = n_children[key->first]            
    );
    
    // Compute GLOSH outlier scores over leaves  
    if (!cl_hierarchy.containsElementNamed((key->first).c_str())) {
      double lambda_max = 1/(eps_birth[key->first]);
      NumericVector glosh = NumericVector(key->second.length(), 1) - (lambda_max/eps[key->first]);
      outlier_scores[key->second - 1] = glosh;
    }
  }
  
  // Compute Salient Clusters
  std::list<int> sc = std::list<int>();
  computeSalientScores(res, "0", sc, cl_hierarchy);
  res.attr("salient_clusters") = wrap(sc);
  
  // Unfold cluster assignments
  IntegerVector cluster = IntegerVector(n, 0);
  for (std::list<int>::iterator it = sc.begin(); it != sc.end(); it++) {
    IntegerVector child_cl = all_children(cl_hierarchy, *it);

    // If at a leaf, use not necessary to recursively get point indices, else need to traverse hierarchy
    if (child_cl.length() == 0){
      List cl = res[patch::to_string(*it)]; 
      cluster[as<IntegerVector>(cl["contains"]) - 1] = *it;
    } else {
      List cl = res[patch::to_string(*it)]; 
      cluster[as<IntegerVector>(cl["contains"]) - 1] = *it; 
      for (IntegerVector::iterator child_cid = child_cl.begin(); child_cid != child_cl.end(); ++child_cid){
        cl = res[patch::to_string(*child_cid)];
        IntegerVector child_contains = as<IntegerVector>(cl["contains"]);
        if (child_contains.length() > 0){
          cluster[child_contains - 1] = *it;
        }
      }
    }
  }
  res.attr("cl_hierarchy") = cl_hierarchy;  // Stores parent/child structure 
  res.attr("cluster") = cluster; 
  res.attr("glosh") = outlier_scores; 
  return(res);
}




