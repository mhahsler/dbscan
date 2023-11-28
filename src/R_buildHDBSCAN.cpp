//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// C++ includes
#include <unordered_map>
#include <stack>
#include <queue>
#include <string> // std::atoi

// Helper functions
#include "utilities.h"

// Macros
#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

// Given a dist vector of "should-link" (1), "should-not-link" (-1), and "don't care" (0)
// constraints in the form of integers, convert constraints to a more compact adjacency list
// representation.
// [[Rcpp::export]]
List distToAdjacency(IntegerVector constraints, const int N){
  std::unordered_map<int, std::vector<int> > key_map = std::unordered_map<int, std::vector<int> >();
  for (int i = 0; i < N; ++i){
    for (int j = 0; j < N; ++j){
      if (i == j) continue;
      int index = i > j ? INDEX_TF(N, j, i) : INDEX_TF(N, i, j);
      int crule = constraints.at(index);
      if (crule != 0){
        if (key_map.count(i+1) != 1){ key_map[i+1] = std::vector<int>(); } // add 1 for base 1
        key_map[i+1].push_back(crule < 0 ? - (j + 1) : j + 1); // add 1 for base 1
      }
    }
  }
  return(wrap(key_map));
}

// Given an hclust object, convert to a dendrogram object (but much faster).
// [[Rcpp::export]]
List buildDendrogram(List hcl) {

  // Extract hclust info
  IntegerMatrix merge = hcl["merge"];
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
      left.attr("members") = (int) 1;
      left.attr("height") = (double) 0.f;
      left.attr("label") = labels.at(-(lm + 1));
      left.attr("leaf") = true;

      // Right
      IntegerVector right = IntegerVector::create(-rm);
      right.attr("members") = (int) 1;
      right.attr("height") = (double) 0.f;
      right.attr("label") = labels.at(-(rm + 1));
      right.attr("leaf") = true;

      // Merge
      new_br = List::create(left, right);
      new_br.attr("members") = 2;
      new_br.attr("midpoint") = 0.5;
    }
    // Second case: 1 is a singleton, the other is a branch
    else if (any(m < 0).is_true()){
      bool isL = lm < 0;

      // Create the leaf from the negative entry
      IntegerVector leaf = IntegerVector::create(isL ? -lm : -rm);
      leaf.attr("members") = 1;
      leaf.attr("height") = 0;
      leaf.attr("label") = labels.at(isL ? -(lm + 1) : -(rm + 1));
      leaf.attr("leaf") = true;

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

// Simple function to iteratively get the sub-children of a nested integer-hierarchy
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

// Extract 'flat' assignments
IntegerVector getSalientAssignments(List cl_tree, List cl_hierarchy, std::list<int> sc, const int n){
  IntegerVector cluster = IntegerVector(n, 0);
  for (std::list<int>::iterator it = sc.begin(); it != sc.end(); it++) {
    IntegerVector child_cl = all_children(cl_hierarchy, *it);

    // If at a leaf, its not necessary to recursively get point indices, else need to traverse hierarchy
    if (child_cl.length() == 0){
      List cl = cl_tree[patch::to_string(*it)];
      cluster[as<IntegerVector>(cl["contains"]) - 1] = *it;
    } else {
      List cl = cl_tree[patch::to_string(*it)];
      cluster[as<IntegerVector>(cl["contains"]) - 1] = *it;
      for (IntegerVector::iterator child_cid = child_cl.begin(); child_cid != child_cl.end(); ++child_cid){
        cl = cl_tree[patch::to_string(*child_cid)];
        IntegerVector child_contains = as<IntegerVector>(cl["contains"]);
        if (child_contains.length() > 0){
          cluster[child_contains - 1] = *it;
        }
      }
    }
  }
  return(cluster);
}

// Retrieve node (x, y) positions in a cluster tree
// [[Rcpp::export]]
NumericMatrix node_xy(List cl_tree, List cl_hierarchy, int cid = 0){

  // Initialize
  if (cid == 0){
    cl_tree["node_xy"] = NumericMatrix(all_children(cl_hierarchy, 0).size()+1, 2);
    cl_tree["leaf_counter"] = 0;
    cl_tree["row_counter"] = 0;
  }

  // Retrieve/set variables
  std::string cid_str = patch::to_string(cid);
  NumericMatrix node_xy_ = cl_tree["node_xy"];
  List cl = cl_tree[cid_str];

  // Increment row index every time
  int row_index = (int) cl_tree["row_counter"];
  cl_tree["row_counter"] = row_index+1;

  // base case
  if (!cl_hierarchy.containsElementNamed(cid_str.c_str())){
    int leaf_index = (int) cl_tree["leaf_counter"];
    node_xy_(row_index, _) = NumericVector::create((double) ++leaf_index, (double) cl["eps_death"]);
    cl_tree["leaf_counter"] = leaf_index;
    NumericMatrix res = NumericMatrix(1, 1);
    res[0] = row_index;
    return(res);
  } else {
    IntegerVector children = cl_hierarchy[cid_str];
    int l_row = (int) node_xy(cl_tree, cl_hierarchy, children.at(0))[0]; // left
    int r_row = (int) node_xy(cl_tree, cl_hierarchy, children.at(1))[0]; // right
    double lvalue = (double) (node_xy_(l_row, 0) + node_xy_(r_row, 0)) / 2;
    node_xy_(row_index, _) = NumericVector::create(lvalue, (double) cl["eps_death"]);

    if (cid != 0){
      NumericMatrix res = NumericMatrix(1, 1);
      res[0] = row_index;
      return(res);
    }
  }

  // Cleanup
  if (cid == 0){
    cl_tree["leaf_counter"] = R_NilValue;
    cl_tree["row_counter"] = R_NilValue;
  }
  return (node_xy_);
}

// Given a cluster tree, convert to a simplified dendrogram
// [[Rcpp::export]]
List simplifiedTree(List cl_tree) {

  // Hierarchical information
  List cl_hierarchy = cl_tree.attr("cl_hierarchy");
  IntegerVector all_childs = all_children(cl_hierarchy, 0);

  // To keep track of members and midpoints
  std::unordered_map<std::string, int> members = std::unordered_map<std::string, int>();
  std::unordered_map<std::string, float> mids = std::unordered_map<std::string, float>();

  // To keep track of where we are
  std::stack<int> cid_stack = std::stack<int>();
  cid_stack.push(0);

  // Iteratively build the hierarchy
  List dendrogram = List();

  // Premake children
  for (IntegerVector::iterator it = all_childs.begin(); it != all_childs.end(); ++it){
    std::string cid_label = patch::to_string(*it);
    List cl = cl_tree[cid_label];
    if (!cl_hierarchy.containsElementNamed(cid_label.c_str())){
      // Create leaf
      IntegerVector leaf = IntegerVector::create(*it);
      leaf.attr("label") = cid_label;
      leaf.attr("members") = 1;
      leaf.attr("height") = cl["eps_death"];
      leaf.attr("midpoint") = 0;
      leaf.attr("leaf") = true;
      dendrogram[cid_label] = leaf;
      members[cid_label] = 1;
      mids[cid_label] = 0;
    }
  }

  // Building the dendrogram bottom-up
  while(!cid_stack.empty()) {
    int cid = cid_stack.top();
    std::string cid_label = patch::to_string(cid);
    List cl = cl_tree[cid_label];

    // Recursive calls
    IntegerVector local_children = cl_hierarchy[cid_label];

    // Members and midpoint extraction
    std::string l_str = patch::to_string(local_children.at(0)), r_str = patch::to_string(local_children.at(1));
    // Rcout << "Comparing: " << l_str << ", " << r_str << std::endl;
    if (!dendrogram.containsElementNamed(l_str.c_str())){ cid_stack.push(local_children.at(0)); continue; }
    if (!dendrogram.containsElementNamed(r_str.c_str())){ cid_stack.push(local_children.at(1)); continue; }

    // Continue building up the hierarchy
    List left = dendrogram[l_str], right = dendrogram[r_str];

    int l_members = members[l_str], r_members = members[r_str];
    float l_mid = mids[l_str], r_mid = mids[r_str];

    // Make the new branch
    List new_branch = List::create(dendrogram[l_str], dendrogram[r_str]);
    new_branch.attr("label") = cid_label;
    new_branch.attr("members") = l_members + r_members;
    new_branch.attr("height") = (float) cl["eps_death"];
    new_branch.attr("class") = "dendrogram";

    // Midpoint calculation
    bool isL = (bool) !cl_hierarchy.containsElementNamed(l_str.c_str()); // is left a leaf
    if (!isL && cl_hierarchy.containsElementNamed(r_str.c_str())){ // is non-singleton merge
      new_branch.attr("midpoint") = (l_members + l_mid + r_mid) / 2;
    } else { // contains a leaf
      int sub_members = isL ? r_members : l_members;
      float mid_pt = isL ? r_mid : l_mid;
      new_branch.attr("midpoint") = ((isL ? 1 : sub_members) + mid_pt) / 2;
    }

    // Save info for later
    members[cid_label] = l_members + r_members;
    mids[cid_label] = (float) new_branch.attr("midpoint");
    dendrogram[cid_label] = new_branch;

    // Done with this node
    cid_stack.pop();
  }
  return(dendrogram["0"]);
}

/* Main processing step to compute all the relevent information in the form of the
 * 'cluster tree' for FOSC. The cluster stability scores are computed via the tree traversal rely on a separate function
 * Requires information associated with hclust elements. See ?hclust in R for more info.
 * 1. merge := an (n-1) x d matrix representing the MST computed from any arbitrary similarity matrix
 * 2. height := the (linkage) distance each new set of clusters forms from the MST
 * 3. order := the point indices of the original data the negative entries in merge refer to
 * Notation: eps is used to arbitrarily refer to the dissimilarity distance used
*/
// [[Rcpp::export]]
List computeStability(const List hcl, const int minPts, bool compute_glosh = false){
  // Extract hclust info
  NumericMatrix merge = hcl["merge"];
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"];
  int n = merge.nrow() + 1, k;

  //  Which cluster does each merge step represent (after the merge, or before the split)
  IntegerVector cl_tracker = IntegerVector(n-1 , 0),
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
  std::unordered_map<std::string, bool> processed = std::unordered_map<std::string, bool>();

  // First pass: Agglomerate up the hierarchy, recording member sizes.
  // This enables a dynamic programming strategy to improve performance below.
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

  // Initialize root (unknown size, might be 0, so don't initialize length)
  std::string root_str = "0";
  contains[root_str] = NumericVector();
  eps[root_str] = NumericVector();
  eps_birth[root_str] = eps_dist.at(eps_dist.length()-1);

  int global_cid = 0;
  // Second pass: Divisively split the hierarchy, recording the epsilon and point index values as needed
  for (k = n-2; k >= 0; --k){
    // Current Merge
    int lm = merge(k, 0), rm = merge(k, 1), cid = cl_tracker.at(k);
    IntegerVector m = IntegerVector::create(lm, rm);
    std::string cl_cid = patch::to_string(cid);

    // Trivial case: split into singletons, record eps, contains, and ensure eps_death is minimal
    if (all(m < 0).is_true()){
      contains[cl_cid].push_back(-lm), contains[cl_cid].push_back(-rm);
      double noise_eps = processed[cl_cid] ? eps_death[cl_cid] : eps_dist.at(k);
      eps[cl_cid].push_back(noise_eps), eps[cl_cid].push_back(noise_eps);
      eps_death[cl_cid] = processed[cl_cid] ? eps_death[cl_cid] : std::min((double) eps_dist.at(k), (double) eps_death[cl_cid]);
    } else if (any(m < 0).is_true()) {
      // Record new point info and mark the non-singleton with the cluster id
      contains[cl_cid].push_back(-(lm < 0 ? lm : rm));
      eps[cl_cid].push_back(processed[cl_cid] ? eps_death[cl_cid] : eps_dist.at(k));
      cl_tracker.at((lm < 0 ? rm : lm) - 1) = cid;
    } else {
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];

      // The minPts step
      if (merge_size1 >= minPts && merge_size2 >= minPts){
        // Record death of current cluster
        eps_death[cl_cid] = eps_dist.at(k);
        processed[cl_cid] = true;

        // Mark the lower merge steps as new clusters
        cl_hierarchy[cl_cid] = IntegerVector::create(global_cid+1, global_cid+2);
        std::string l_index = patch::to_string(global_cid+1), r_index = patch::to_string(global_cid+2);
        cl_tracker.at(lm - 1) = ++global_cid, cl_tracker.at(rm - 1) = ++global_cid;

        // Record the distance the new clusters appeared and initialize containers
        contains[l_index] = IntegerVector(), contains[r_index] = IntegerVector();
        eps[l_index] = NumericVector(), eps[r_index] = NumericVector(); ;
        eps_birth[l_index] = eps_dist.at(k), eps_birth[r_index] = eps_dist.at(k);
        eps_death[l_index] = eps_dist.at(lm - 1), eps_death[r_index] = eps_dist.at(rm - 1);
        processed[l_index] = false, processed[r_index] = false;
        n_children[cl_cid] = merge_size1 + merge_size2;
      } else {
        // Inherit cluster identity
        cl_tracker.at(lm - 1) = cid,  cl_tracker.at(rm - 1) = cid;
      }
    }
  }

  // Aggregate data into a returnable list
  // NOTE: the 'contains' element will be empty for all inner nodes w/ minPts == 1, else
  // it will contain only the objects that were considered 'noise' at that hierarchical level
  List res = List();
  NumericVector outlier_scores;
  if (compute_glosh) { outlier_scores = NumericVector( n, -1.0); }
  for (std::unordered_map<std::string, IntegerVector>::iterator key = contains.begin(); key != contains.end(); ++key){
    int nc = n_children[key->first];
    res[key->first] = List::create(
      _["contains"] = key->second,
      _["eps"] = eps[key->first],
      _["eps_birth"] = eps_birth[key->first],
      _["eps_death"] = eps_death[key->first],
      _["stability"] = sum(1/eps[key->first] - 1/eps_birth[key->first]) + (nc * 1/eps_death[key->first] - nc * 1/eps_birth[key->first]),
      //_["_stability"] = 1/eps[key->first] - 1/eps_birth[key->first],
      _["n_children"] = n_children[key->first]
    );

    // Compute GLOSH outlier scores (HDBSCAN only)
    if (compute_glosh){
      if (eps[key->first].size() > 0){ // contains noise points
        double eps_max = std::numeric_limits<double>::infinity();
        IntegerVector leaf_membership = all_children(cl_hierarchy, atoi(key->first.c_str()), true);
        if (leaf_membership.length() == 0){ // is itself a leaf
          eps_max = eps_death[key->first];
        } else {
          for (IntegerVector::iterator it = leaf_membership.begin(); it != leaf_membership.end(); ++it){
            eps_max = std::min(eps_max, eps_death[patch::to_string(*it)]);
          }
        }
        NumericVector eps_max_vec =  NumericVector(eps[key->first].size(), eps_max) / as<NumericVector>(eps[key->first]);
        NumericVector glosh = Rcpp::rep(1.0, key->second.length()) - eps_max_vec;
        outlier_scores[key->second - 1] = glosh;
      }
        // MFH: If the point is never an outlier (0/0) then set GLOSH to 0
        outlier_scores[is_nan(outlier_scores)] = 0.0;
    }
  }

  // Store meta-data as attributes
  res.attr("n") = n; // number of points in the original data
  res.attr("cl_hierarchy") = cl_hierarchy;  // Stores parent/child structure
  res.attr("cl_tracker") = cl_tracker; // stores cluster id formation for each merge step, used for cluster extraction
  res.attr("minPts") = minPts; // needed later
  // res.attr("root") = minPts == 1; // needed later to ensure root is not captured as a cluster
  if (compute_glosh){ res.attr("glosh") = outlier_scores; } // glosh outlier scores (hdbscan only)
  return(res);
}

// Validates a given list of instance-level constraints for symmetry. Since the number of
// constraints might change dramatically based on the problem, and initial loop is performed
// to figure out whether it would be faster to check via an adjacencty list or matrix
// [[Rcpp::export]]
List validateConstraintList(List& constraints, int n){
  std::vector< std::string > keys = as< std::vector< std::string > >(constraints.names());
  bool is_valid = true, tmp_valid, use_matrix = false;

  int n_constraints = 0;
  for (List::iterator it = constraints.begin(); it != constraints.end(); ++it){
    n_constraints += as<IntegerVector>(*it).size();
  }

  // Sparsity check: if the constraints make up a sufficiently large amount of
  // the solution space, use matrix to check validity
  if (n_constraints/(n*n) > 0.20){
    use_matrix = true;
  }

  // Check using adjacency matrix
  if (use_matrix){
    IntegerMatrix adj_matrix = IntegerMatrix(Dimension(n, n));
    int from, to;
    for (std::vector< std::string >::iterator it = keys.begin(); it != keys.end(); ++it){
      // Get constraints
      int cid = atoi(it->c_str()); // to base-0
      IntegerVector cs_ = constraints[*it];

      // Positive "should-link" constraints
      IntegerVector pcons = as<IntegerVector>(cs_[cs_ > 0]);
      for (IntegerVector::iterator pc = pcons.begin(); pc != pcons.end(); ++pc){
        from = (*pc < cid ? *pc : cid) - 1;
        to = (*pc > cid ? *pc : cid) - 1;
        adj_matrix(from, to) = 1;
      }

      // Negative "should-not-link" constraints
      IntegerVector ncons = -(as<IntegerVector>(cs_[cs_ < 0]));
      for (IntegerVector::iterator nc = ncons.begin(); nc != ncons.end(); ++nc){
        from = (*nc < cid ? *nc : cid) - 1;
        to = (*nc > cid ? *nc : cid) - 1;
        adj_matrix(from, to) = -1;
      }
    }

    // Check symmetry
    IntegerVector lower = lowerTri(adj_matrix);
    IntegerMatrix adj_t = Rcpp::transpose(adj_matrix);
    IntegerVector lower_t = lowerTri(adj_t);
    LogicalVector valid_check = lower == lower_t;
    is_valid = all(valid_check == TRUE).is_true();

    // Try to merge the two
    if (!is_valid){
      int sum = 0;
      for (int i = 0; i < lower.size(); ++i){
        sum = lower.at(i) + lower_t.at(i);
        lower[i] = sum > 0 ? 1 : sum < 0 ? -1 : 0;
      }
    }
    constraints = distToAdjacency(lower, n);
  }
  // Else check using given adjacency list
  else {
    for (std::vector< std::string >::iterator it = keys.begin(); it != keys.end(); ++it){
      // Get constraints
      int cid = atoi(it->c_str());
      IntegerVector cs_ = constraints[*it];

      // Positive "should-link" constraints
      IntegerVector pcons = as<IntegerVector>(cs_[cs_ > 0]);
      for (IntegerVector::iterator pc = pcons.begin(); pc != pcons.end(); ++pc){
        int ic = *pc < 0 ? -(*pc) : *pc;
        std::string ic_str = patch::to_string(ic);
        bool exists = constraints.containsElementNamed(ic_str.c_str());
        tmp_valid = exists ? contains(as<IntegerVector>(constraints[ic_str]), cid) : false;
        if (!tmp_valid){
          if (!exists){
            constraints[ic_str] = IntegerVector::create(cid);
          } else {
            IntegerVector con_vec = constraints[ic_str];
            con_vec.push_back(cid);
            constraints[ic_str] = con_vec;
          }
          is_valid = false;
        }
      }

      // Negative "should-not-link" constraints
      IntegerVector ncons = -(as<IntegerVector>(cs_[cs_ < 0]));
      for (IntegerVector::iterator nc = ncons.begin(); nc != ncons.end(); ++nc){
        int ic = *nc < 0 ? -(*nc) : *nc;
        std::string ic_str = patch::to_string(ic);
        bool exists = constraints.containsElementNamed(ic_str.c_str());
        tmp_valid = exists ? contains(as<IntegerVector>(constraints[ic_str]), cid) : false;
        if (!tmp_valid){
          if (!exists){
            constraints[ic_str] = IntegerVector::create(-cid);
          } else {
            IntegerVector con_vec = constraints[ic_str];
            con_vec.push_back(-cid);
            constraints[ic_str] = con_vec;
          }
          is_valid = false;
        }
      }
    }
  }
  // Produce warning if asymmetric constraints detected; return attempt at fixing constraints.
  if (!is_valid){
      warning("Incomplete (asymmetric) constraints detected. Populating constraint list.");
  }
  return(constraints);
}

// [[Rcpp::export]]
double computeVirtualNode(IntegerVector noise, List constraints){
  if (noise.length() == 0) return(0);
  if (Rf_isNull(constraints)) return(0);

  // Semi-supervised extraction
  int satisfied_constraints = 0;
  // Rcout << "Starting constraint based optimization" << std::endl;
  for (IntegerVector::iterator it = noise.begin(); it != noise.end(); ++it){
    std::string cs_str = patch::to_string(*it);
    if (constraints.containsElementNamed(cs_str.c_str())){
      // Get constraints
      IntegerVector cs_ = constraints[cs_str];

      // Positive "should-link" constraints
      IntegerVector pcons = as<IntegerVector>(cs_[cs_ > 0]);
      for (IntegerVector::iterator pc = pcons.begin(); pc != pcons.end(); ++pc){
        satisfied_constraints += contains(noise, *pc);
      }

      // Negative "should-not-link" constraints
      IntegerVector ncons = -(as<IntegerVector>(cs_[cs_ < 0]));
      for (IntegerVector::iterator nc = ncons.begin(); nc != ncons.end(); ++nc){
        satisfied_constraints += (1 - contains(noise, *nc));
      }
    }
  }
  return(satisfied_constraints);
}


// Framework for Optimal Selection of Clusters (FOSC)
// Traverses a cluster tree hierarchy to compute a flat solution, maximizing the:
// - Unsupervised soln: the 'most stable' clusters following the give linkage criterion
// - SS soln w/ instance level Constraints: constraint-based w/ unsupervised tiebreaker
// - SS soln w/ mixed objective function: maximizes J = α JU + (1 − α) JSS
// [[Rcpp::export]]
NumericVector fosc(List cl_tree, std::string cid, std::list<int>& sc, List cl_hierarchy,
                   bool prune_unstable_leaves=false, // whether to prune -very- unstable subbranches
                   const double alpha = 0, // mixed objective case
                   bool useVirtual = false, // return virtual node as well
                   const int n_constraints = 0, // number of constraints
                   List constraints = R_NilValue) // instance-level constraints
{
  // Base case: at a leaf
  if (!cl_hierarchy.containsElementNamed(cid.c_str())){
    List cl = cl_tree[cid];
    sc.push_back(std::atoi(cid.c_str())); // assume the leaf will be a salient cluster until proven otherwise
    return(NumericVector::create((double) cl["stability"],
                                 (double) useVirtual ? cl["vscore"] : 0));
  } else {
    // Non-base case: at a merge of clusters, determine which to keep
    List cl = cl_tree[cid];

    // Get child stability/constraint scores
    NumericVector scores, stability_scores = NumericVector(), constraint_scores = NumericVector();
    IntegerVector child_ids = cl_hierarchy[cid];
    for (int i = 0, clen = child_ids.length(); i < clen; ++i){
      int child_id = child_ids.at(i);
      scores = fosc(cl_tree, patch::to_string(child_id), sc, cl_hierarchy, prune_unstable_leaves, alpha, useVirtual, n_constraints, constraints);
      stability_scores.push_back(scores.at(0));
      constraint_scores.push_back(scores.at(1));
    }

    // If semisupervised scenario, normalizing should be stored in 'total_stability'
    double total_stability = (contains(cl_tree.attributeNames(),"total_stability") ? (double) cl_tree.attr("total_stability") : 1.0);

    // Compare and update stability scores
    double old_stability_score = (double) cl["stability"] / total_stability;
    double new_stability_score = (double) sum(stability_scores) / total_stability;

    // Compute instance-level constraints if necessary
    double old_constraint_score = 0, new_constraint_score = 0;
    if (useVirtual){
      // Rcout << "old constraint score for " << cid << ": " << (double) cl["vscore"] << std::endl;
      old_constraint_score = (double) cl["vscore"];
      new_constraint_score = (double) sum(constraint_scores) + (double) computeVirtualNode(cl["contains"], constraints)/n_constraints;
    }

    bool keep_children = true;
    // If the score is unchanged, remove the children and add parent
    if (useVirtual){
      if (old_constraint_score < new_constraint_score && cid != "0"){
        // Children satisfies more constraints
        cl["vscore"] = new_constraint_score;
        cl["score"] = alpha * new_stability_score + (1 - alpha) * new_constraint_score;
        // Rcout << "1: score for " << cid << ":" << (double) cl["score"] << std::endl;
        // Rcout << "(old constraint): " << old_constraint_score << ", (new constraint): " << new_constraint_score << std::endl;
      } else if (old_constraint_score > new_constraint_score && cid != "0"){
        // Parent satisfies more constraints
        cl["vscore"] = old_constraint_score;
        cl["score"] = alpha * old_stability_score + (1 - alpha) * old_constraint_score;
        // Rcout << "2: score for " << cid << ":" << (double) cl["score"] << std::endl;
        keep_children = false;
      } else {
        // Resolve tie using unsupervised, stability-based approach
        if (old_stability_score < new_stability_score){
          // Children are more stable
          cl["score"] = new_stability_score / total_stability;
          // Rcout << "3: score for " << cid << ":" << (double) cl["score"] << std::endl;
        } else {
          // Parent is more stable
          cl["score"] = old_stability_score / total_stability;
          // Rcout << "4: score for " << cid << ":" << (double) cl["score"] << std::endl;
          // Rcout << "(old stability): " << old_stability_score << ", (total stability): " << total_stability << std::endl;
          keep_children = false;
        }
        cl["vscore"] = old_constraint_score;
      }
    } else {
      // Use unsupervised, stability-based approach only
      if (old_stability_score < new_stability_score){
        cl["score"] = new_stability_score; // keep children
      } else {
        cl["score"] = old_stability_score;
        keep_children = false;
      }
    }


    // Prune children and add parent (cid) if need be
    if (!keep_children && cid != "0") {
      IntegerVector children = all_children(cl_hierarchy, std::atoi(cid.c_str())); // use all_children to prune subtrees
      for (int i = 0, clen = children.length(); i < clen; ++i){
        sc.remove(children.at(i)); // use list for slightly better random deletion performance
      }
      sc.push_back(std::atoi(cid.c_str()));
    } else if (keep_children && prune_unstable_leaves){
      // If flag passed, prunes leaves with insignificant stability scores
      // this can happen in cases where one leaf has a stability score significantly greater
      // than both its siblings and its parent (or other ancestors), causing sibling branches
      // to be considered as clusters even though they may nto be significantly more stable than their parent
      if (all(stability_scores < old_stability_score).is_false()){
        for (int i = 0, clen = child_ids.length(); i < clen; ++i){
          if (stability_scores.at(i) < old_stability_score){
            IntegerVector to_prune = all_children(cl_hierarchy, child_ids.at(i)); // all sub members
            for (IntegerVector::iterator it = to_prune.begin(); it != to_prune.end(); ++it){
              //Rcout << "Pruning: " << *it << std::endl;
              sc.remove(*it);
            }
          }
        }
      }
    }

    // Save scores for traversal up and for later
    cl_tree[cid] = cl;

    // Return this sub trees score
    return(NumericVector::create((double) cl["score"], useVirtual ? (double) cl["vscore"] : 0));
  }
}

// Given a cluster tree object with computed stability precomputed scores from computeStability,
// extract the 'most stable' or salient flat cluster assignments. The large number of derivable
// arguments due to fosc being a recursive function
// [[Rcpp::export]]
List extractUnsupervised(List cl_tree, bool prune_unstable = false){
  // Compute Salient Clusters
  std::list<int> sc = std::list<int>();
  List cl_hierarchy = cl_tree.attr("cl_hierarchy");
  int n = as<int>(cl_tree.attr("n"));
  fosc(cl_tree, "0", sc, cl_hierarchy, prune_unstable); // Assume root node is always id == 0

  // Store results as attributes
  cl_tree.attr("cluster") = getSalientAssignments(cl_tree, cl_hierarchy, sc, n); // Flat assignments
  cl_tree.attr("salient_clusters") = wrap(sc); // salient clusters
  return(cl_tree);
}

// [[Rcpp::export]]
List extractSemiSupervised(List cl_tree, List constraints, float alpha = 0, bool prune_unstable_leaves = false){
  // Rcout << "Starting semisupervised extraction..." << std::endl;
  List root = cl_tree["0"];
  List cl_hierarchy = cl_tree.attr("cl_hierarchy");
  int n = as<int>(cl_tree.attr("n"));

  // Compute total number of constraints
  int n_constraints = 0;
  for (int i = 0, n = constraints.length(); i < n; ++i){
    IntegerVector cl_constraints = constraints.at(i);
    n_constraints += cl_constraints.length();
  }

  // Initialize root
  List cl = cl_tree["0"];
  cl["vscore"] = 0;
  cl_tree["0"] = cl; // replace to keep changes

  // Compute initial gamma values or "virtual nodes" for both leaf and internal nodes
  IntegerVector cl_ids = all_children(cl_hierarchy, 0);
  for (IntegerVector::iterator it = cl_ids.begin(); it != cl_ids.end(); ++it){
    if (*it != 0){
      std::string cid_str = patch::to_string(*it);
      List cl = cl_tree[cid_str];

      // Store the initial fraction of constraints satisfied for each node as 'vscore'
      // NOTE: leaf scores represent \hat{gamma}, internal represent virtual node scores
      if (cl_hierarchy.containsElementNamed(cid_str.c_str())){
        // Extract the point indices the cluster contains
        IntegerVector child_cl = all_children(cl_hierarchy, *it), child_ids;
        List cl_container = List();
        for (IntegerVector::iterator ch_id = child_cl.begin(); ch_id != child_cl.end(); ++ch_id){
          List ch_cl = cl_tree[patch::to_string(*ch_id)];
          //child_ids = combine(child_ids, ch_cl["contains"]);
          cl_container.push_back(as<IntegerVector>(ch_cl["contains"]));
        }
        cl_container.push_back(as<IntegerVector>(cl["contains"]));
        child_ids = concat_int(cl_container);
        cl["vscore"] = computeVirtualNode(child_ids, constraints)/n_constraints;
      } else { // is leaf node
        cl["vscore"] = computeVirtualNode(cl["contains"], constraints)/n_constraints;
      }
      cl_tree[cid_str] = cl; // replace to keep changes
    }
  }

  // First pass: compute unsupervised soln as a means of extracting normalizing constant J_U^*
  cl_tree = extractUnsupervised(cl_tree, false);
  IntegerVector stable_sc = cl_tree.attr("salient_clusters");
  double total_stability = 0.0f;
  for (IntegerVector::iterator it = stable_sc.begin(); it != stable_sc.end(); ++it){
    List cl = cl_tree[patch::to_string(*it)];
    total_stability += (double) cl["stability"];
  }
  cl_tree.attr("total_stability") = total_stability;
  // Rcout << "Total stability: " << total_stability << std::endl;

  // Compute stable clusters w/ instance-level constraints
  std::list<int> sc = std::list<int>();
  fosc(cl_tree, "0", sc, cl_hierarchy, prune_unstable_leaves,
       alpha, true, n_constraints, constraints); // semi-supervised parameters

  // Store results as attributes and return
  cl_tree.attr("salient_clusters") = wrap(sc);
  cl_tree.attr("cluster") = getSalientAssignments(cl_tree, cl_hierarchy, sc, n);
  return(cl_tree);
}




