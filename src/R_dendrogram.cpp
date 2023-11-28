//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#include <Rcpp.h>
#include <sstream>
#include <string>

#include "union_find.h"
using namespace Rcpp;

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
// Ditto with atoi!
int fast_atoi( const char * str )
{
  int val = 0;
  while( *str ) {
    val = val*10 + (*str++ - '0');
  }
  return val;
}

int which_int(IntegerVector x, int target) {
  int size = (int) x.size();
  for (int i = 0; i < size; ++i) {
    if (x(i) == target) return(i);
  }
  return(-1);
}


// [[Rcpp::export]]
List reach_to_dendrogram(const Rcpp::List reachability, const NumericVector pl_order) {

  // Set up sorted reachability distance
  NumericVector pl = Rcpp::clone(as<NumericVector>(reachability["reachdist"])).sort();

  // Get 0-based order
  IntegerVector order = Rcpp::clone(as<IntegerVector>(reachability["order"])) - 1;

  /// Initialize disjoint-set structure
  int n_nodes = order.size();
  UnionFind uf((size_t) n_nodes);

  // Create leaves
  List dendrogram(n_nodes);
  for (int i = 0; i < n_nodes; ++i) {
    IntegerVector leaf = IntegerVector();
    leaf.push_back(i+1);
    leaf.attr("label") = patch::to_string(i + 1);
    leaf.attr("members") = 1;
    leaf.attr("height") = 0;
    leaf.attr("leaf") = true;
    dendrogram.at(i) = leaf;
  }

  // Precompute the q order
  IntegerVector q_order(n_nodes);
  for (int i = 0; i < n_nodes - 1; ++i) {
    q_order.at(i) = order(which_int(order, pl_order(i)) - 1);
  }

  // Get the index of the point with next smallest reach dist and its neighbor
  IntegerVector members(n_nodes, 1);
  int insert = 0, p = 0, q = 0, p_i = 0, q_i = 0;
  for (int i = 0; i < (n_nodes-1); ++i) {
    p = pl_order(i);
    q = q_order(i);  // left neighbor in ordering
    if (q == -1) { stop("Left neighbor not found"); }

    // Get the actual index of the branch(es) containing the p and q
    p_i = uf.Find(p), q_i = uf.Find(q);
    List branch = List::create(dendrogram.at(q_i), dendrogram.at(p_i));

    // generic proxy blocks attr access for mixed types, so keep track of members manually!
    branch.attr("members") = members.at(p_i) + members.at(q_i);
    branch.attr("height") = pl(i);
    branch.attr("class") = "dendrogram";

    // Merge the two, retrieving the new index
    uf.Union(p_i, q_i);
    insert = uf.Find(q_i); // q because q_branch is first in the new branch

    // Update members reference and insert the branch
    members.at(insert) = branch.attr("members");
    dendrogram.at(insert) = branch;
  }
  return(dendrogram.at(insert));
}

int DFS(List d, List& rp, int pnode, NumericVector stack) {
  if (d.hasAttribute("leaf")) { // If at a leaf node, compare to previous node
    std::string leaf_label = as<std::string>( d.attr("label") );
    rp[leaf_label] = stack; // Record the ancestors reachability values
    std::string pnode_label = patch::to_string(pnode);
    double new_reach = 0.0f;
    if(!rp.containsElementNamed(pnode_label.c_str())) { // 1st time seeing this point
      new_reach = INFINITY;
    } else { // Smallest Common Ancestor
      NumericVector reachdist_p = rp[pnode_label];
      new_reach = min(intersect(stack, reachdist_p));
    }
    NumericVector reachdist = rp["reachdist"];
    IntegerVector order = rp["order"];
    reachdist.push_back(new_reach);
    int res = fast_atoi(leaf_label.c_str());
    order.push_back(res);
    rp["order"] = order;
    rp["reachdist"] = reachdist;
    return(res);
  } else {
    double cheight = d.attr("height");
    stack.push_back(cheight);
    List left = d[0];
    // Recursively go left, recording the reachability distances on the stack
    pnode = DFS(left, rp, pnode, stack);
    if (d.length() > 1) {
      for (int sub_branch = 1; sub_branch < d.length(); ++sub_branch)  {
        pnode = DFS(d[sub_branch], rp, pnode, stack); // pnode;
      }
    }
    return(pnode);
  }
}

// [[Rcpp::export]]
List dendrogram_to_reach(const Rcpp::List x) {
  Rcpp::List rp = List::create(_["order"] = IntegerVector::create(),
                               _["reachdist"] = NumericVector::create());
  NumericVector stack = NumericVector::create();
  DFS(x, rp, 0, stack);
  List res = List::create(_["reachdist"] = rp["reachdist"], _["order"] = rp["order"]);
  res.attr("class") = "reachability";
  return(res);
}

// [[Rcpp::export]]
List mst_to_dendrogram(const NumericMatrix mst) {

  // Set up sorted vector values
  NumericVector p_order = mst(_, 0);
  NumericVector q_order = mst(_, 1);
  NumericVector dist = mst(_, 2);
  int n_nodes = p_order.length() + 1;

  // Make sure to clone so as to not make changes by reference
  p_order = Rcpp::clone(p_order);
  q_order = Rcpp::clone(q_order);

  // UnionFind data structure for fast agglomerative building
  UnionFind uf((size_t) n_nodes);

  // Create leaves
  List dendrogram(n_nodes);
  for (int i = 0; i < n_nodes; ++i) {
    IntegerVector leaf = IntegerVector();
    leaf.push_back(i+1);
    leaf.attr("label") = patch::to_string(i + 1);
    leaf.attr("members") = 1;
    leaf.attr("height") = 0;
    leaf.attr("leaf") = true;
    dendrogram.at(i) = leaf;
  }

  // Get the index of the point with next smallest reach dist and its neighbor
  IntegerVector members(n_nodes, 1);
  int insert = 0, p = 0, q = 0, p_i = 0, q_i = 0;
  for (int i = 0; i < (n_nodes-1); ++i) {
    p = p_order(i), q = q_order(i);

    // Get the actual index of the branch(es) containing the p and q
    p_i = uf.Find(p), q_i = uf.Find(q);

    // Merge the two, retrieving the new index
    uf.Union(p_i, q_i);
    List branch = List::create(dendrogram.at(q_i), dendrogram.at(p_i));

    insert = uf.Find(q_i); // q because q_branch is first in the new branch

    // Update members in the branch
    int tmp_members = members.at(p_i) + members.at(q_i);

    // Branches with equivalent distances are merged simultaneously
    while((i + 1) < (n_nodes-1) && dist(i + 1) == dist(i)){
      i += 1;
      p = p_order(i), q = q_order(i);
      p_i = uf.Find(p), q_i = uf.Find(q);

      // Merge the branches, update current insert index
      int insert2 = uf.Find(q_i);
      branch.push_back(insert == insert2 ? dendrogram.at(p_i) : dendrogram.at(q_i));
      tmp_members += insert == insert2 ? members.at(p_i) : members.at(q_i);
      uf.Union(p_i, q_i);
      insert = uf.Find(q_i);

    }
    // Generic proxy blocks attr access for mixed types, so need to keep track of members manually!
    branch.attr("height") = dist(i);
    branch.attr("class") = "dendrogram";
    branch.attr("members") = tmp_members;

    // Update members reference and insert the branch
    members.at(insert) = branch.attr("members");
    dendrogram.at(insert) = branch;

  }
  return(dendrogram.at(insert));
}


