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
  size_t size = x.size();
  for (unsigned int i = 0; i < size; ++i) {
    if (x(i) == target) return(i);
  }
  return(-1);
}



// [[Rcpp::export]]
List reach_to_dendrogram(const Rcpp::List reachability, NumericVector pl, NumericVector pl_order) {
  Rcpp::NumericVector reachdist = reachability["reachdist"];
  Rcpp::IntegerVector order = reachability["order"];
  int n_nodes = order.size();
  UnionFind uf((size_t) n_nodes);

  // Create leaves
  std::vector<List> dendrogram(n_nodes);
  for (int i = 0; i < n_nodes; ++i) {
    List leaf = List();
    leaf.push_back(i+1); 
    std::string label = patch::to_string(i + 1);
    leaf.attr("label") = label;
    leaf.attr("members") = 1;
    leaf.attr("height") = 0;
    leaf.attr("leaf") = true;
    dendrogram.at(i) = leaf;
  }

  // Get the index of the point with next smallest reach dist and its neighbor
  int insert = 0, prev_insert = 0, tmp = 0;
  for (int i = 0; i < (n_nodes-1); ++i) {
    int p = pl_order(i); 
    int q = order(which_int(order, p) - 1); // left neighbor in ordering
    if (q == -1) { stop("Left neighbor not found"); }

    // Get the actual index of the branch(es) containing the p and q
    int p_i = uf.Find(p), q_i = uf.Find(q);
    List p_branch = dendrogram.at(p_i), q_branch = dendrogram.at(q_i);
    int p_members = p_branch.attr("members"), q_members = q_branch.attr("members");
    List branch = List::create(q_branch, p_branch);
    branch.attr("members") = p_members + q_members;
    branch.attr("height") = pl(i);
    branch.attr("class") = "dendrogram";
    // if (i < 5 ) Rcpp::Rcout << "p_i: " << p_i << ", q_i: " << q_i << std::endl;

    // Merge the two, retrieving the new index and deleting the old
    uf.Union(p_i, q_i);
    tmp = uf.Find(q_i); // q because q_branch is first in the new branch
    if (tmp != insert) { prev_insert = insert; } 
    insert = tmp;
    dendrogram.at(insert) = branch;
    dendrogram.at(insert == p_i ? q_i : p_i) = NULL;
  }
  Rcpp::List result = dendrogram.at(insert);
  return(result);
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

  
