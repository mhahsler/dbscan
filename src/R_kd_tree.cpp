//----------------------------------------------------------------------
//                     KD Tree class definition
// File:                   R_kd_tree.cpp
//----------------------------------------------------------------------
// Copyright (c) 2017 Michael Hahsler, Matthew Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
#include "ANN/ANN.h"
#include "ANN/ANNperf.h"
#include "R_kNN.h"

void finalize_kdtree(ANNkd_tree* kd_tree){
  //Rcout << "Finalizer called!\n";
  // cleanup
  if(kd_tree != NULL){
    ANNpointArray dataPoints = kd_tree->thePoints();
    annDeallocPts(dataPoints);
    // delete kd_tree; // Don't need to delete tree; it standard finalizer will do that part
    // annClose(); // remove trivial node
  }
}
typedef XPtr<ANNkd_tree,PreserveStorage,finalize_kdtree> XPtrKdTree;

// Return KDtree as SEXP pointer
// [[Rcpp::export]]
SEXP kd_tree_int(NumericMatrix data, int bucketSize, int splitRule) {
  // copy data
  int nrow = data.nrow();
  int ncol = data.ncol();
  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = data(i, j);
    }
  }
  // create kd-tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrow, ncol, bucketSize, (ANNsplitRule) splitRule);
  XPtrKdTree kdtree_xptr((ANNkd_tree*) kdTree); // register finalizer 
  return(kdtree_xptr);
}

ANNkd_tree* getKdTree(SEXP kdtree_ptr){
  ANNkd_tree* kd_tree;
  if( TYPEOF(kdtree_ptr) == EXTPTRSXP){
    Rcpp::XPtr< ANNkd_tree > kdtree_xptr((SEXP) kdtree_ptr);
    kd_tree = (ANNkd_tree*) &(*kdtree_xptr);
    return kd_tree;
  } else {
    Rcpp::stop("Recieved invalid or non-external pointer.");
  }
}

// [[Rcpp::export]]
void printKdTree(SEXP kdtree_ptr, bool with_pts){
  ANNkd_tree* kd_tree = getKdTree(kdtree_ptr);
  if (kd_tree != NULL){
    kd_tree->Print((ANNbool) with_pts);
  }
}

// [[Rcpp::export]]
void printKdTreeStats(SEXP kdtree_ptr){
  ANNkdStats tree_stats = ANNkdStats(); 
  ANNkd_tree* kd_tree = getKdTree(kdtree_ptr);
  if (kd_tree != NULL){
    kd_tree->getStats(tree_stats);
    Rcout << "Tree built from " << tree_stats.n_pts << " points in " << tree_stats.dim << " dimensional space" << std::endl; 
    Rcout << "Split nodes = " << tree_stats.n_spl; 
    Rcout << ", leaf nodes = " << tree_stats.n_lf; 
    Rcout << ", Avg. leaf aspect ratio = " << tree_stats.avg_ar << std::endl;  
  }
}

// [[Rcpp::export]]
List test_KDtree(SEXP kdtree_ptr, int k, double approx = 0){
  
  ANNkd_tree* kd_tree = getKdTree(kdtree_ptr);
  const int n = kd_tree->nPoints(); 
  NumericMatrix d(n, k); // knn distances
  IntegerMatrix id(n, k); // knn point ids 
  
  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];
  
  ANNpointArray dataPts = kd_tree->thePoints();
  for (int i=0; i<n; i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();
    
    ANNpoint queryPt = dataPts[i];
    
    kd_tree->annkSearch(queryPt, k+1, nnIdx, dists, approx);
    // remove self match
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;
    
    NumericVector ndists = NumericVector(dists, dists+k+1)[take];
    d(i, _) = sqrt(ndists);
  }
  
  // cleanup - finalizer takes care of tree and original data points tree was built with
  delete [] dists;
  delete [] nnIdx;
  
  // prepare results
  List ret;
  ret["dist"] = d;
  ret["id"] = id;
  ret["k"] = k;
  return ret;
}




