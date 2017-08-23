//----------------------------------------------------------------------
//                    KD Tree class definition
// File:                    R_kd_tree.h
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


// Finalizer 
void finalize_kdtree(ANNkd_tree* kd_tree);
typedef XPtr<ANNkd_tree,PreserveStorage,finalize_kdtree> XPtrKdTree;

// Build KD tree 
SEXP kd_tree_int(NumericMatrix data, int bucketSize, int splitRule);

// Get kd tree object 
ANNkd_tree* getKdTree(SEXP kdtree_ptr);

// Utility functions
void printKdTree(SEXP kdtree_ptr, bool with_pts = false);
void printKdTreeStats(SEXP kdtree_ptr);


