//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.h
//----------------------------------------------------------------------
// Copyright (c) 2016 Michael Hahsler, Matt Piekenbrock. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include <Rcpp.h>
using namespace Rcpp;

class UnionFind
{
  Rcpp::IntegerVector parent;
  Rcpp::IntegerVector rank;
  
  public:
  UnionFind(const int size);
  ~UnionFind();
  void Union(const int x, const int y); 
  const int Find(const int x); 
  
}; // class UnionFind