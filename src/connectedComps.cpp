//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#include <Rcpp.h>
using namespace Rcpp;

// Find connected components in kNN and frNN objects.

// [[Rcpp::export]]
IntegerVector comps_kNN(IntegerMatrix nn, bool mutual) {
  R_xlen_t n = nn.nrow();

  // create label vector
  std::vector<int> label(n);
  std::iota(std::begin(label), std::end(label), 1); // Fill with 1, 2, ..., n.
  //iota is C++11 only
  //int value = 1;
  //std::vector<int>::iterator first = label.begin(), last = label.end();
  //while(first != last) *first++ = value++;

  // create sorted sets so we can use set operations
  std::vector< std::set<int> > nn_set(n);
  IntegerVector r;
  std::vector<int> s;
  for(int i = 0; i < n; ++i) {
    r = na_omit(nn(i,_));
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::set<int>::iterator it;
  R_xlen_t i, j;
  int newlabel, oldlabel;

  for(i = 0; i < n; ++i) {
    // check all neighbors of i
    for (it = nn_set[i].begin(); it != nn_set[i].end(); ++it) {
      j = *it-1; // index in nn starts with 1

      // edge was already checked
      //if(j<i) continue;

      // already in the same cluster
      if(label[i] == label[j]) continue;

      // check if points are in each others nn list (i is already in j)
      if(!mutual || nn_set[j].find(i+1) != nn_set[j].end()) {
        if(label[i] > label[j]) {
          newlabel = label[j]; oldlabel = label[i];
        }else{
          newlabel = label[i]; oldlabel = label[j];
        }

        // relabel
        for(int k = 0; k < n; ++k) {
          if(label[k] == oldlabel) label[k] = newlabel;
        }
      }
    }
  }

  return wrap(label);
}

// [[Rcpp::export]]
IntegerVector comps_frNN(List nn, bool mutual) {
  R_xlen_t n = nn.length();

  // create label vector
  std::vector<int> label(n);
  std::iota(std::begin(label), std::end(label), 1); // Fill with 1, 2, ..., n.
  //iota is C++11 only
  //int value = 1;
  //std::vector<int>::iterator first = label.begin(), last = label.end();
  //while(first != last) *first++ = value++;

  // create sorted sets so we can use set operations
  std::vector< std::set<int> > nn_set(n);
  IntegerVector r;
  std::vector<int> s;
  for(R_xlen_t i = 0; i < n; ++i) {
    r = nn[i];
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::set<int>::iterator it;
  R_xlen_t i, j;
  int newlabel, oldlabel;

  for(i = 0; i < n; ++i) {
    // check all neighbors of i
    for (it = nn_set[i].begin(); it != nn_set[i].end(); ++it) {
      j = *it-1; // index in nn starts with 1

      // edge was already checked
      //if(j<i) continue;

      // already in the same cluster
      if(label[i] == label[j]) continue;

      // check if points are in each others nn list (i is already in j)
      if(!mutual || nn_set[j].find(i+1) != nn_set[j].end()) {
        if(label[i] > label[j]) {
          newlabel = label[j]; oldlabel = label[i];
        }else{
          newlabel = label[i]; oldlabel = label[j];
        }

        // relabel
        for(int k = 0; k < n; ++k) {
          if(label[k] == oldlabel) label[k] = newlabel;
        }
      }
    }
  }

  return wrap(label);
}
