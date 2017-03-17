//----------------------------------------------------------------------
//                  Jarvis-Patrick Clustering
//----------------------------------------------------------------------
// Copyright (c) 2017 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector JP_int(IntegerMatrix nn, unsigned int kt) {
  int n = nn.nrow();

  // create label vector
  std::vector<int> label(n);
  //iota is C++11 only
  //std::iota(std::begin(label), std::end(label), 1); // Fill with 1, 2, ..., n.
  int value = 1;
  std::vector<int>::iterator first = label.begin(), last = label.end();
  while(first != last) *first++ = value++;

  // create sorted sets so we can use set operations
  std::vector< std::set<int> > nn_set(nn.nrow());
  IntegerVector r;
  std::vector<int> s;
  for(int i = 0; i < n; ++i) {
    r = nn(i,_);
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::vector<int> z;
  std::set<int>::iterator it;
  int i, j, newlabel, oldlabel;

  for(i = 0; i < n; ++i) {
    // check all neighbors of i
    for (it = nn_set[i].begin(); it != nn_set[i].end(); ++it) {
      j = *it-1;

      // edge was already checked
      if(j<i) continue;

      // already in the same cluster
      if(label[i] == label[j]) continue;

      // check if points are in each others snn list (is is already in j)
      if(nn_set[j].find(i+1) != nn_set[j].end()) {

        // calculate link strength as the number of shared points
        z.clear();
        std::set_intersection(nn_set[i].begin(), nn_set[i].end(),
          nn_set[j].begin(), nn_set[j].end(),
          std::back_inserter(z));

        // this could be done faster with set union
        // +1 since i is in j
        if(z.size()+1 >= kt) {
          // update labels
          if(label[i] > label[j]) {
            newlabel = label[j]; oldlabel = label[i];
          }else{
            newlabel = label[i]; oldlabel = label[j];
          }

          for(int k = 0; k < n; ++k) {
            if(label[k] == oldlabel) label[k] = newlabel;
          }
        }
      }
    }
  }

  return wrap(label);
}

// [[Rcpp::export]]
IntegerMatrix SNN_sim_int(IntegerMatrix nn) {
  int n = nn.nrow();
  int k = nn.ncol();

  IntegerMatrix snn(n, k);

  // create sorted sets so we can use set operations
  std::vector< std::set<int> > nn_set(n);
  IntegerVector r;
  std::vector<int> s;
  for(int i = 0; i < n; ++i) {
    r = nn(i,_);
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::vector<int> z;
  int j;

  for(int i = 0; i < n; ++i) {
    // check all neighbors of i
    for (int j_ind = 0; j_ind < k; ++j_ind) {
      j = nn(i, j_ind)-1;

      // edge was already checked
      //if(j<i) continue;

      // check if points are in each others snn list (i is already in j)
      if(nn_set[j].find(i+1) != nn_set[j].end()) {

        // calculate link strength as the number of shared points
        z.clear();
        std::set_intersection(nn_set[i].begin(), nn_set[i].end(),
          nn_set[j].begin(), nn_set[j].end(),
          std::back_inserter(z));

        // +1 for i being in j
        snn(i, j_ind) = z.size()+1;
      }
    }
  }

  return snn;
}
