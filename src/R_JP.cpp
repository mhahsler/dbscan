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
  R_xlen_t n = nn.nrow();

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
  for(R_xlen_t i = 0; i < n; ++i) {
    r = nn(i,_);
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::vector<int> z;
  std::set<int>::iterator it;
  R_xlen_t i, j;
  int newlabel, oldlabel;

  for(i = 0; i < n; ++i) {
    // check all neighbors of i
    for (it = nn_set[i].begin(); it != nn_set[i].end(); ++it) {
      j = *it-1; // index in nn starts with 1

      // edge was already checked
      if(j<i) continue;

      // already in the same cluster
      if(label[i] == label[j]) continue;

      // check if points are in each others snn list (i is already in j)
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


// jp == true: use the definition by Jarvis-Patrick: A link is created between a pair of
// points, p and q, if and only if p and q have each other  in their k-nearest neighbor lists.
// jp == false: just count the shared NNs
// [[Rcpp::export]]
IntegerMatrix SNN_sim_int(IntegerMatrix nn, LogicalVector jp) {
  R_xlen_t n = nn.nrow();
  R_xlen_t k = nn.ncol();

  IntegerMatrix snn(n, k);

  // create sorted sets so we can use set operations
  std::vector< std::set<int> > nn_set(n);
  IntegerVector r;
  std::vector<int> s;
  for(R_xlen_t i = 0; i < n; ++i) {
    r = nn(i,_);
    s =  as<std::vector<int> >(r);
    nn_set[i].insert(s.begin(), s.end());
  }

  std::vector<int> z;
  int j;

  for(R_xlen_t i = 0; i < n; ++i) {
    // check all neighbors of i
    for (R_xlen_t j_ind = 0; j_ind < k; ++j_ind) {
      j = nn(i, j_ind)-1;

      bool i_in_j = (nn_set[j].find(i+1) != nn_set[j].end());

      if(is_false(all(jp)) || i_in_j) {
        // calculate link strength as the number of shared points
        z.clear();
        std::set_intersection(nn_set[i].begin(), nn_set[i].end(),
          nn_set[j].begin(), nn_set[j].end(),
          std::back_inserter(z));
        snn(i, j_ind) = z.size();
        // +1 if i is in j
        if(i_in_j) snn(i, j_ind)++;

      }else snn(i, j_ind) = 0;

    }
  }

  return snn;
}
