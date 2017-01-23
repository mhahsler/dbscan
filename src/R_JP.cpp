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
  std::iota(std::begin(label), std::end(label), 1); // Fill with 1, 2, ..., n.

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
  int newlabel, oldlabel;

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(label[i] == label[j]) continue;

      if(nn_set[j].find(i+1) != nn_set[j].end()) {

        z.clear();
        std::set_intersection(nn_set[i].begin(), nn_set[i].end(),
          nn_set[j].begin(), nn_set[j].end(),
          std::back_inserter(z));

        // this could be done faster with set unions
        if(z.size() >= kt) {
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
