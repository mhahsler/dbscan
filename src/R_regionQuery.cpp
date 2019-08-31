//----------------------------------------------------------------------
//                              Region Query
// File:                        regionQuery.cpp
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)



#include <Rcpp.h>
#include "R_regionQuery.h"

using namespace Rcpp;

// Note: Region query returns self-matches!

// these function takes an id for the points in the k-d tree
nn regionQueryDist(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx) {

  // find fixed radius nearest neighbors
  ANNpoint queryPt = dataPts[id];
  std::pair< std::vector<int>, std::vector<double> > ret =
    kdTree->annkFRSearch2(queryPt, eps2, approx);
  // Note: the points are not sorted by distance!

  return(ret);
}

std::vector<int> regionQuery(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx = 0) {

  // find fixed radius nearest neighbors
  ANNpoint queryPt = dataPts[id];
  std::pair< std::vector<int>, std::vector<double> > ret =
    kdTree->annkFRSearch2(queryPt, eps2, approx);
  // Note: the points are not sorted by distance!

  return(ret.first);
}


// these function takes an query point not in the tree
nn regionQueryDist_point(ANNpoint queryPt, ANNpointArray dataPts,
	ANNpointSet* kdTree, double eps2, double approx) {

  // find fixed radius nearest neighbors
  std::pair< std::vector<int>, std::vector<double> > ret =
    kdTree->annkFRSearch2(queryPt, eps2, approx);
  // Note: the points are not sorted by distance!

  return(ret);
}

std::vector<int> regionQuery_point(ANNpoint queryPt, ANNpointArray dataPts,
	ANNpointSet* kdTree, double eps2, double approx = 0) {

  // find fixed radius nearest neighbors
  std::pair< std::vector<int>, std::vector<double> > ret =
    kdTree->annkFRSearch2(queryPt, eps2, approx);
  // Note: the points are not sorted by distance!

  return(ret.first);
}

