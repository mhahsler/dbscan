//----------------------------------------------------------------------
//                              Region Query
// File:                        R_regionQuery.h
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#ifndef REGIONQUERY_H
#define REGIONQUERY_H

#include <Rcpp.h>
#include "ANN/ANN.h"
using namespace Rcpp;

// pair of ids and dists
typedef std::pair< std::vector<int>, std::vector<double> > nn ;

// Note: Region query returns self-matches!

// these function takes an id for the points in the k-d tree
nn regionQueryDist(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);

std::vector<int> regionQuery(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);

// these function takes an query point not in the tree
nn regionQueryDist_point(ANNpoint queryPt, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);

std::vector<int> regionQuery_point(ANNpoint queryPt, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);

#endif
