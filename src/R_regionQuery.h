//----------------------------------------------------------------------
//                              Region Query
// File:                        regionQuery.h
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

typedef std::pair< std::vector<int>, std::vector<double> > nn ;

nn regionQueryDist(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);

std::vector<int> regionQuery(int id, ANNpointArray dataPts, ANNpointSet* kdTree,
  double eps2, double approx);



#endif
