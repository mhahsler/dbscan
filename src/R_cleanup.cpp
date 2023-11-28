//----------------------------------------------------------------------
//              R interface to dbscan using the ANN library
//----------------------------------------------------------------------
// Copyright (c) 2015 Michael Hahsler. All Rights Reserved.
//
// This software is provided under the provisions of the
// GNU General Public License (GPL) Version 3
// (see: http://www.gnu.org/licenses/gpl-3.0.en.html)

#include <Rcpp.h>
#include "ANN/ANN.h"

using namespace Rcpp;

// [[Rcpp::export]]
void ANN_cleanup() {
  annClose();
}
