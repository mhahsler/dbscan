#ifndef DB_UTIL_H
#define DB_UTIL_H

#include <Rcpp.h>
using namespace Rcpp;

#define INDEX_TF(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1)

// std::to_string is apparently a c++11 only thing that crashes appveyor, so using ostringstream it is!
namespace patch { template < typename T > std::string to_string( const T& n ); }

template <typename T, typename C> bool contains (const T& container, const C& key)
{
  if (std::find(container.begin(), container.end(), key) != container.end()){
    return true; 
  } else {
    return false; 
  }
}

IntegerVector lowerTri(IntegerMatrix m);
IntegerVector which_cpp( NumericVector x, double value);
IntegerVector which_cpp( IntegerVector x, int value);
IntegerVector which_geq( IntegerVector x, int value);
NumericVector combine(const NumericVector& t1, const NumericVector& t2);
IntegerVector combine(const IntegerVector& t1, const IntegerVector& t2);
IntegerVector concat_int (List const& container); 

#endif