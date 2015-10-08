# dbscan - Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms - R package

[![CRAN version](http://www.r-pkg.org/badges/version/dbscan)](http://cran.r-project.org/web/packages/dbscan/index.html)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/dbscan)](http://cran.r-project.org/web/packages/dbscan/index.html)
[![Travis-CI Build Status](https://travis-ci.org/mhahsler/dbscan.svg?branch=master)](https://travis-ci.org/mhahsler/dbscan)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/dbscan?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/dbscan)

 This R package provides a fast C++ reimplementation of several density-based algorithms of the DBSCAN 
 family for spatial data. 
 Includes the DBSCAN (density-based spatial clustering of applications with noise) and 
 OPTICS (ordering points to identify the clustering structure) clustering algorithms and the 
 LOF (local outlier factor) algorithm. The implementations uses the kd-tree data 
 structure (from library ANN) for faster k-nearest neighbor search. 
 An R interface to fast kNN and fixed-radius NN search is also provided.

This implementation is typically faster than the native R implementaiton in package `fpc`, or the 
implementations in WEKA, ELKI and Python's scikit-learn.

## Simple Example
```
install.packages("dbscan")
library("dbscan")
data("iris")
iris <- as.matrix(iris[,1:4])
 
## run DBSCAN
res <- dbscan(iris, eps = .4, minPts = 4)
res

## visualize results
pairs(iris, col = res$cluster + 1L)
```
