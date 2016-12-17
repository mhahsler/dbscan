# dbscan - Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms - R package

[![CRAN version](http://www.r-pkg.org/badges/version/dbscan)](https://cran.r-project.org/package=dbscan)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/dbscan)](https://cran.r-project.org/package=dbscan)
[![Travis-CI Build Status](https://travis-ci.org/mhahsler/dbscan.svg?branch=master)](https://travis-ci.org/mhahsler/dbscan)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/dbscan?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/dbscan)

 This R package provides a fast C++ reimplementation of several density-based algorithms of the DBSCAN 
 family for spatial data. 
 Includes the __DBSCAN__ (density-based spatial clustering of applications with noise) and 
 __OPTICS/OPTICSXi__ (ordering points to identify the clustering structure) clustering algorithms and the 
 __LOF__ (local outlier factor) algorithm. The implementations uses the kd-tree data 
 structure (from library ANN) for faster k-nearest neighbor search. 
 An R interface to __fast kNN and fixed-radius NN search__ is also provided.

This implementation is typically faster than the native R implementation in package `fpc`, or the 
implementations in WEKA, ELKI and Python's scikit-learn.

## Installation

* __Stable CRAN version:__ install from within R.
* __Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/dbscan/build/artifacts) or install via `install_github("mhahsler/dbscan")` (requires devtools) 

## Examples
```R
library("dbscan")

## use the numeric variables in the iris dataset
data("iris")
x <- as.matrix(iris[, 1:4])
 
## DBSCAN
db <- dbscan(x, eps = .4, minPts = 4)
db
## visualize results (noise is shown in black)
pairs(x, col = db$cluster + 1L)

## LOF (local outlier factor) 
lof <- lof(x, k = 4)
## larger bubbles in the visualization have a larger LOF
pairs(x, cex = lof)

## OPTICS
opt <- optics(x, eps = 1, minPts = 4)
opt

## extract DBSCAN-like clustering 
opt <- extractDBSCAN(opt, eps_cl = .4)

## create a reachability plot (extracted DBSCAN clusters at eps_cl=.4 are colored)
plot(opt)

## plot the extracted DBSCAN clustering
pairs(x, col = opt$cluster + 1L)

## extract a hierarchical clustering using the Xi method (captures clusters of varying density)
opt <- extractXi(opt, xi = .05)
opt
plot(opt)
```

## License 
The dbscan package is licensed under the [GNU General Public License (GPL) Version 3](http://www.gnu.org/licenses/gpl-3.0.en.html). The __OPTICSXi__ R implementation was directly ported from [ELKI](http://elki.dbs.ifi.lmu.de/) frameworks available Java source code (GNU AGPLv3), with explicit permission granted by the original author, [Erich Schubert](http://www.dbs.ifi.lmu.de/cms/Erich_Schubert).  


## Further Information

* Development version of [dbscan on github](https://github.com/mhahsler/dbscan).
* List of changes from [NEWS.md](https://github.com/mhahsler/dbscan/blob/master/NEWS.md)
* [dbscan reference manual](http://cran.r-project.org/web/packages/dbscan/dbscan.pdf)

_Maintainer:_ [Michael Hahsler](http://michael.hahsler.net)


