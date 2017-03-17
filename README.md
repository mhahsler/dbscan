# dbscan - Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms - R package

[![CRAN version](http://www.r-pkg.org/badges/version/dbscan)](https://cran.r-project.org/package=dbscan)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/dbscan)](https://cran.r-project.org/package=dbscan)
[![Travis-CI Build Status](https://travis-ci.org/mhahsler/dbscan.svg?branch=master)](https://travis-ci.org/mhahsler/dbscan)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/dbscan?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/dbscan)

 This R package provides a fast C++ reimplementation of several density-based algorithms of the DBSCAN 
 family for spatial data. 
 The package includes: 
 
* __DBSCAN:__ Density-based spatial clustering of applications with noise.
* __OPTICS/OPTICSXi:__ Ordering points to identify the clustering structure clustering algorithms.
* __HDBSCAN:__  Hierarchical DBSCAN with simplified hierarchy extraction.
* __LOF:__ Local outlier factor algorithm. 

The implementations uses the kd-tree data 
 structure (from library ANN) for faster k-nearest neighbor search. 
 An R interface to __fast kNN and fixed-radius NN search__ is provided along with __Jarvis-Patrick clustering__ and __Shared Nearest Neighbor Clustering.__

The implementations are typically faster than the native R implementations (e.g., dbscan in package `fpc`), or the 
implementations in [WEKA](http://www.cs.waikato.ac.nz/ml/weka/), [ELKI](https://elki-project.github.io/) and [Python's scikit-learn](http://scikit-learn.org/).

## Installation

__Stable CRAN version:__ install from within R with
```R
install.packages("dbscan")
```
__Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/dbscan/build/artifacts) or install from GitHub (needs devtools).
```R 
install_git("mhahsler/dbscan")
```

## Usage

Load the package and use the numeric variables in the iris dataset
```R
library("dbscan")

data("iris")
x <- as.matrix(iris[, 1:4])
```

Run DBSCAN
```R
db <- dbscan(x, eps = .4, minPts = 4)
db
```

```
DBSCAN clustering for 150 objects.
Parameters: eps = 0.4, minPts = 4
The clustering contains 4 cluster(s) and 25 noise points.

 0  1  2  3  4 
25 47 38 36  4 

Available fields: cluster, eps, minPts
```

Visualize results (noise is shown in black)
```R
pairs(x, col = db$cluster + 1L)
```


Calculate LOF (local outlier factor) 
and visualize (larger bubbles in the visualization have a larger LOF)
```R
lof <- lof(x, k = 4)
pairs(x, cex = lof)
```

Run OPTICS
```R
opt <- optics(x, eps = 1, minPts = 4)
opt
```

```
OPTICS clustering for 150 objects.
Parameters: minPts = 4, eps = 1, eps_cl = NA, xi = NA
Available fields: order, reachdist, coredist, predecessor, minPts, eps, eps_cl, xi
```

Extract DBSCAN-like clustering from OPTICS 
and create a reachability plot (extracted DBSCAN clusters at eps_cl=.4 are colored)
```R
opt <- extractDBSCAN(opt, eps_cl = .4)
plot(opt)
```

Extract a hierarchical clustering using the Xi method (captures clusters of varying density)
```R
opt <- extractXi(opt, xi = .05)
opt
plot(opt)
```

## License 
The dbscan package is licensed under the [GNU General Public License (GPL) Version 3](http://www.gnu.org/licenses/gpl-3.0.en.html). The __OPTICSXi__ R implementation was directly ported from the [ELKI](http://elki.dbs.ifi.lmu.de/) framework's Java implementation (GNU AGPLv3), with explicit permission granted by the original author, [Erich Schubert](http://www.dbs.ifi.lmu.de/cms/Erich_Schubert).  


## Further Information

* Development version of [dbscan on github](https://github.com/mhahsler/dbscan).
* List of changes from [NEWS.md](https://github.com/mhahsler/dbscan/blob/master/NEWS.md)
* [dbscan reference manual](https://cran.r-project.org/package=dbscan/dbscan.pdf)

_Maintainer:_ [Michael Hahsler](http://michael.hahsler.net)


