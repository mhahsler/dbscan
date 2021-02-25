# dbscan - Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms - R package

[![CRAN version](https://www.r-pkg.org/badges/version/dbscan)](https://cran.r-project.org/package=dbscan)
[![Rdoc](https://www.rdocumentation.org/badges/version/dbscan)](https://www.rdocumentation.org/packages/dbscan)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/dbscan)](https://cran.r-project.org/package=dbscan)
[![R build status](https://github.com/mhahsler/dbscan/workflows/R-CMD-check/badge.svg)](https://github.com/mhahsler/dbscan/actions)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/dbscan?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/dbscan)

This R package provides a fast C++ (re)implementation of several density-based algorithms with a focus on the DBSCAN family for clustering spatial data.
The package includes: 
 
__Clustering__

- __DBSCAN:__ Density-based spatial clustering of applications with noise.
- __HDBSCAN:__  Hierarchical DBSCAN with simplified hierarchy extraction.
- __OPTICS/OPTICSXi:__ Ordering points to identify the clustering structure clustering algorithms.
- __FOSC:__ Framework for Optimal Selection of Clusters for unsupervised and semisupervised clustering of hierarchical cluster tree.
- __Jarvis-Patrick clustering__
- __SNN Clustering__: Shared Nearest Neighbor Clustering.

__Outlier Detection__

- __LOF:__ Local outlier factor algorithm. 
- __GLOSH:__ Global-Local Outlier Score from Hierarchies algorithm. 

__Fast Nearest-Neighbor Search (using kd-trees)__

- __kNN search__
- __Fixed-radius NN search__

The implementations use the kd-tree data structure (from library ANN) for faster k-nearest neighbor search, and are typically faster than the native R implementations (e.g., dbscan in package `fpc`), or the 
implementations in [WEKA](https://www.cs.waikato.ac.nz/ml/weka/), [ELKI](https://elki-project.github.io/) and [Python's scikit-learn](https://scikit-learn.org/).

## Installation

__Stable CRAN version:__ install from within R with
```R
install.packages("dbscan")
```
__Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/dbscan/build/artifacts) or install from GitHub (needs devtools).
```R 
library("devtools")
install_github("mhahsler/dbscan")
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

Run HDBSCAN (captures stable clusters)
```R
hdb <- hdbscan(x, minPts = 4)
hdb
```

```
HDBSCAN clustering for 150 objects.
Parameters: minPts = 4
The clustering contains 2 cluster(s) and 0 noise points.

  1   2 
100  50 

Available fields: cluster, minPts, cluster_scores, membership_prob, outlier_scores, hc
```

Visualize the results as a simplified tree 
```R
plot(hdb, show_flat = T)
```

See how well each point corresponds to the clusters found by the model used
```R
  colors <- mapply(function(col, i) adjustcolor(col, alpha.f = hdb$membership_prob[i]), 
                   palette()[hdb$cluster+1], seq_along(hdb$cluster))
  plot(x, col=colors, pch=20)
```

## License 
The dbscan package is licensed under the [GNU General Public License (GPL) Version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). The __OPTICSXi__ R implementation was directly ported from the ELKI framework's Java implementation (GNU AGPLv3), with explicit permission granted by the original author, Erich Schubert.  


## Further Information
* List of changes from [NEWS.md](https://github.com/mhahsler/dbscan/blob/master/NEWS.md)

