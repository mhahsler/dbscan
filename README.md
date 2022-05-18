R package dbscan - Density Based Clustering of Applications with Noise
(DBSCAN) and Related Algorithms
================

[![CRAN
version](https://www.r-pkg.org/badges/version/dbscan)](https://cran.r-project.org/package=dbscan)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/dbscan)](https://cran.r-project.org/package=dbscan)
[![R build
status](https://github.com/mhahsler/dbscan/workflows/R-CMD-check/badge.svg)](https://github.com/mhahsler/dbscan/actions)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/dbscan?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/dbscan)

This R package provides a fast C++ (re)implementation of several
density-based algorithms with a focus on the DBSCAN family for
clustering spatial data. The package includes:

**Clustering**

-   **DBSCAN:** Density-based spatial clustering of applications with
    noise.
-   **HDBSCAN:** Hierarchical DBSCAN with simplified hierarchy
    extraction.
-   **OPTICS/OPTICSXi:** Ordering points to identify the clustering
    structure clustering algorithms.
-   **FOSC:** Framework for Optimal Selection of Clusters for
    unsupervised and semisupervised clustering of hierarchical cluster
    tree.
-   **Jarvis-Patrick clustering**
-   **SNN Clustering**: Shared Nearest Neighbor Clustering.

**Outlier Detection**

-   **LOF:** Local outlier factor algorithm.
-   **GLOSH:** Global-Local Outlier Score from Hierarchies algorithm.

**Fast Nearest-Neighbor Search (using kd-trees)**

-   **kNN search**
-   **Fixed-radius NN search**

The implementations use the kd-tree data structure (from library ANN)
for faster k-nearest neighbor search, and are typically faster than the
native R implementations (e.g., dbscan in package `fpc`), or the
implementations in [WEKA](https://www.cs.waikato.ac.nz/ml/weka/),
[ELKI](https://elki-project.github.io/) and [Python’s
scikit-learn](https://scikit-learn.org/).

## Installation

**Stable CRAN version:** install from within R with

``` r
install.packages("dbscan")
```

**Current development version:** Install from GitHub (needs devtools and
a installed R development environment).

``` r
devtools::install_github("mhahsler/dbscan")
```

## Usage

Load the package and use the numeric variables in the iris dataset

``` r
library("dbscan")

data("iris")
x <- as.matrix(iris[, 1:4])
```

DBSCAN

``` r
db <- dbscan(x, eps = 0.4, minPts = 4)
db
```

    ## DBSCAN clustering for 150 objects.
    ## Parameters: eps = 0.4, minPts = 4
    ## The clustering contains 4 cluster(s) and 25 noise points.
    ## 
    ##  0  1  2  3  4 
    ## 25 47 38 36  4 
    ## 
    ## Available fields: cluster, eps, minPts

Visualize the resulting clustering (noise points are shown in black).

``` r
pairs(x, col = db$cluster + 1L)
```

![](inst/README_files/dbscan-1.png)<!-- -->

OPTICS

``` r
opt <- optics(x, eps = 1, minPts = 4)
opt
```

    ## OPTICS ordering/clustering for 150 objects.
    ## Parameters: minPts = 4, eps = 1, eps_cl = NA, xi = NA
    ## Available fields: order, reachdist, coredist, predecessor, minPts, eps,
    ##                   eps_cl, xi

Extract DBSCAN-like clustering from OPTICS and create a reachability
plot (extracted DBSCAN clusters at eps_cl=.4 are colored)

``` r
opt <- extractDBSCAN(opt, eps_cl = 0.4)
plot(opt)
```

![](inst/README_files/OPTICS_extractDBSCAN-1.png)<!-- -->

HDBSCAN

``` r
hdb <- hdbscan(x, minPts = 4)
hdb
```

    ## HDBSCAN clustering for 150 objects.
    ## Parameters: minPts = 4
    ## The clustering contains 2 cluster(s) and 0 noise points.
    ## 
    ##   1   2 
    ## 100  50 
    ## 
    ## Available fields: cluster, minPts, coredist, cluster_scores,
    ##                   membership_prob, outlier_scores, hc

Visualize the hierarchical clustering as a simplified tree. HDBSCAN
finds 2 stable clusters.

``` r
plot(hdb, show_flat = TRUE)
```

![](inst/README_files/hdbscan-1.png)<!-- -->

## License

The dbscan package is licensed under the [GNU General Public License
(GPL) Version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). The
**OPTICSXi** R implementation was directly ported from the ELKI
framework’s Java implementation (GNU AGPLv3), with explicit permission
granted by the original author, Erich Schubert.

## Changes

-   List of changes from
    [NEWS.md](https://github.com/mhahsler/dbscan/blob/master/NEWS.md)

## References

-   Hahsler M, Piekenbrock M, Doran D (2019). [dbscan: Fast
    Density-Based Clustering with
    R](https://doi.org/10.18637/jss.v091.i01). *Journal of Statistical
    Software*, *91*(1), 1-30. doi: 10.18637/jss.v091.i01.
-   Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu (1996). A
    Density-Based Algorithm for Discovering Clusters in Large Spatial
    Databases with Noise. Institute for Computer Science, University of
    Munich. Proceedings of 2nd International Conference on Knowledge
    Discovery and Data Mining (KDD-96), 226-231.
    <https://dl.acm.org/doi/10.5555/3001460.3001507>
-   Breunig, M., Kriegel, H., Ng, R., and Sander, J. (2000). LOF:
    identifying density-based local outliers. In ACM Int. Conf. on
    Management of Data, pages 93-104. doi:
    <https://doi.org/10.1145/335191.335388>
-   Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander
    (1999). OPTICS: Ordering Points To Identify the Clustering
    Structure. ACM SIGMOD international conference on Management of
    data. ACM Press. pp. doi: <https://doi.org/10.1145/304181.304187>
-   Campello RJGB, Moulavi D, Zimek A, Sander J (2015). Hierarchical
    density estimates for data clustering, visualization, and outlier
    detection. ACM Transactions on Knowledge Discovery from Data (TKDD),
    10(5):1-51. doi: <https://doi.org/10.1145/2733381>
-   R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a Similarity
    Measure Based on Shared Near Neighbors. IEEE Trans. Comput. 22, 11
    (November 1973), 1025-1034. doi:
    <https://doi.org/10.1109/T-C.1973.223640>
-   Levent Ertoz, Michael Steinbach, Vipin Kumar, Finding Clusters of
    Different Sizes, Shapes, and Densities in Noisy, High Dimensional
    Data, SIAM International Conference on Data Mining, 2003, 47-59.
    doi: <https://doi.org/10.1137/1.9781611972733.5>
