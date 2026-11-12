# Package index

## Clustering Algorithms

Cluster data with DBSCAN, OPTICS, HDBSCAN, shared nearest neighbor, and
Jarvis-Patrick algorithms.

### DBSCAN

Density-based spatial clustering of applications with noise.

- [`dbscan()`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`is.corepoint()`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`predict(`*`<dbscan_fast>`*`)`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  : Density-based Spatial Clustering of Applications with Noise (DBSCAN)

### OPTICS

Ordering points to identify the clustering structure.

- [`optics()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`print(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`plot(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.reachability(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.dendrogram(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractDBSCAN()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractXi()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`predict(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  : Ordering Points to Identify the Clustering Structure (OPTICS)
- [`print(`*`<reachability>`*`)`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  [`plot(`*`<reachability>`*`)`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  [`as.reachability()`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  : Reachability Distances

### HDBSCAN

Hierarchical DBSCAN, its minimum spanning tree hierarchy, and the
extraction of clusters from hierarchies.

- [`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`print(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`plot(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`coredist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`mrdist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`predict(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  : Hierarchical DBSCAN (HDBSCAN)
- [`extractFOSC()`](https://michael.hahsler.net/dbscan/reference/extractFOSC.md)
  : Framework for the Optimal Extraction of Clusters from Hierarchies
- [`mst_dist()`](https://michael.hahsler.net/dbscan/reference/mst_dist.md)
  [`mst_mrdist()`](https://michael.hahsler.net/dbscan/reference/mst_dist.md)
  : Minimum Spanning Tree Hierarchies

### Shared Nearest Neighbor Clustering

Clustering based on shared nearest neighbor graphs.

- [`sNNclust()`](https://michael.hahsler.net/dbscan/reference/sNNclust.md)
  : Shared Nearest Neighbor Clustering
- [`jpclust()`](https://michael.hahsler.net/dbscan/reference/jpclust.md)
  : Jarvis-Patrick Clustering

## Nearest Neighbor Search

Find k-nearest, fixed-radius, and shared nearest neighbors and work with
nearest-neighbor graphs.

- [`kNN()`](https://michael.hahsler.net/dbscan/reference/kNN.md)
  [`sort(`*`<kNN>`*`)`](https://michael.hahsler.net/dbscan/reference/kNN.md)
  [`adjacencylist(`*`<kNN>`*`)`](https://michael.hahsler.net/dbscan/reference/kNN.md)
  [`print(`*`<kNN>`*`)`](https://michael.hahsler.net/dbscan/reference/kNN.md)
  : Find the k Nearest Neighbors
- [`frNN()`](https://michael.hahsler.net/dbscan/reference/frNN.md)
  [`sort(`*`<frNN>`*`)`](https://michael.hahsler.net/dbscan/reference/frNN.md)
  [`adjacencylist(`*`<frNN>`*`)`](https://michael.hahsler.net/dbscan/reference/frNN.md)
  [`print(`*`<frNN>`*`)`](https://michael.hahsler.net/dbscan/reference/frNN.md)
  : Find the Fixed Radius Nearest Neighbors
- [`sNN()`](https://michael.hahsler.net/dbscan/reference/sNN.md)
  [`sort(`*`<sNN>`*`)`](https://michael.hahsler.net/dbscan/reference/sNN.md)
  [`print(`*`<sNN>`*`)`](https://michael.hahsler.net/dbscan/reference/sNN.md)
  : Find Shared Nearest Neighbors
- [`kNNdist()`](https://michael.hahsler.net/dbscan/reference/kNNdist.md)
  [`kNNdistplot()`](https://michael.hahsler.net/dbscan/reference/kNNdist.md)
  : Calculate and Plot k-Nearest Neighbor Distances
- [`adjacencylist()`](https://michael.hahsler.net/dbscan/reference/NN.md)
  [`sort(`*`<NN>`*`)`](https://michael.hahsler.net/dbscan/reference/NN.md)
  [`plot(`*`<NN>`*`)`](https://michael.hahsler.net/dbscan/reference/NN.md)
  : NN — Nearest Neighbors Superclass
- [`comps()`](https://michael.hahsler.net/dbscan/reference/comps.md) :
  Find Connected Components in a Nearest-neighbor Graph

## Outlier Detection

Calculate local and hierarchical outlier scores and estimate local point
density.

- [`lof()`](https://michael.hahsler.net/dbscan/reference/lof.md) : Local
  Outlier Factor Score
- [`glosh()`](https://michael.hahsler.net/dbscan/reference/glosh.md) :
  Global-Local Outlier Score from Hierarchies
- [`pointdensity()`](https://michael.hahsler.net/dbscan/reference/pointdensity.md)
  : Calculate Local Density at Each Data Point
- [`kNNdist()`](https://michael.hahsler.net/dbscan/reference/kNNdist.md)
  [`kNNdistplot()`](https://michael.hahsler.net/dbscan/reference/kNNdist.md)
  : Calculate and Plot k-Nearest Neighbor Distances

## Cluster Validation

Summarize, evaluate, and visualize clusterings, including reachability
plots and cluster hierarchies.

- [`ncluster()`](https://michael.hahsler.net/dbscan/reference/ncluster.md)
  [`nnoise()`](https://michael.hahsler.net/dbscan/reference/ncluster.md)
  : Number of Clusters, Noise Points, and Observations
- [`dbcv()`](https://michael.hahsler.net/dbscan/reference/dbcv.md) :
  Density-Based Clustering Validation Index (DBCV)
- [`hullplot()`](https://michael.hahsler.net/dbscan/reference/hullplot.md)
  [`clplot()`](https://michael.hahsler.net/dbscan/reference/hullplot.md)
  : Plot Clusters
- [`as.dendrogram()`](https://michael.hahsler.net/dbscan/reference/dendrogram.md)
  : Coercions to Dendrogram

## Tidyverse Tidiers

Turn clustering objects into tidy tibbles and augment the original data.

- [`tidy()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`augment()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`glance()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  : Turn a dbscan clustering object into a tidy tibble

## Data Sets

Example and benchmark data sets for density-based clustering.

- [`moons`](https://michael.hahsler.net/dbscan/reference/moons.md) :
  Moons Data
- [`DS3`](https://michael.hahsler.net/dbscan/reference/DS3.md) : DS3:
  Spatial data with arbitrary shapes
- [`DBCV_datasets`](https://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_1`](https://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_2`](https://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_3`](https://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_4`](https://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  : DBCV Paper Datasets
