# Package index

## Clustering Algorithms

Cluster data with DBSCAN, HDBSCAN, OPTICS, shared nearest neighbor, and
Jarvis-Patrick algorithms.

- [`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`is.corepoint()`](http://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`predict(`*`<dbscan_fast>`*`)`](http://michael.hahsler.net/dbscan/reference/dbscan.md)
  : Density-based Spatial Clustering of Applications with Noise (DBSCAN)
- [`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`print(`*`<hdbscan>`*`)`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`plot(`*`<hdbscan>`*`)`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`coredist()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`mrdist()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`predict(`*`<hdbscan>`*`)`](http://michael.hahsler.net/dbscan/reference/hdbscan.md)
  : Hierarchical DBSCAN (HDBSCAN)
- [`optics()`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`print(`*`<optics>`*`)`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`plot(`*`<optics>`*`)`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.reachability(`*`<optics>`*`)`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.dendrogram(`*`<optics>`*`)`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractDBSCAN()`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractXi()`](http://michael.hahsler.net/dbscan/reference/optics.md)
  [`predict(`*`<optics>`*`)`](http://michael.hahsler.net/dbscan/reference/optics.md)
  : Ordering Points to Identify the Clustering Structure (OPTICS)
- [`sNNclust()`](http://michael.hahsler.net/dbscan/reference/sNNclust.md)
  : Shared Nearest Neighbor Clustering
- [`jpclust()`](http://michael.hahsler.net/dbscan/reference/jpclust.md)
  : Jarvis-Patrick Clustering
- [`extractFOSC()`](http://michael.hahsler.net/dbscan/reference/extractFOSC.md)
  : Framework for the Optimal Extraction of Clusters from Hierarchies
- [`ncluster()`](http://michael.hahsler.net/dbscan/reference/ncluster.md)
  [`nnoise()`](http://michael.hahsler.net/dbscan/reference/ncluster.md)
  : Number of Clusters, Noise Points, and Observations

## Nearest Neighbor Search

Find k-nearest, fixed-radius, and shared nearest neighbors and work with
nearest-neighbor graphs.

- [`kNN()`](http://michael.hahsler.net/dbscan/reference/kNN.md)
  [`sort(`*`<kNN>`*`)`](http://michael.hahsler.net/dbscan/reference/kNN.md)
  [`adjacencylist(`*`<kNN>`*`)`](http://michael.hahsler.net/dbscan/reference/kNN.md)
  [`print(`*`<kNN>`*`)`](http://michael.hahsler.net/dbscan/reference/kNN.md)
  : Find the k Nearest Neighbors
- [`frNN()`](http://michael.hahsler.net/dbscan/reference/frNN.md)
  [`sort(`*`<frNN>`*`)`](http://michael.hahsler.net/dbscan/reference/frNN.md)
  [`adjacencylist(`*`<frNN>`*`)`](http://michael.hahsler.net/dbscan/reference/frNN.md)
  [`print(`*`<frNN>`*`)`](http://michael.hahsler.net/dbscan/reference/frNN.md)
  : Find the Fixed Radius Nearest Neighbors
- [`sNN()`](http://michael.hahsler.net/dbscan/reference/sNN.md)
  [`sort(`*`<sNN>`*`)`](http://michael.hahsler.net/dbscan/reference/sNN.md)
  [`print(`*`<sNN>`*`)`](http://michael.hahsler.net/dbscan/reference/sNN.md)
  : Find Shared Nearest Neighbors
- [`kNNdist()`](http://michael.hahsler.net/dbscan/reference/kNNdist.md)
  [`kNNdistplot()`](http://michael.hahsler.net/dbscan/reference/kNNdist.md)
  : Calculate and Plot k-Nearest Neighbor Distances
- [`adjacencylist()`](http://michael.hahsler.net/dbscan/reference/NN.md)
  [`sort(`*`<NN>`*`)`](http://michael.hahsler.net/dbscan/reference/NN.md)
  [`plot(`*`<NN>`*`)`](http://michael.hahsler.net/dbscan/reference/NN.md)
  : NN — Nearest Neighbors Superclass
- [`comps()`](http://michael.hahsler.net/dbscan/reference/comps.md) :
  Find Connected Components in a Nearest-neighbor Graph

## Outlier Detection

Calculate local and hierarchical outlier scores and estimate local point
density.

- [`lof()`](http://michael.hahsler.net/dbscan/reference/lof.md) : Local
  Outlier Factor Score
- [`glosh()`](http://michael.hahsler.net/dbscan/reference/glosh.md) :
  Global-Local Outlier Score from Hierarchies
- [`pointdensity()`](http://michael.hahsler.net/dbscan/reference/pointdensity.md)
  : Calculate Local Density at Each Data Point
- [`kNNdist()`](http://michael.hahsler.net/dbscan/reference/kNNdist.md)
  [`kNNdistplot()`](http://michael.hahsler.net/dbscan/reference/kNNdist.md)
  : Calculate and Plot k-Nearest Neighbor Distances

## Clustering Evaluation

Evaluate density-based clusterings with the Density-Based Clustering
Validation index.

- [`dbcv()`](http://michael.hahsler.net/dbscan/reference/dbcv.md) :
  Density-Based Clustering Validation Index (DBCV)

## Cluster Visualization and Hierarchies

Plot clusters and work with reachability plots and cluster hierarchies.

- [`hullplot()`](http://michael.hahsler.net/dbscan/reference/hullplot.md)
  [`clplot()`](http://michael.hahsler.net/dbscan/reference/hullplot.md)
  : Plot Clusters
- [`print(`*`<reachability>`*`)`](http://michael.hahsler.net/dbscan/reference/reachability.md)
  [`plot(`*`<reachability>`*`)`](http://michael.hahsler.net/dbscan/reference/reachability.md)
  [`as.reachability()`](http://michael.hahsler.net/dbscan/reference/reachability.md)
  : Reachability Distances
- [`as.dendrogram()`](http://michael.hahsler.net/dbscan/reference/dendrogram.md)
  : Coercions to Dendrogram

## Tidiers

Turn clustering objects into tidy tibbles and augment the original data.

- [`tidy()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`augment()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`glance()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  : Turn an dbscan clustering object into a tidy tibble

## Data Sets

Example and benchmark data sets for density-based clustering.

- [`moons`](http://michael.hahsler.net/dbscan/reference/moons.md) :
  Moons Data
- [`DS3`](http://michael.hahsler.net/dbscan/reference/DS3.md) : DS3:
  Spatial data with arbitrary shapes
- [`DBCV_datasets`](http://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_1`](http://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_2`](http://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_3`](http://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  [`Dataset_4`](http://michael.hahsler.net/dbscan/reference/DBCV_datasets.md)
  : DBCV Paper Datasets
