# Package index

## Clustering Algorithms

Cluster data with DBSCAN, HDBSCAN, OPTICS, shared nearest neighbor, and
Jarvis-Patrick algorithms.

- [`dbscan()`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`is.corepoint()`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  [`predict(`*`<dbscan_fast>`*`)`](https://michael.hahsler.net/dbscan/reference/dbscan.md)
  : Density-based Spatial Clustering of Applications with Noise (DBSCAN)
- [`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`print(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`plot(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`coredist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`mrdist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  [`predict(`*`<hdbscan>`*`)`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  : Hierarchical DBSCAN (HDBSCAN)
- [`optics()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`print(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`plot(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.reachability(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`as.dendrogram(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractDBSCAN()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`extractXi()`](https://michael.hahsler.net/dbscan/reference/optics.md)
  [`predict(`*`<optics>`*`)`](https://michael.hahsler.net/dbscan/reference/optics.md)
  : Ordering Points to Identify the Clustering Structure (OPTICS)
- [`sNNclust()`](https://michael.hahsler.net/dbscan/reference/sNNclust.md)
  : Shared Nearest Neighbor Clustering
- [`jpclust()`](https://michael.hahsler.net/dbscan/reference/jpclust.md)
  : Jarvis-Patrick Clustering
- [`extractFOSC()`](https://michael.hahsler.net/dbscan/reference/extractFOSC.md)
  : Framework for the Optimal Extraction of Clusters from Hierarchies
- [`ncluster()`](https://michael.hahsler.net/dbscan/reference/ncluster.md)
  [`nnoise()`](https://michael.hahsler.net/dbscan/reference/ncluster.md)
  : Number of Clusters, Noise Points, and Observations

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

## Clustering Evaluation

Evaluate density-based clusterings with the Density-Based Clustering
Validation index.

- [`dbcv()`](https://michael.hahsler.net/dbscan/reference/dbcv.md) :
  Density-Based Clustering Validation Index (DBCV)

## Cluster Visualization and Hierarchies

Plot clusters and work with reachability plots and cluster hierarchies.

- [`hullplot()`](https://michael.hahsler.net/dbscan/reference/hullplot.md)
  [`clplot()`](https://michael.hahsler.net/dbscan/reference/hullplot.md)
  : Plot Clusters
- [`print(`*`<reachability>`*`)`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  [`plot(`*`<reachability>`*`)`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  [`as.reachability()`](https://michael.hahsler.net/dbscan/reference/reachability.md)
  : Reachability Distances
- [`as.dendrogram()`](https://michael.hahsler.net/dbscan/reference/dendrogram.md)
  : Coercions to Dendrogram

## Tidiers

Turn clustering objects into tidy tibbles and augment the original data.

- [`tidy()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`augment()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  [`glance()`](https://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
  : Turn an dbscan clustering object into a tidy tibble

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
