# dbscan 1.1-10 (2022-01-14)

## New Features
* is.corepoint() for DBSCAN.
* coredist() and mrdist() for HDBSCAN.
* find connected components with comps().

# Changes
* reachability plot now shows all undefined distances as a dashed line.

## Bugfix
* memory leak in mrd calculation fixed.

# dbscan 1.1-9 (2022-01-10)

## Changes
* We use now roxygen2.  

## New Features
* Added predict for hdbscan (as suggested by moredatapls)

# dbscan 1.1-8 (2021-04-26)

## Bugfixes
* LOF: fixed numerical issues with k-nearest neighbor distance on Solaris.

# dbscan 1.1-7 (2021-04-21)

## Bugfixes
* Fixed description of k in knndistplot and added minPts argument.
* Fixed bug for tied distances in lof (reported by sverchkov).

## Changes
* lof: the density parameter was changes to minPts to be consistent with the original paper and dbscan. Note that minPts = k + 1.

# dbscan 1.1-6 (2021-02-24)

## Improvements 
* Improved speed of LOF for large ks (following suggestions by eduardokapp). 
* kNN: results is now not sorted again for kd-tree queries which is much faster (by a factor of 10).
* ANN library: annclose() is now only called once when the package is unloaded. This is in preparation to support persistent kd-trees using external pointers.
* hdbscan lost parameter xdist.

## Bugfixes
* removed dependence on methods.
* fixed problem in hullplot for singleton clusters (reported by Fernando Archuby).
* GLOSH now also accepts data.frames.
* GLOSH returns now 0 instead of NaN if we have k duplicate points in the data.

# dbscan 1.1-5 (2019-10-22)

## New Features
* kNN and frNN gained parameter query to query neighbors for points not in the data.
* sNN gained parameter jp to decide if the shared NN should be counted using the definition by Jarvis and Patrick.


# dbscan 1.1-4 (2019-08-05)

## New Features
* kNNdist gained parameter all to indicate if a matrix with the distance to all 
  nearest neighbors up to k should be returned.

## Bugfixes
* kNNdist now correctly returns the distances to the kth neighbor 
  (reported by zschuster).
* dbscan: check eps and minPts parameters to avoid undefined results (reported by ArthurPERE).


# dbscan 1.1-3 (2018-11-12)

## Bugfix
* pointdensity was double counting the query point (reported by Marius Hofert).

# dbscan 1.1-2 (2018-05-18)

## New Features
* OPTICS now calculates eps if it is omitted.

## Bugfix
* Example now only uses igraph conditionally since it is unavailable 
  on Solaris (reported by B. Ripley).

# dbscan 1.1-1 (2017-03-19)

## Bugfix

* Fixed problem with constant name on Solaris in ANN code (reported by B. Ripley).

# dbscan 1.1-0 (2017-03-18)

## New Features

* HDBSCAN was added.
* extractFOSC (optimal selection of clusters for HDBSCAN) was added.
* GLOSH outlier score was added.
* hullplot uses now filled polygons as the default.
* hullplot now used PCA if the data has more than 2 dimensions.
* Added NN superclass for kNN and frNN with plot and with adjacencylist().
* Added shared nearest neighbor clustering as sNNclust() and sNN to calculate
  the number of shared nearest neighbors.
* Added pointdensity function.
* Unsorted kNN and frNN can now be sorted using sort().
* kNN and frNN now also accept kNN and frNN objects, respectively. This can 
  be used to create a new kNN (frNN) with a reduced k or eps.
* Datasets added: DS3 and moon.

## Interface Changes

* Improved interface for dbscan() and optics(): ... it now passed on to frNN.
* OPTICS clustering extraction methods are now called extractDBSCAN and 
  extractXi.
* kNN and frNN are now objects with a print function.
* dbscan now also accepts a frNN object as input.
* jpclust and sNNclust now return a list instead of just the 
  cluster assignments.

# dbscan 1.0-0 (2017-02-02)

## New Features

* The package has now a vignette.
* Jarvis-Patrick clustering is now available as jpclust().
* Improved interface for dbscan() and optics(): ... is now passed on to frNN.
* OPTICS clustering extraction methods are now called extractDBSCAN and 
  extractXi.
* hullplot uses now filled polygons as the default.
* hullplot now used PCA if the data has more than 2 dimensions.
* kNN and frNN are now objects with a print function.
* dbscan now also accepts a frNN object as input.


# dbscan 0.9-8 (2016-08-05)

## New Features

* Added hullplot to plot a scatter plot with added convex cluster hulls.
* OPTICS: added a predecessor correction step that is used by 
    the ELKI implementation (Matt Piekenbrock).  

## Bugfixes

* Fixed a memory problem in frNN (reported by Yilei He).

# dbscan 0.9-7 (2016-04-14)

* OPTICSXi is now implemented (thanks to Matt Piekenbrock).
* DBSCAN now also accepts MinPts (with a capital M) to be
    compatible with the fpc version.
* DBSCAN objects are now also of class db scan_fast to avoid clashes with fpc.
* DBSCAN and OPTICS have now predict functions.
* Added test for unhandled NAs.
* Fixed LOF for more than k duplicate points (reported by Samneet Singh).

# dbscan 0.9-6 (2015-12-14)

* OPTICS: fixed second bug reported by Di Pang
* all methods now also accept dist objects and have a search
    method "dist" which precomputes distances.

# dbscan 0.9-5 (2015-10-04)

* OPTICS: fixed bug with first observation reported by Di Pang
* OPTICS: clusterings can now be extracted using optics_cut

# dbscan 0.9-4 (2015-09-17)

* added tests (testthat).
* input data is now checked if it can safely be coerced into a
    numeric matrix (storage.mode double).
* fixed self matches in kNN and frNN (now returns the first NN correctly).

# dbscan 0.9-3 (2015-9-2)

* Added weights to DBSCAN.

# dbscan 0.9-2 (2015-08-11)

* Added kNN interface.
* Added frNN (fixed radius NN) interface.
* Added LOF.
* Added OPTICS.
* All algorithms check now for interrupt (CTRL-C/Esc).
* DBSCAN now returns a list instead of a numeric vector.

# dbscan 0.9-1 (2015-07-21)

* DBSCAN: Improved speed by avoiding repeated sorting of point ids.
* Added linear NN search option.
* Added fast calculation for kNN distance.
* fpc and microbenchmark are now used conditionally in the examples.

# dbscan 0.9-0 (2015-07-15)

* initial release
