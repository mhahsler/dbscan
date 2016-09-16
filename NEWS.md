# Changes in version 0.9-8.1 (2016-xx-x)

* Improved interface for dbscan() and optics(): ... is not passed on to frNN.
* OPTICS clustering extraction methods are now called extractDBSCAN and 
  extractXi.
* hullplot uses now filled polygons as the default.

# Changes in version 0.9-8 (2016-08-05)

* OPTICS: added a predecessor correction step that is used by 
    the ELKI implementation (Matt Piekenbrock).  
* Added hullplot to plot a scatter plot with added convex cluster hulls.
* Fixed a memory problem in frNN (reported by Yilei He).

# Changes in version 0.9-7 (2016-04-14)

* OPTICSXi is now implemented (thanks to Matt Piekenbrock).
* DBSCAN now also accepts MinPts (with a capital M) to be
    compatible with the fpc version.
* DBSCAN objects are now also of class db scan_fast to avoid clashes with fpc.
* DBSCAN and OPTICS have now predict functions.
* Added test for unhandled NAs.
* Fixed LOF for more than k duplicate points (reported by Samneet Singh).

# Changes in version 0.9-6 (2015-12-14)

* OPTICS: fixed second bug reported by Di Pang
* all methods now also accept dist objects and have a search
    method "dist" which precomputes distances.

# Changes in version 0.9-5 (2015-10-04)

* OPTICS: fixed bug with first observation reported by Di Pang
* OPTICS: clusterings can now be extracted using optics_cut

# Changes in version 0.9-4 (2015-09-17)

* added tests (testthat).
* input data is now checked if it can safely be coerced into a
    numeric matrix (storage.mode double).
* fixed self matches in kNN and frNN (now returns the first NN correctly).

# Changes in version 0.9-3 (2015-9-2)

* Added weights to DBSCAN.

# Changes in version 0.9-2 (2015-08-11)

* Added kNN interface.
* Added frNN (fixed radius NN) interface.
* Added LOF.
* Added OPTICS.
* All algorithms check now for interrupt (CTRL-C/Esc).
* DBSCAN now returns a list instead of a numeric vector.

# Changes in version 0.9-1 (2015-07-21)

* DBSCAN: Improved speed by avoiding repeated sorting of point ids.
* Added linear NN search option.
* Added fast calculation for kNN distance.
* fpc and microbenchmark are now used conditionally in the examples.

# Initial Release version 0.9-0 (2015-07-15)
