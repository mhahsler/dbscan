# dbscan: Density-Based Spatial Clustering of Applications with Noise (DBSCAN) and Related Algorithms

A fast reimplementation of several density-based algorithms of the
DBSCAN family. Includes the clustering algorithms DBSCAN (density-based
spatial clustering of applications with noise) and HDBSCAN (hierarchical
DBSCAN), the ordering algorithm OPTICS (ordering points to identify the
clustering structure), shared nearest neighbor clustering, and the
outlier detection algorithms LOF (local outlier factor) and GLOSH
(global-local outlier score from hierarchies). The implementations use
the kd-tree data structure (from library ANN) for faster k-nearest
neighbor search. An R interface to fast kNN and fixed-radius NN search
is also provided. Hahsler, Piekenbrock and Doran (2019)
[doi:10.18637/jss.v091.i01](https://doi.org/10.18637/jss.v091.i01) .

## Key functions

- Clustering:
  [`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md),
  [`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md),
  [`optics()`](http://michael.hahsler.net/dbscan/reference/optics.md),
  [`jpclust()`](http://michael.hahsler.net/dbscan/reference/jpclust.md),
  [`sNNclust()`](http://michael.hahsler.net/dbscan/reference/sNNclust.md)

- Outliers:
  [`lof()`](http://michael.hahsler.net/dbscan/reference/lof.md),
  [`glosh()`](http://michael.hahsler.net/dbscan/reference/glosh.md),
  [`pointdensity()`](http://michael.hahsler.net/dbscan/reference/pointdensity.md)

- Nearest Neighbors:
  [`kNN()`](http://michael.hahsler.net/dbscan/reference/kNN.md),
  [`frNN()`](http://michael.hahsler.net/dbscan/reference/frNN.md),
  [`sNN()`](http://michael.hahsler.net/dbscan/reference/sNN.md)

## References

Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based
Clustering with R. Journal of Statistical Software, 91(1), 1-30.
[doi:10.18637/jss.v091.i01](https://doi.org/10.18637/jss.v091.i01)

## See also

Useful links:

- <https://github.com/mhahsler/dbscan>

- Report bugs at <https://github.com/mhahsler/dbscan/issues>

## Author

**Maintainer**: Michael Hahsler <mhahsler@lyle.smu.edu>
([ORCID](https://orcid.org/0000-0003-2716-1405)) \[copyright holder\]

Authors:

- Michael Hahsler <mhahsler@lyle.smu.edu>
  ([ORCID](https://orcid.org/0000-0003-2716-1405)) \[copyright holder\]

- Matthew Piekenbrock \[copyright holder\]

Other contributors:

- Sunil Arya \[contributor, copyright holder\]

- David Mount \[contributor, copyright holder\]

- Claudia Malzer \[contributor\]
