# Minimum Spanning Tree Hierarchies

Computes the minimum spanning tree (MST) of a graph with Prim's
algorithm and returns it as a single-linkage hierarchy
([hclust](https://rdrr.io/r/stats/hclust.html) object).

## Usage

``` r
mst_dist(x)

mst_mrdist(x, minPts, coredist = NULL)
```

## Arguments

- x:

  for `mst_dist()` a [dist](https://rdrr.io/r/stats/dist.html) object
  with finite distances for at least two objects. For `mst_mrdist()` a
  data matrix (or an object that can be converted to a numeric matrix)
  with finite values.

- minPts:

  integer; minimum number of points used to calculate the core distances
  (see
  [`coredist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)).

- coredist:

  optional vector with the core distance of each point. If `NULL`, then
  it is calculated with
  [`coredist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md).

## Value

an [hclust](https://rdrr.io/r/stats/hclust.html) object representing the
single-linkage hierarchy of the MST. Note: the rows of the `merge`
matrix keep the order in which the MST was built and are **not**
normalized to the convention documented in
[hclust](https://rdrr.io/r/stats/hclust.html) (singletons before
clusters, singletons in increasing index order, clusters in increasing
step order). The raw order is kept so that cluster ids stay consistent
with
[`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md).
Code that relies on this convention (e.g., other packages reading
`$merge` directly) may fail or give unexpected results. The `method`
component is `"robust single"`. `dist.method` is taken from the
`"method"` attribute of `x` for `mst_dist()` (it is omitted if `x` has
none) and is `"mutual reachability (euclidean)"` for `mst_mrdist()`.
`mst_mrdist()` also stores `minPts` in the result.

## Details

- `mst_dist()` computes the MST of the complete graph defined by a
  [dist](https://rdrr.io/r/stats/dist.html) object.

- `mst_mrdist()` computes the MST of the graph of mutual reachability
  distances (see
  [`mrdist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md))
  directly from a data matrix and `minPts`. The core distances are
  calculated first (see
  [`coredist()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)).
  The mutual reachability distances are computed on the fly, so **no
  complete distance matrix is created**. This saves memory (\\O(n d)\\
  instead of \\O(n^2)\\) and is used by
  [`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)
  for data matrices. The result is the same as
  `mst_dist(mrdist(x, minPts))`. Only Euclidean distances are supported.

MSTs are not unique if there are ties in the distances. Many ties are
very likely for mutual reachability distances constructed as
`mrd(a, b) = max(coredist(a), coredist(b), dist(a, b))`. Ties in this
implementation are broken in order of the rows in the input data which
is consistent with many other implementations.

## See also

Other HDBSCAN functions:
[`extractFOSC()`](https://michael.hahsler.net/dbscan/reference/extractFOSC.md),
[`glosh()`](https://michael.hahsler.net/dbscan/reference/glosh.md),
[`hdbscan()`](https://michael.hahsler.net/dbscan/reference/hdbscan.md)

## Examples

``` r
data("moons")

## MST of a distance matrix
hc <- mst_dist(dist(moons))
hc
#> 
#> Call:
#> mst_dist(x = dist(moons))
#> 
#> Cluster method   : robust single 
#> Distance         : euclidean 
#> Number of objects: 100 
#> 

## MST of the mutual reachability graph without creating the distance matrix.
hc_mr <- mst_mrdist(moons, minPts = 5)

## same as
hc_mr2 <- mst_dist(mrdist(moons, minPts = 5))
all.equal(hc_mr$merge, hc_mr2$merge)
#> [1] TRUE

## extract clusters (this is HDBSCAN)
cl <- extractFOSC(hc_mr)
hullplot(moons, cl$cluster, main = "HDBSCAN")
```
