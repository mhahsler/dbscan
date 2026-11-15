# Find the k Nearest Neighbors

This function finds all k nearest neighbors in a data matrix (including
distances).

## Usage

``` r
kNN(
  x,
  k,
  query = NULL,
  sort = TRUE,
  search = "kdtree",
  bucketSize = 10,
  splitRule = "suggest",
  shrinkRule = "suggest",
  approx = 0
)

# S3 method for class 'kNN'
sort(x, decreasing = FALSE, ...)

# S3 method for class 'kNN'
adjacencylist(x, ...)

# S3 method for class 'kNN'
print(x, ...)
```

## Arguments

- x:

  a numeric data matrix (or data.frame), a
  [dist](https://rdrr.io/r/stats/dist.html) object or a kNN object.

- k:

  number of neighbors to find.

- query:

  a data matrix with the points to query. If query is not specified, the
  NN for all the points in `x` is returned. If query is specified then
  `x` needs to be a data matrix.

- sort:

  sort the neighbors of each point by distance (ties are broken by id)?
  `sort = FALSE` skips this step and returns the neighbors in the order
  found by the search, which is faster. kNN objects can be sorted later
  using [`sort()`](https://rdrr.io/r/base/sort.html).

- search:

  nearest neighbor search strategy (one of `"kdtree"`, `"linear"`,
  `"dist"` or `"bdtree"`). See the section on the search parameters.

- bucketSize:

  max size of the k-d or bd tree leafs.

- splitRule:

  rule to split the k-d or bd tree. One of `"std"`, `"midpt"`, `"fair"`,
  `"sl_midpt"`, `"sl_fair"` or `"suggest"` (SL stands for sliding). The
  default value `"suggest"` uses ANNs best guess.

- shrinkRule:

  shrinking rule used to build the box-decomposition tree (only used for
  `search = "bdtree"`). One of `"none"` (no shrinking, equivalent to a
  k-d tree), `"simple"`, `"centroid"` or `"suggest"`. The default value
  `"suggest"` uses ANNs best guess.

- approx:

  approximation factor for the tree-based searches. The j-th neighbor
  found is at most `(1 + approx)` times as far away as the true j-th
  nearest neighbor. `0` (default) gives exact results; larger values are
  faster, but some true nearest neighbors may be omitted. Ignored by
  `search = "linear"` and `"dist"`.

- decreasing:

  for [`sort()`](https://rdrr.io/r/base/sort.html), sort the neighbors
  by decreasing distance instead of increasing distance?

- ...:

  further arguments (ignored).

## Value

An object of class `kNN` (subclass of
[NN](https://michael.hahsler.net/dbscan/reference/NN.md)) containing a
list with the following components:

- dist :

  a matrix with distances (one row per point in `x` or `query`, one
  column per neighbor).

- id :

  a matrix with the row numbers in `x` of the neighbors.

- k :

  number `k` used.

- sort :

  logical; are the neighbors in each row sorted by increasing distance
  (with ties broken by id)?

- metric :

  used distance metric. This is `"euclidean"` for data matrices. For a
  [dist](https://rdrr.io/r/stats/dist.html) object, it is the method
  stored in the object (which can be `NULL`).

## Details

To speed up k nearest neighbors, a k-d tree is used for data matrices by
default. Box-decomposition trees, linear search, and precomputing
distance matrices are also available.

### Using a kNN object as `x`

If `x` is a `kNN` object, then the first `k` neighbors are returned
(sorted by distance). All other arguments are ignored. An error is
raised if `x` contains fewer than `k` neighbors.

### Sorting, ties and self-matches

- **Sorting:** With `sort = TRUE` (default), the neighbors of each point
  are ordered by distance and ties are broken by id. With
  `sort = FALSE`, the neighbors are returned in the order found by the
  search (the tree-based searches return them by increasing distance,
  but tied neighbors are not ordered by id), which avoids the cost of
  sorting. A `kNN` object can be sorted later using
  [`sort()`](https://rdrr.io/r/base/sort.html).

- **Ties:** If the kth and the (k+1)th nearest neighbor are tied, then
  the neighbor found first is returned and the other one is ignored.

- **Self-matches:** If no `query` is specified, then self-matches are
  removed. If `query` is specified, then a query point that is also in
  `x` is returned as its own nearest neighbor (at distance 0).

### Input requirements

- `x` can be a numeric matrix, a data.frame with only numeric columns, a
  [dist](https://rdrr.io/r/stats/dist.html) object or a `kNN` object.

- `query` needs to be a numeric matrix or data.frame with the same
  number of columns as `x`, and can only be used if `x` is a data
  matrix.

- `k` needs to be a whole number with `1 <= k < nrow(x)`.

### Non-Euclidean distances

To use non-Euclidean distances, a distance matrix needs to be
precomputed by the user and supplied as `x`. Note that this leads to
restrictions of the number of observations that can be used since the
distance matrix grows with \\O(n^2)\\.

### Details on the search parameters:

- `search` controls what method is used if `x` is a data matrix.
  Available methods are:

  - `"kdtree"` uses a Euclidean k-d tree to speed up NN search. This is
    the default with good overall performance.

  - `"bdtree"` uses a box-decomposition tree. Its spatial boxes can, in
    addition to splitting, also shrink. According to the ANN authors,
    this tree type may be faster than a k-d tree on tightly clustered
    data.

  - `"linear"` performs linear search directly on the data.

  - `"dist"` precomputes Euclidean distances first using
    [`dist()`](https://rdrr.io/r/stats/dist.html). While fast for small
    data sets, the memory requirements for the distance matrix grow with
    \\O(n^2)\\.

  The trees and linear search use the highly optimized ANN library
  (Mount and Arya, 2010) for Euclidean distances.

  If `x` is a distance matrix, then `search` is ignored and `"dist"` is
  used.

- Tree specific parameters

  - `bucketSize`, `splitRule` and `shrinkRule` control how the tree is
    built. `bucketSize` and `splitRule` apply to k-d and bd trees;
    `shrinkRule` applies only to bd trees (`"none"` gives a tree that
    behaves like a k-d tree).

  - `approx` uses the approximate nearest neighbor search implemented in
    ANN. If `approx > 0`, the distance to the j-th returned neighbor is
    at most `(1 + approx)` times the distance to the true j-th nearest
    neighbor, so some true nearest neighbors may be omitted. This can
    lead to a significant speedup. The default `approx = 0` gives exact
    results. `approx` is only used by the tree-based searches. For more
    details see Mount and Arya (2010).

### Handling of `NA`s:

The tree-based searches and `search = "linear"` cannot handle `NA`,
`NaN` or infinite values in the data (in `x` or `query`); an error is
raised. (The ANN library does not check its input and would silently
return wrong results.) With `search = "dist"`, the distances are
calculated by [`dist()`](https://rdrr.io/r/stats/dist.html), which
leaves out coordinates that are `NA` for a pair of points and scales up
the remaining sum. These distances are not exactly Euclidean distances.
The resulting distance matrix (and a
[dist](https://rdrr.io/r/stats/dist.html) object supplied as `x`) cannot
contain `NA`s, which happens, for example, if a point has only `NA`
values. Otherwise, remove or impute missing values before performing the
search.

## References

David M. Mount and Sunil Arya (2010). ANN: A Library for Approximate
Nearest Neighbor Searching, <https://www.cs.umd.edu/~mount/ANN/>.

## See also

Other NN functions:
[`NN`](https://michael.hahsler.net/dbscan/reference/NN.md),
[`comps()`](https://michael.hahsler.net/dbscan/reference/comps.md),
[`frNN()`](https://michael.hahsler.net/dbscan/reference/frNN.md),
[`kNNdist()`](https://michael.hahsler.net/dbscan/reference/kNNdist.md),
[`sNN()`](https://michael.hahsler.net/dbscan/reference/sNN.md)

## Author

Michael Hahsler

## Examples

``` r
data(iris)
x <- iris[, -5]

# Example 1: finding kNN for all points in a data matrix (using a k-d tree)
nn <- kNN(x, k = 5)
nn
#> k-nearest neighbors for 150 objects (k=5).
#> Distance metric: euclidean 
#> 
#> Available fields: dist, id, k, sort, metric

# explore neighborhood of point 10
i <- 10
nn$id[i,]
#>  1  2  3  4  5 
#> 35  2 31 13 26 
plot(x, col = ifelse(seq_len(nrow(iris)) %in% nn$id[i,], "red", "black"))


# visualize the 5 nearest neighbors
plot(nn, x)


# visualize a reduced 2-NN graph
plot(kNN(nn, k = 2), x)


# Example 2: find kNN for query points
q <- x[c(1,100),]
nn <- kNN(x, k = 10, query = q)

# we plot the first two features for iris
plot(nn, x[, 1:2], col = "grey")
points(q, pch = 3, lwd = 2)


# Example 3: find kNN using a precomputed (non-Euclidean) distance matrix
d <- dist(x, method = "manhattan")
nn <- kNN(d, k = 1)
plot(nn, x)
```
