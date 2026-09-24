# Find the Fixed Radius Nearest Neighbors

This function uses a kd-tree to find the fixed radius nearest neighbors
(including distances) fast.

## Usage

``` r
frNN(
  x,
  eps,
  query = NULL,
  sort = TRUE,
  search = "kdtree",
  bucketSize = 10,
  splitRule = "suggest",
  approx = 0
)

# S3 method for class 'frNN'
sort(x, decreasing = FALSE, ...)

# S3 method for class 'frNN'
adjacencylist(x, ...)

# S3 method for class 'frNN'
print(x, ...)
```

## Arguments

- x:

  a data matrix, a dist object or a frNN object.

- eps:

  neighbors radius.

- query:

  a data matrix with the points to query. If query is not specified, the
  NN for all the points in `x` is returned. If query is specified then
  `x` needs to be a data matrix.

- sort:

  sort the neighbors by distance? This is expensive and can be done
  later using [`sort()`](https://rdrr.io/r/base/sort.html).

- search:

  nearest neighbor search strategy (one of `"kdtree"`, `"linear"` or
  `"dist"`).

- bucketSize:

  max size of the kd-tree leafs.

- splitRule:

  rule to split the kd-tree. One of `"STD"`, `"MIDPT"`, `"FAIR"`,
  `"SL_MIDPT"`, `"SL_FAIR"` or `"SUGGEST"` (SL stands for sliding).
  `"SUGGEST"` uses ANNs best guess.

- approx:

  use approximate nearest neighbors. All NN up to a distance of a factor
  of `1 + approx` eps may be used. Some actual NN may be omitted leading
  to spurious clusters and noise points. However, the algorithm will
  enjoy a significant speedup.

- decreasing:

  sort in decreasing order?

- ...:

  further arguments

## Value

`frNN()` returns an object of class frNN (subclass of
[NN](http://michael.hahsler.net/dbscan/reference/NN.md)) containing a
list with the following components:

- id :

  a list of integer vectors. Each vector contains the ids (row numbers)
  of the fixed radius nearest neighbors.

- dist :

  a list with distances (same structure as `id`).

- eps :

  neighborhood radius `eps` that was used.

- metric :

  used distance metric.

[`adjacencylist()`](http://michael.hahsler.net/dbscan/reference/NN.md)
returns a list with one entry per data point in `x`. Each entry contains
the id of the nearest neighbors.

## Details

If `x` is specified as a data matrix, then Euclidean distances and fast
nearest neighbor lookup using a kd-tree are used.

To create a frNN object from scratch, you need to supply at least the
elements `id` with a list of integer vectors with the nearest neighbor
ids for each point and `eps` (see below).

**Self-matches:** Self-matches are not returned!

## References

David M. Mount and Sunil Arya (2010). ANN: A Library for Approximate
Nearest Neighbor Searching, <http://www.cs.umd.edu/~mount/ANN/>.

## See also

Other NN functions:
[`NN`](http://michael.hahsler.net/dbscan/reference/NN.md),
[`comps()`](http://michael.hahsler.net/dbscan/reference/comps.md),
[`kNN()`](http://michael.hahsler.net/dbscan/reference/kNN.md),
[`kNNdist()`](http://michael.hahsler.net/dbscan/reference/kNNdist.md),
[`sNN()`](http://michael.hahsler.net/dbscan/reference/sNN.md)

## Author

Michael Hahsler

## Examples

``` r
data(iris)
x <- iris[, -5]

# Example 1: Find fixed radius nearest neighbors for each point
nn <- frNN(x, eps = .5)
nn
#> fixed radius nearest neighbors for 150 objects (eps=0.5).
#> Distance metric: euclidean 
#> 
#> Available fields: dist, id, eps, metric, sort

# Number of neighbors
hist(lengths(adjacencylist(nn)),
  xlab = "k", main="Number of Neighbors",
  sub = paste("Neighborhood size eps =", nn$eps))


# Explore neighbors of point i = 10
i <- 10
nn$id[[i]]
#>  [1] 35  2 31 13 26 50 30 46  3  4  8 36 12 48 40 29 27  1  7 18 41
nn$dist[[i]]
#>  [1] 0.1000000 0.1732051 0.1732051 0.1732051 0.2000000 0.2645751 0.2645751
#>  [8] 0.2645751 0.3162278 0.3162278 0.3316625 0.3464102 0.3464102 0.3464102
#> [15] 0.3741657 0.4472136 0.4472136 0.4690416 0.4795832 0.5000000 0.5000000
plot(x, col = ifelse(seq_len(nrow(iris)) %in% nn$id[[i]], "red", "black"))


# get an adjacency list
head(adjacencylist(nn))
#> [[1]]
#>  [1] 18  5 40 28 29 41  8 50 38 22 49 27 20 47 36 12 11 32 37 21 35 44 10 24
#> 
#> [[2]]
#>  [1] 35 46 13 10 26 31 36  3 50  4 30 48  8 40 12 29 27
#> 
#> [[3]]
#>  [1] 48  4  7 13 46 43 30 35  2 36 10 50 31 39 12  8 41  9 38 26 40
#> 
#> [[4]]
#>  [1] 48 30 31  3 13 46 39 43  9 35 10  7  2 12 26 50
#> 
#> [[5]]
#>  [1] 38  1 18 41  8 40 28 20 22 29 47 50 49 27 12 11 36 44  7
#> 
#> [[6]]
#>  [1] 19 11 49 45 20 47 17 22 33 34
#> 

# plot the fixed radius neighbors (and then reduced to a radius of .3)
plot(nn, x)

plot(frNN(nn, eps = .3), x)


## Example 2: find fixed-radius NN for query points
q <- x[c(1,100),]
nn <- frNN(x, eps = .5, query = q)

plot(nn, x, col = "grey")
points(q, pch = 3, lwd = 2)
```
