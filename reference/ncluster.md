# Number of Clusters, Noise Points, and Observations

Extract the number of clusters or the number of noise points for a
clustering. This function works with any clustering result that contains
a list element named `cluster` with a clustering vector. In addition,
`nobs` (see [`stats::nobs()`](https://rdrr.io/r/stats/nobs.html)) is
also available to retrieve the number of clustered points.

## Usage

``` r
ncluster(object, ...)

nnoise(object, ...)
```

## Arguments

- object:

  a clustering result object containing a `cluster` element.

- ...:

  additional arguments are unused.

## Value

returns the number if clusters or noise points.

## See also

Other clustering functions:
[`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md),
[`extractFOSC()`](http://michael.hahsler.net/dbscan/reference/extractFOSC.md),
[`hdbscan()`](http://michael.hahsler.net/dbscan/reference/hdbscan.md),
[`jpclust()`](http://michael.hahsler.net/dbscan/reference/jpclust.md),
[`optics()`](http://michael.hahsler.net/dbscan/reference/optics.md),
[`sNNclust()`](http://michael.hahsler.net/dbscan/reference/sNNclust.md)

## Examples

``` r
data(iris)
iris <- as.matrix(iris[, 1:4])

res <- dbscan(iris, eps = .7, minPts = 5)
res
#> DBSCAN clustering for 150 objects.
#> Parameters: eps = 0.7, minPts = 5
#> Using euclidean distances and borderpoints = TRUE
#> The clustering contains 2 cluster(s) and 3 noise points.
#> 
#>  0  1  2 
#>  3 50 97 
#> 
#> Available fields: cluster, eps, minPts, metric, borderPoints

ncluster(res)
#> [1] 2
nnoise(res)
#> [1] 3
nobs(res)
#> [1] 150

# the functions also work with kmeans and other clustering algorithms.
cl <- kmeans(iris, centers = 3)
ncluster(cl)
#> [1] 3
nnoise(cl)
#> [1] 0
nobs(res)
#> [1] 150
```
