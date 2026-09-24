# Coercions to Dendrogram

Provides a new generic function to coerce objects to dendrograms with
[`stats::as.dendrogram()`](https://rdrr.io/r/stats/dendrogram.html) as
the default. Additional methods for
[hclust](https://rdrr.io/r/stats/hclust.html),
[hdbscan](http://michael.hahsler.net/dbscan/reference/hdbscan.md) and
[reachability](http://michael.hahsler.net/dbscan/reference/reachability.md)
objects are provided.

## Usage

``` r
as.dendrogram(object, ...)

# Default S3 method
as.dendrogram(object, ...)

# S3 method for class 'hclust'
as.dendrogram(object, ...)

# S3 method for class 'hdbscan'
as.dendrogram(object, ...)

# S3 method for class 'reachability'
as.dendrogram(object, ...)
```

## Arguments

- object:

  the object

- ...:

  further arguments

## Details

Coercion methods for [hclust](https://rdrr.io/r/stats/hclust.html),
[hdbscan](http://michael.hahsler.net/dbscan/reference/hdbscan.md) and
[reachability](http://michael.hahsler.net/dbscan/reference/reachability.md)
objects to dendrogram are provided.

The coercion from `hclust` is a faster C++ reimplementation of the
coercion in package `stats`. The original implementation can be called
using
[`stats::as.dendrogram()`](https://rdrr.io/r/stats/dendrogram.html).

The coercion from
[hdbscan](http://michael.hahsler.net/dbscan/reference/hdbscan.md) builds
the non-simplified HDBSCAN hierarchy as a dendrogram object.
