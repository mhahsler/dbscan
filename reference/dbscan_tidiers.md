# Turn an dbscan clustering object into a tidy tibble

Provides [tidy()](https://generics.r-lib.org/reference/tidy.html),
[augment()](https://generics.r-lib.org/reference/augment.html), and
[glance()](https://generics.r-lib.org/reference/glance.html) verbs for
clusterings created with algorithms in package `dbscan` to work with
[tidymodels](https://www.tidymodels.org/).

## Usage

``` r
tidy(x, ...)

# S3 method for class 'dbscan'
tidy(x, ...)

# S3 method for class 'hdbscan'
tidy(x, ...)

# S3 method for class 'general_clustering'
tidy(x, ...)

augment(x, ...)

# S3 method for class 'dbscan'
augment(x, data = NULL, newdata = NULL, ...)

# S3 method for class 'hdbscan'
augment(x, data = NULL, newdata = NULL, ...)

# S3 method for class 'general_clustering'
augment(x, data = NULL, newdata = NULL, ...)

glance(x, ...)

# S3 method for class 'dbscan'
glance(x, ...)

# S3 method for class 'hdbscan'
glance(x, ...)

# S3 method for class 'general_clustering'
glance(x, ...)
```

## Arguments

- x:

  An `dbscan` object returned from
  [`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md).

- ...:

  further arguments are ignored without a warning.

- data:

  The data used to create the clustering.

- newdata:

  New data to predict cluster labels for.

## See also

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`generics::augment()`](https://generics.r-lib.org/reference/augment.html),
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html),
[`dbscan()`](http://michael.hahsler.net/dbscan/reference/dbscan.md)

## Examples

``` r

data(iris)
x <- scale(iris[, 1:4])

## dbscan
db <- dbscan(x, eps = .9, minPts = 5)
db
#> DBSCAN clustering for 150 objects.
#> Parameters: eps = 0.9, minPts = 5
#> Using euclidean distances and borderpoints = TRUE
#> The clustering contains 2 cluster(s) and 4 noise points.
#> 
#>  0  1  2 
#>  4 49 97 
#> 
#> Available fields: cluster, eps, minPts, metric, borderPoints

# summarize model fit with tidiers
tidy(db)
#> # A tibble: 3 × 3
#>   cluster  size noise
#>   <fct>   <int> <lgl>
#> 1 0           4 TRUE 
#> 2 1          49 FALSE
#> 3 2          97 FALSE
glance(db)
#> # A tibble: 1 × 3
#>    nobs n.clusters nexcluded
#>   <int>      <int>     <int>
#> 1   150          2         4

# augment for this model needs the original data
augment(db, x)
#> # A tibble: 150 × 6
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster noise
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>    <lgl>
#>  1       -0.898      1.02          -1.34       -1.31 1        FALSE
#>  2       -1.14      -0.132         -1.34       -1.31 1        FALSE
#>  3       -1.38       0.327         -1.39       -1.31 1        FALSE
#>  4       -1.50       0.0979        -1.28       -1.31 1        FALSE
#>  5       -1.02       1.25          -1.34       -1.31 1        FALSE
#>  6       -0.535      1.93          -1.17       -1.05 1        FALSE
#>  7       -1.50       0.786         -1.34       -1.18 1        FALSE
#>  8       -1.02       0.786         -1.28       -1.31 1        FALSE
#>  9       -1.74      -0.361         -1.34       -1.31 1        FALSE
#> 10       -1.14       0.0979        -1.28       -1.44 1        FALSE
#> # ℹ 140 more rows

# to augment new data, the original data is also needed
augment(db, x, newdata = x[1:5, ])
#> # A tibble: 5 × 6
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster noise
#>          <dbl>       <dbl>        <dbl>       <dbl> <fct>    <lgl>
#> 1       -0.898      1.02          -1.34       -1.31 1        FALSE
#> 2       -1.14      -0.132         -1.34       -1.31 1        FALSE
#> 3       -1.38       0.327         -1.39       -1.31 1        FALSE
#> 4       -1.50       0.0979        -1.28       -1.31 1        FALSE
#> 5       -1.02       1.25          -1.34       -1.31 1        FALSE

## hdbscan
hdb <- hdbscan(x, minPts = 5)

# summarize model fit with tidiers
tidy(hdb)
#> # A tibble: 3 × 4
#>   cluster  size cluster_score noise
#>   <fct>   <int>         <dbl> <lgl>
#> 1 0           2          NA   TRUE 
#> 2 1          98         118.  FALSE
#> 3 2          50          87.4 FALSE
glance(hdb)
#> # A tibble: 1 × 3
#>    nobs n.clusters nexcluded
#>   <int>      <int>     <int>
#> 1   150          2         2

# augment for this model needs the original data
augment(hdb, x)
#> # A tibble: 150 × 8
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster .coredist
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>        <dbl>
#>  1       -0.898      1.02          -1.34       -1.31 2            0.236
#>  2       -1.14      -0.132         -1.34       -1.31 2            0.236
#>  3       -1.38       0.327         -1.39       -1.31 2            0.310
#>  4       -1.50       0.0979        -1.28       -1.31 2            0.283
#>  5       -1.02       1.25          -1.34       -1.31 2            0.291
#>  6       -0.535      1.93          -1.17       -1.05 2            0.463
#>  7       -1.50       0.786         -1.34       -1.18 2            0.496
#>  8       -1.02       0.786         -1.28       -1.31 2            0.248
#>  9       -1.74      -0.361         -1.34       -1.31 2            0.551
#> 10       -1.14       0.0979        -1.28       -1.44 2            0.270
#> # ℹ 140 more rows
#> # ℹ 2 more variables: .membership_prob <dbl>, .outlier_scores <dbl>

# to augment new data, the original data is also needed
augment(hdb, x, newdata = x[1:5, ])
#> # A tibble: 5 × 8
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster .coredist
#>          <dbl>       <dbl>        <dbl>       <dbl> <fct>        <dbl>
#> 1       -0.898      1.02          -1.34       -1.31 2               NA
#> 2       -1.14      -0.132         -1.34       -1.31 2               NA
#> 3       -1.38       0.327         -1.39       -1.31 2               NA
#> 4       -1.50       0.0979        -1.28       -1.31 2               NA
#> 5       -1.02       1.25          -1.34       -1.31 2               NA
#> # ℹ 2 more variables: .membership_prob <dbl>, .outlier_scores <dbl>

## Jarvis-Patrick clustering
cl <- jpclust(x, k = 20, kt = 15)

# summarize model fit with tidiers
tidy(cl)
#> # A tibble: 7 × 3
#>   cluster  size noise
#>   <fct>   <int> <lgl>
#> 1 0           0 TRUE 
#> 2 1          49 FALSE
#> 3 2           1 FALSE
#> 4 3          13 FALSE
#> 5 4          36 FALSE
#> 6 5           2 FALSE
#> 7 6          49 FALSE
glance(cl)
#> # A tibble: 1 × 3
#>    nobs n.clusters nexcluded
#>   <int>      <int>     <int>
#> 1   150          6         0

# augment for this model needs the original data
augment(cl, x)
#> # A tibble: 150 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>   
#>  1       -0.898      1.02          -1.34       -1.31 1       
#>  2       -1.14      -0.132         -1.34       -1.31 1       
#>  3       -1.38       0.327         -1.39       -1.31 1       
#>  4       -1.50       0.0979        -1.28       -1.31 1       
#>  5       -1.02       1.25          -1.34       -1.31 1       
#>  6       -0.535      1.93          -1.17       -1.05 1       
#>  7       -1.50       0.786         -1.34       -1.18 1       
#>  8       -1.02       0.786         -1.28       -1.31 1       
#>  9       -1.74      -0.361         -1.34       -1.31 1       
#> 10       -1.14       0.0979        -1.28       -1.44 1       
#> # ℹ 140 more rows

## Shared Nearest Neighbor clustering
cl <- sNNclust(x, k = 20, eps = 0.8, minPts = 15)

# summarize model fit with tidiers
tidy(cl)
#> # A tibble: 3 × 3
#>   cluster  size noise
#>   <fct>   <int> <lgl>
#> 1 0          10 TRUE 
#> 2 1          45 FALSE
#> 3 2          95 FALSE
glance(cl)
#> # A tibble: 1 × 3
#>    nobs n.clusters nexcluded
#>   <int>      <int>     <int>
#> 1   150          2        10

# augment for this model needs the original data
augment(cl, x)
#> # A tibble: 150 × 5
#>    Sepal.Length Sepal.Width Petal.Length Petal.Width .cluster
#>           <dbl>       <dbl>        <dbl>       <dbl> <fct>   
#>  1       -0.898      1.02          -1.34       -1.31 1       
#>  2       -1.14      -0.132         -1.34       -1.31 1       
#>  3       -1.38       0.327         -1.39       -1.31 1       
#>  4       -1.50       0.0979        -1.28       -1.31 1       
#>  5       -1.02       1.25          -1.34       -1.31 1       
#>  6       -0.535      1.93          -1.17       -1.05 1       
#>  7       -1.50       0.786         -1.34       -1.18 1       
#>  8       -1.02       0.786         -1.28       -1.31 1       
#>  9       -1.74      -0.361         -1.34       -1.31 1       
#> 10       -1.14       0.0979        -1.28       -1.44 1       
#> # ℹ 140 more rows
```
