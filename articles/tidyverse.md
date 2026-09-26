# Using dbscan with tidyverse

The **dbscan** package provides
[`tidy()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md),
[`augment()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md),
and
[`glance()`](http://michael.hahsler.net/dbscan/reference/dbscan_tidiers.md)
methods for its clustering algorithms, making them easy to use with
tidyverse, ggplot2, and
[tidymodels](https://www.tidymodels.org/learn/statistics/k-means/).

Load the packages and prepare the numeric variables from the iris data:

``` r

library(dbscan)
library(tidyverse)

x <- iris[, 1:4]
db <- x %>% dbscan(eps = .42, minPts = 5)
```

Get cluster statistics as a tibble:

``` r

tidy(db)
#> # A tibble: 4 × 3
#>   cluster  size noise
#>   <fct>   <int> <lgl>
#> 1 0          29 TRUE 
#> 2 1          48 FALSE
#> 3 2          37 FALSE
#> 4 3          36 FALSE
```

Visualize the clustering with ggplot2, using an x for noise points:

``` r

augment(db, x) %>%
  ggplot(aes(x = Petal.Length, y = Petal.Width)) +
  geom_point(aes(color = .cluster, shape = noise)) +
  scale_shape_manual(values = c(19, 4))
```

![DBSCAN clusters in the iris
data](tidyverse_files/figure-html/plot-1.png)
