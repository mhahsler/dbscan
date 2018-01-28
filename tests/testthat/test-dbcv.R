library("dbscan")
library("testthat")

context("DBCV")

## Test that it can handle are singleton clusters 
x <- iris[, 1:4]
distx <- dist(x)
sl <- hclust(distx, method = "single")
res <- dbscan::dbcv(x, cl = cutree(sl, h = sl$height[20]), xdist = distx)


testthat::expect_true(class(res) == "numeric")