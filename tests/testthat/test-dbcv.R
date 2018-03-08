library("dbscan")
library("testthat")

context("DBCV")

## Test that it can handle are singleton clusters 
x <- iris[, 1:4]
distx <- dist(x)
sl <- hclust(distx, method = "single")
res <- dbscan::dbcv(x, cl = cutree(sl, h = sl$height[20]), xdist = distx)


testthat::expect_true(class(res) == "numeric")



x <- cbind(rnorm(5000), rnorm(5000))
distx <- dist(x)
cl <- sample(x = 1:10, size = 1000, replace = TRUE)
res <- dbscan::dbcv(x = x, cl = cl, xdist = distx)
microbenchmark::microbenchmark(dbscan::dbcv(x = x, cl = cl, xdist = distx), times = 5L)
