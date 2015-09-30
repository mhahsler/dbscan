library("dbscan")
library("testthat")

context("DBSCAN")

data("iris")

## Species is a factor
expect_error(dbscan::dbscan(iris))

iris <- as.matrix(iris[,1:4])

res <- dbscan::dbscan(iris, eps = .4, minPts = 4)

expect_identical(length(res$cluster), nrow(iris))

## expected result of table(res$cluster) is:
expect_equivalent(table(res$cluster),
    as.table(c("0" = 25L, "1" = 47L, "2" = 38L, "3" = 36L, "4" = 4L)))

## compare with dbscan from package fpc (only if installed)
if (requireNamespace("fpc", quietly = TRUE)) {
    res2 <- fpc::dbscan(iris, eps = .4, MinPts = 4)

    expect_equivalent(res$cluster, res2$cluster)
}


## compare on example data from fpc
set.seed(665544)
n <- 600
x <- cbind(
	   x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
	   y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
	   )

res <- dbscan::dbscan(x, eps = .2, minPts = 4)

expect_identical(length(res$cluster), nrow(x))

## compare with dbscan from package fpc (only if installed)
if (requireNamespace("fpc", quietly = TRUE)) {
    res2 <- fpc::dbscan(x, eps = .2, MinPts = 4)

    expect_equivalent(res$cluster, res2$cluster)
}

