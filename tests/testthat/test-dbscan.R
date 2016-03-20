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

## compare with dist-based versions
res_d <- dbscan::dbscan(dist(x), eps = .2, minPts = 4)
expect_identical(res, res_d)
res_d2 <- dbscan::dbscan(x, eps = .2, minPts = 4, search = "dist")
expect_identical(res, res_d2)

## compare with dbscan from package fpc (only if installed)
if (requireNamespace("fpc", quietly = TRUE)) {
    res2 <- fpc::dbscan(x, eps = .2, MinPts = 4)

    expect_equivalent(res$cluster, res2$cluster)
}

## missing values, but distances are fine
x_na <- x
x_na[c(1,3,5), 1] <- NA
expect_error(dbscan::dbscan(x_na, eps = .2, minPts = 4), regexp = "NA")
res_d1 <- dbscan::dbscan(x_na, eps = .2, minPts = 4, search = "dist")
res_d2 <- dbscan::dbscan(dist(x_na), eps = .2, minPts = 4)
expect_equal(res_d1, res_d2)

## introduce NAs into dist
x_na[c(1,3,5), 2] <- NA
expect_error(dbscan::dbscan(x_na, eps = .2, minPts = 4), regexp = "NA")
expect_error(dbscan::dbscan(x_na, eps = .2, minPts = 4, search = "dist"),
  regexp = "NA")
expect_error(dbscan::dbscan(dist(x_na), eps = .2, minPts = 4), regexp = "NA")

