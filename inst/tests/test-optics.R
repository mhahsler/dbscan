library("dbscan")
library("testthat")

context("OPTICS")

set.seed(2)
n <- 400

x <- cbind(
  x = runif(4, 0, 1) + rnorm(n, sd=0.1),
  y = runif(4, 0, 1) + rnorm(n, sd=0.1)
)

### run OPTICS
eps <- 1
minPts <- 10
res <- optics(x, eps = eps,  minPts = minPts)

expect_identical(length(res$order), nrow(x))
expect_identical(length(res$reachdist), nrow(x))
expect_identical(length(res$coredist), nrow(x))
expect_identical(res$eps, eps)
expect_identical(res$minPts, minPts)

### FIXME: more tests
