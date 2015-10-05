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
eps_optics <- 1
#eps <- .06
eps <- .1
minPts <- 10
res <- optics(x, eps = eps_optics,  minPts = minPts)

expect_identical(length(res$order), nrow(x))
expect_identical(length(res$reachdist), nrow(x))
expect_identical(length(res$coredist), nrow(x))
expect_identical(res$eps, eps_optics)
expect_identical(res$minPts, minPts)

### compare result with DBSCAN
### "clustering created from a cluster-ordered is nearly indistinguishable
### from a clustering created by DBSCAN. Only some border objects may
### be missed"

res <- optics_cut(res, eps_cl = eps)
#plot(res)

db <- dbscan(x, minPts = minPts, eps = eps)
#plot(x, col = res$cluster+1L)
#plot(x, col = db$cluster+1L)

# match clusters (get rid of border points which might differ)
pure <- sapply(split(db$cluster, res$cluster),
  FUN = function(x) length(unique(x)))

expect_true(all(pure[names(pure) != "0"] == 1L))
