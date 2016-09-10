library("dbscan")
library("testthat")

context("OPTICS-XI")

load(system.file("test_data/test_data.rda", package = "dbscan"))
load(system.file("test_data/elki_optics.rda", package = "dbscan"))
load(system.file("test_data/elki_optics_xi.rda", package = "dbscan"))

### run OPTICS XI with parameters: xi=0.01, eps=1.0, minPts=5
x <- test_data
res <- dbscan::optics(x, eps = 1.0,  minPts = 5)
res <- dbscan::extractXi(res, xi=0.10, minimum=F)

### Check to make sure ELKI results match R
testthat::expect_equivalent(elki_optics_xi, res$clusters_xi[, c("start", "end")])
