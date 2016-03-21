library("dbscan")
library("testthat")

context("OPTICS-XI")

load(system.file("test_data/test_data.rda", package = "dbscan"))
load(system.file("test_data/elki_optics.rda", package = "dbscan"))
load(system.file("test_data/elki_optics_xi.rda", package = "dbscan"))

### run OPTICS XI with parameters: xi=0.01, eps=1.0, minPts=5
x <- test_data
res <- optics(x, eps = 1.0,  minPts = 5)
res <- opticsXi(res, xi=0.01, minimum=F)

### Check to make sure ELKI results match R 
testthat::expect_equal(elki_optics_xi, res$clusters_xi)

