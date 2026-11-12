test_that("OPTICS-XI", {
  load(test_path("fixtures", "test_data.rda"))
  load(test_path("fixtures", "elki_optics_xi.rda"))

  ### The reference results were created with ELKI 0.8.0
  ### (OPTICSXi with OPTICSHeap, eps = 1.0 and minPts = 5)
  x <- test_data
  res <- optics(x, eps = 1.0, minPts = 5)

  ### the cluster order needs to be the same
  expect_equal(res$order, elki$ID)
  expect_equal(res$predecessor[res$order], elki$predecessor)

  xi_clusters <- function(xi, ...) {
    expected <- extractXi(res, xi = xi, minimum = FALSE, ...)$clusters_xi[, c("start", "end")]
    class(expected) <- "data.frame"
    expected
  }

  ### xi = 0.10
  expect_identical(elki_optics_xi, xi_clusters(0.10))

  ### xi = 0.05
  expect_identical(elki_optics_xi_005, xi_clusters(0.05))

  ### without the predecessor correction (ELKI option -opticsxi.nocorrect)
  expect_identical(elki_optics_xi_005_nocorrect,
                   xi_clusters(0.05, correctPredecessors = FALSE))
})

test_that("OPTICS-XI cluster labels", {
  load(test_path("fixtures", "test_data.rda"))
  res <- optics(test_data, eps = 1.0, minPts = 5)

  res_xi <- extractXi(res, xi = 0.10)
  expect_length(res_xi$cluster, nrow(test_data))
  expect_true(all(res_xi$clusters_xi$end - res_xi$clusters_xi$start + 1 >= 5))

  # minimum = TRUE only keeps non-overlapping clusters
  res_min <- extractXi(res, xi = 0.10, minimum = TRUE)
  expect_lte(nrow(res_min$clusters_xi), nrow(res_xi$clusters_xi))
  expect_false(is.unsorted(res_min$clusters_xi$start))

  expect_error(extractXi(res, xi = 1))
  expect_error(extractXi(res, xi = 0))
})
