test_that("OPTICS-XI", {
  load(test_path("fixtures", "test_data.rda"))
  load(test_path("fixtures", "elki_optics.rda"))
  load(test_path("fixtures", "elki_optics_xi.rda"))

  ### run OPTICS XI with parameters: xi=0.01, eps=1.0, minPts=5
  x <- test_data
  res <- optics(x, eps = 1.0,  minPts = 5)
  res <- extractXi(res, xi = 0.10, minimum = FALSE)

  ### Check to make sure ELKI results match R
  expected <- res$clusters_xi[, c("start", "end")]
  class(expected) <- "data.frame"
  expect_identical(elki_optics_xi, expected)
})
