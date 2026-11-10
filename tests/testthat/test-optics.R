test_that("OPTICS", {
  load(test_path("fixtures", "test_data.rda"))
  load(test_path("fixtures", "elki_optics.rda"))

  x <- test_data

  ### run OPTICS
  eps <- .1
  #eps <- .06
  eps_cl <- .1
  minPts <- 10
  res <- optics(x, eps = eps,  minPts = minPts)

  expect_length(res$order, nrow(x))
  expect_type(res$order, "integer")
  expect_length(res$reachdist, nrow(x))
  expect_length(res$coredist, nrow(x))
  expect_identical(res$eps, eps)
  expect_identical(res$minPts, as.integer(minPts))
  expect_type(res$minPts, "integer")

  ### compare with distance based version!
  res_d <- optics(dist(x), eps = eps,  minPts = minPts)
  expect_equal(res, res_d)

  #plot(res)
  #plot(res_d)

  ### compare with elki's result
  expect_equal(res$order, elki$ID)
  expect_equal(round(res$reachdist[res$order], 3), round(elki$reachability, 3))

  ### compare result with DBSCAN
  ### "clustering created from a cluster-ordered is nearly indistinguishable
  ### from a clustering created by DBSCAN. Only some border objects may
  ### be missed"

  # extract DBSCAN clustering
  res <- extractDBSCAN(res, eps_cl = eps_cl)
  #plot(res)

  # are there any clusters with only border points?
  frnn <- frNN(x, eps_cl)
  good <- vapply(frnn$id, function(x) (length(x) + 1L) >= minPts, logical(1L))
  #plot(x, col = (res$cluster+1L))
  c_good <- res$cluster[good]
  c_notgood <- res$cluster[!good]
  expect_false(setdiff(c_notgood, c_good) != 0L)

  # compare with DBSCAN
  db <- dbscan(x, minPts = minPts, eps = eps)
  #plot(x, col = res$cluster+1L)
  #plot(x, col = db$cluster+1L)

  # match clusters (get rid of border points which might differ)
  pure <- vapply(
    split(db$cluster, res$cluster), function(x) length(unique(x)), integer(1L)
  )

  expect_true(all(pure[names(pure) != "0"] == 1L))

  ## missing values, but distances are fine
  x_na <- x
  x_na[c(1,3,5), 1] <- NA
  expect_error(optics(x_na, eps = .2, minPts = 4), regexp = "NA")
  res_d1 <- optics(x_na, eps = .2, minPts = 4, search = "dist")
  res_d2 <- optics(dist(x_na), eps = .2, minPts = 4)
  expect_equal(res_d1, res_d2)

  ## introduce NAs into dist
  x_na[c(1,3,5), 2] <- NA
  expect_error(optics(x_na, eps = .2, minPts = 4), regexp = "NA")
  expect_error(optics(x_na, eps = .2, minPts = 4, search = "dist"),
    regexp = "NA")
  expect_error(optics(dist(x_na), eps = .2, minPts = 4), regexp = "NA")

  ## Create OPTICS-converted and single-linkage dendrograms
  res <- optics(test_data, eps = Inf, minPts = 2)
  res_dend <- as.dendrogram(res)
  reference <- as.dendrogram(hclust(dist(test_data), method = "single"))

  ## Test dendrogram ordering
  expect_equal(as.integer(unlist(res_dend)), res$order)

  ## Test Single Linkage with minPts=2, eps=Inf for strict equivalence
  ## Note: Reordering needed to correct for isomorphisms
  ref_order <- order.dendrogram(reference)
  reference <- reorder(reference, ref_order, agglo.FUN = mean)
  expect_equal(reference, reorder(res_dend, ref_order, agglo.FUN = mean))

  # Make sure any epsilon that queries the entire neighborhood works,
  # error otherwise
  max_rd <- max(res$reachdist[!is.infinite(res$reachdist)], na.rm = TRUE)
  expect_error(as.dendrogram(optics(test_data, eps = max_rd-1e-7,  minPts = 2)), regexp = "Eps")
  expect_error(as.dendrogram(optics(test_data, eps = max_rd, minPts = nrow(test_data) + 1)), regexp = "'minPts'")

  ## Test symmetric relation between reachability <-> dendrogram structures
  expect_equal(as.reachability(as.dendrogram(res))$reachdist, res$reachdist)
  expect_equal(as.reachability(as.dendrogram(res))$order, res$order)
})

test_that("optics with the estimated eps gives every point a core distance", {
  # The estimated eps equals the largest core distance. Rounding errors in the
  # squared distance comparison used to leave one point without a core distance.
  data("moons")
  x_iris <- unique(as.matrix(iris[, 1:4]))

  for (x in list(as.matrix(moons), x_iris)) {
    for (minPts in c(3, 5, 10)) {
      res <- optics(x, minPts = minPts)
      expect_true(all(is.finite(res$coredist)))
      expect_equal(res$coredist, kNNdist(x, k = minPts - 1))
    }
  }
})

test_that("optics with eps = Inf uses an unlimited neighborhood", {
  x <- unique(as.matrix(iris[, 1:4]))
  minPts <- 5

  # the default is an estimated finite eps
  res_est <- optics(x, minPts = minPts)
  expect_true(is.finite(res_est$eps))
  expect_equal(res_est$eps, max(kNNdist(x, k = minPts - 1)), tolerance = 1e-6)

  # Inf is used as is and gives the same result as a very large eps
  res_inf <- optics(x, eps = Inf, minPts = minPts)
  res_large <- optics(x, eps = 1e6, minPts = minPts)
  expect_identical(res_inf$eps, Inf)
  expect_identical(res_inf$order, res_large$order)
  expect_equal(res_inf$reachdist, res_large$reachdist)
  expect_equal(res_inf$coredist, res_large$coredist)

  # only the first reachability distance is undefined
  expect_equal(sum(is.infinite(res_inf$reachdist)), 1L)

  # the estimated eps can leave more undefined reachability distances
  expect_gt(sum(is.infinite(res_est$reachdist)), 1L)

  # distance matrix and distance search give the same result
  res_dist <- optics(dist(x), eps = Inf, minPts = minPts)
  expect_identical(res_dist$order, res_inf$order)
  expect_equal(res_dist$reachdist, res_inf$reachdist)
  res_search <- optics(x, eps = Inf, minPts = minPts, search = "dist")
  expect_equal(res_search$reachdist, res_inf$reachdist)

  # the complete structure can be converted to a dendrogram
  expect_s3_class(as.dendrogram(res_inf), "dendrogram")
})
