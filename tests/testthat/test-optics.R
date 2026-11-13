test_that("OPTICS", {
  load(test_path("fixtures", "test_data.rda"))
  load(test_path("fixtures", "elki_optics.rda"))

  x <- test_data

  ### run OPTICS
  eps <- .1
  #eps <- .06
  eps_cl <- .1
  minPts <- 10
  # eps is small on purpose (undefined reachability distances give a warning)
  res <- suppressWarnings(optics(x, eps = eps,  minPts = minPts))

  expect_length(res$order, nrow(x))
  expect_type(res$order, "integer")
  expect_length(res$reachdist, nrow(x))
  expect_length(res$coredist, nrow(x))
  expect_identical(res$eps, eps)
  expect_identical(res$minPts, as.integer(minPts))
  expect_type(res$minPts, "integer")

  ### compare with distance based version!
  res_d <- suppressWarnings(optics(dist(x), eps = eps,  minPts = minPts))
  expect_equal(res, res_d)

  #plot(res)
  #plot(res_d)

  ### compare with elki's result
  expect_equal(res$order, elki$ID)
  expect_equal(round(res$reachdist[res$order], 3), round(elki$reachability, 3))
  expect_equal(res$predecessor[res$order], elki$predecessor)

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
  expect_error(as.dendrogram(suppressWarnings(optics(test_data, eps = max_rd-1e-7,  minPts = 2))), regexp = "Eps")
  expect_error(as.dendrogram(suppressWarnings(optics(test_data, eps = max_rd, minPts = nrow(test_data) + 1))), regexp = "'minPts'")

  ## Test symmetric relation between reachability <-> dendrogram structures
  expect_equal(as.reachability(as.dendrogram(res))$reachdist, res$reachdist)
  expect_equal(as.reachability(as.dendrogram(res))$order, res$order)
})

test_that("optics with the default eps gives every point a core distance", {
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

  # the default is Inf and NULL is treated as Inf
  res_def <- optics(x, minPts = minPts)
  expect_identical(res_def$eps, Inf)
  res_null <- optics(x, eps = NULL, minPts = minPts)
  expect_identical(res_null$eps, Inf)

  # Inf is used as is and gives the same result as a very large eps
  res_inf <- optics(x, eps = Inf, minPts = minPts)
  expect_identical(res_def$order, res_inf$order)
  expect_equal(res_def$reachdist, res_inf$reachdist)
  expect_identical(res_null$order, res_inf$order)
  res_large <- optics(x, eps = 1e6, minPts = minPts)
  expect_identical(res_inf$eps, Inf)
  # iris has many tied distances, so the order can depend on floating point
  # details (e.g., FMA on arm64). Compare order-independent results.
  expect_identical(sort(res_inf$order), sort(res_large$order))
  expect_equal(sort(res_inf$reachdist), sort(res_large$reachdist))
  expect_equal(res_inf$coredist, res_large$coredist)

  # only the first reachability distance is undefined
  expect_equal(sum(is.infinite(res_inf$reachdist)), 1L)

  # a finite eps (here the largest core distance) can leave more undefined
  # reachability distances
  res_est <- suppressWarnings(
    optics(x, eps = max(res_inf$coredist) * (1 + 1e-8), minPts = minPts))
  expect_gt(sum(is.infinite(res_est$reachdist)), 1L)

  # distance matrix and distance search give the same result
  res_dist <- optics(dist(x), eps = Inf, minPts = minPts)
  # iris has many tied distances. Distances computed from the data matrix and
  # by dist() can differ in the last bits (e.g., FMA on arm64), which changes
  # tie-breaking and thus the order. We compare order-independent results.
  expect_equal(res_dist$coredist, res_inf$coredist)
  expect_equal(sort(res_dist$reachdist), sort(res_inf$reachdist))
  res_search <- optics(x, eps = Inf, minPts = minPts, search = "dist")
  expect_equal(sort(res_search$reachdist), sort(res_inf$reachdist))

  # the complete structure can be converted to a dendrogram
  expect_s3_class(as.dendrogram(res_inf), "dendrogram")
})

test_that("OPTICS breaks ties by lower id", {
  # points 2 and 3 have the same reachability distance from point 1
  x <- matrix(c(0, -1, 1), ncol = 1)
  expect_equal(optics(x, eps = 2, minPts = 2)$order, 1:3)
  expect_equal(optics(dist(x), eps = 2, minPts = 2)$order, 1:3)

  # the point with the lower id is processed first (not the one that is
  # first found): swapping the ids of the two tied points swaps their order
  x2 <- matrix(c(0, 1, -1), ncol = 1)
  expect_equal(optics(x2, eps = 2, minPts = 2)$order, 1:3)

  # tied reachability distances are reported
  res <- optics(x, eps = 2, minPts = 2)
  expect_equal(res$reachdist[res$order], c(Inf, 1, 1))
})

test_that("OPTICS with eps = Inf (exact Prim-style computation)", {
  set.seed(1234)
  data("moons")

  datasets <- list(
    moons = as.matrix(moons),
    blobs = do.call(rbind, lapply(1:3, function(k) matrix(rnorm(120, k * 3), ncol = 2))),
    # Many exact ties. We use integer coordinates so that the distances are
    # computed exactly. For decimal values, ties are only ties up to rounding
    # errors, which can be broken differently by different distance
    # implementations.
    ties = round(matrix(runif(300), ncol = 3) * 10),
    dups = rbind(matrix(0, 10, 2), matrix(1, 10, 2), c(5, 5)),
    oned = matrix(c(0, -1, 1, 5, 6, 7.5, 7.5), ncol = 1)
  )

  for (nm in names(datasets)) {
    x <- datasets[[nm]]
    for (minPts in c(2, 3, 5, 10)) {
      if (minPts > nrow(x)) next

      # eps = Inf uses the exact algorithm; a huge finite eps uses the general
      # neighborhood-based implementation.
      res <- optics(x, eps = Inf, minPts = minPts)
      ref <- optics(x, eps = 1e10, minPts = minPts)

      info <- paste(nm, "minPts =", minPts)
      expect_equal(res$order, ref$order, info = info)
      expect_equal(res$reachdist, ref$reachdist, info = info)
      expect_equal(res$coredist, ref$coredist, info = info)
      expect_equal(res$predecessor, ref$predecessor, info = info)

      # dist input gives the same result
      res_d <- optics(dist(x), eps = Inf, minPts = minPts)
      expect_equal(res_d$order, res$order, info = info)
      expect_equal(res_d$reachdist, res$reachdist, info = info)
      expect_equal(res_d$predecessor, res$predecessor, info = info)
    }
  }

  x <- datasets$moons
  res <- optics(x, eps = Inf, minPts = 5)
  expect_s3_class(res, "optics")
  expect_identical(res$eps, Inf)
  expect_identical(res$minPts, 5L)
  expect_true(is.na(res$eps_cl))
  expect_equal(sort(res$order), seq_len(nrow(x)))
  # only the first point has an undefined reachability distance
  expect_equal(sum(is.infinite(res$reachdist)), 1L)
  expect_true(is.infinite(res$reachdist[res$order[1]]))
  expect_true(is.na(res$predecessor[res$order[1]]))
  expect_false(anyNA(res$predecessor[res$order[-1]]))

  # data frames and other search methods
  expect_equal(optics(moons, eps = Inf, minPts = 5)$order, res$order)
  expect_equal(optics(x, eps = Inf, minPts = 5, search = "linear")$order, res$order)

  # the result can be used to extract clusters and to build a dendrogram
  expect_no_error(extractDBSCAN(res, eps_cl = 0.5))
  expect_no_error(extractXi(res, xi = 0.05))
  expect_no_error(as.dendrogram(res))

  # fewer points than minPts falls back to the general implementation
  expect_equal(optics(x[1:3, ], eps = Inf, minPts = 5)$coredist, rep(Inf, 3))

  # approximate search uses the general implementation
  expect_equal(sort(optics(x, eps = Inf, minPts = 5, approx = 0.1)$order),
               seq_len(nrow(x)))

  # invalid input
  x_na <- x
  x_na[1, 1] <- NA
  expect_error(optics(x_na, eps = Inf, minPts = 5))
})

test_that("optics warns about undefined reachability distances for a fixed eps", {
  data("moons")
  x <- as.matrix(moons)

  # eps is too small: points without a reachable core point
  expect_warning(res <- optics(x, eps = 0.3, minPts = 5),
                 "undefined reachability distance")
  expect_gt(sum(is.infinite(res$reachdist)), 1L)
  expect_warning(optics(dist(x), eps = 0.3, minPts = 5),
                 "undefined reachability distance")

  # the message reports the number of additional undefined values
  n_extra <- sum(is.infinite(res$reachdist)) - 1L
  expect_warning(optics(x, eps = 0.3, minPts = 5),
                 paste0("optics: ", n_extra, " of ", nrow(x), " points"))

  # no warning if only the first reachability distance is undefined
  expect_no_warning(optics(x, eps = Inf, minPts = 5))
  expect_no_warning(optics(x, minPts = 5))
  expect_no_warning(optics(x, eps = NULL, minPts = 5))
  expect_no_warning(optics(x, eps = 10, minPts = 5))

  # no warning for an eps that is large enough
  res_inf <- optics(x, eps = Inf, minPts = 5)
  expect_no_warning(res_ok <- optics(x, eps = max(res_inf$reachdist[is.finite(res_inf$reachdist)]) * 1.01,
                                     minPts = 5))
  expect_equal(sum(is.infinite(res_ok$reachdist)), 1L)
})
