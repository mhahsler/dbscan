test_that("OPTICS", {
  load(system.file("test_data/test_data.rda", package = "dbscan"))
  load(system.file("test_data/elki_optics.rda", package = "dbscan"))

  x <- test_data

  ### run OPTICS
  eps <- .1
  #eps <- .06
  eps_cl <- .1
  minPts <- 10
  res <- optics(x, eps = eps,  minPts = minPts)

  expect_length(res$order, nrow(x))
  expect_length(res$reachdist, nrow(x))
  expect_length(res$coredist, nrow(x))
  expect_identical(res$eps, eps)
  expect_identical(res$minPts, minPts)

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
  res <- optics(test_data, eps = Inf,  minPts = 2)
  res_dend <- as.dendrogram(res)
  reference <- as.dendrogram(hclust(dist(test_data), method = "single"))

  ## Test dendrogram ordering
  expect_equal(as.integer(unlist(res_dend)), res$order)

  ## Test Single Linkage with minPts=2, eps=INF for strict equivalence
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
