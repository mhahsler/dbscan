test_that("frNN", {
  set.seed(665544)
  n <- 1000
  x <- cbind(
    x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

  ## no duplicates first!
  #x <- x[!duplicated(x),]

  rownames(x) <- paste0("Object_", seq_len(nrow(x)))

  eps <- .5
  nn <- frNN(x, eps = eps, sort = TRUE)

  ## check dimensions
  expect_identical(nn$eps, eps)
  expect_length(nn$dist, nrow(x))
  expect_length(nn$id, nrow(x))
  for (id in nn$id) expect_type(id, "integer")

  expect_identical(lengths(nn$dist), lengths(nn$id))

  ## check visually
  #plot(x)
  #points(x[nn$id[[1]],], col="red", lwd=5)
  #points(x[nn$id[[2]],], col="green", lwd=5)
  #points(x[1:2,, drop = FALSE], col="blue", pch="+", cex=2)

  ## compare with manually found NNs
  nn_d <- frNN(dist(x), eps = eps, sort = TRUE)
  expect_equal(nn, nn_d)

  nn_d2 <- frNN(x, eps = eps, sort = TRUE, search = "dist")
  expect_equal(nn, nn_d2)

  ## without sorting
  nn2 <- frNN(x, eps = eps, sort = FALSE)
  expect_identical(lapply(nn$id, sort),
    lapply(nn2$id, sort))

  ## search options
  nn_linear <- frNN(x, eps=eps, search = "linear")
  nn_bd <- frNN(x, eps=eps, search = "bdtree")
  expect_equal(sort(nn_bd)$id, sort(nn_linear)$id)
  expect_equal(nn, nn_linear)

  ## split options
  for (so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
    nn3 <- frNN(x, eps=eps, splitRule = so)
    expect_equal(nn, nn3)
  }

  ## shrink rules (bd-tree)
  for (sr in c("none", "simple", "centroid", "suggest")) {
    nn3 <- frNN(x, eps=eps, search = "bdtree", shrinkRule = sr)
    expect_equal(sort(nn3)$id, sort(nn_linear)$id)
  }
  expect_error(frNN(x, eps=eps, search = "bdtree", shrinkRule = "bogus"))

  ## bucket size
  for (bs in c(5, 10, 15, 100)) {
    nn3 <- frNN(x, eps=eps, bucketSize = bs)
    expect_equal(nn, nn3)
  }


  ## add 100 copied points to check if self match filtering works
  x <- rbind(x, x[sample(seq_len(nrow(x)), 100),])
  rownames(x) <- paste0("Object_", seq_len(nrow(x)))

  eps <- .5
  nn <- frNN(x, eps = eps, sort = TRUE)

  ## compare with manually found NNs
  nn_d <- frNN(x, eps = eps, sort = TRUE, search = "dist")

  expect_equal(nn, nn_d)

  ## sort and frNN to reduce eps
  nn5 <- frNN(x, eps = .5, sort = FALSE)
  expect_false(nn5$sort)

  nn5s <- sort(nn5)
  expect_true(nn5s$sort)
  expect_true(all(vapply(nn5s$dist, function(x) !is.unsorted(x), logical(1L))))

  expect_error(frNN(nn5, eps = 1))
  nn2 <- frNN(nn5, eps = .2)
  expect_true(all(vapply(nn2$dist, function(x) all(x <= 0.2), logical(1L))))


  ## test with simple data
  x <- data.frame(x = 1:10, row.names = LETTERS[1:10], check.names = FALSE)
  nn <- frNN(x, eps = 2)
  expect_identical(nn$id[[1]], 2:3)
  expect_identical(nn$id[[5]], c(4L, 6L, 3L, 7L))
  expect_identical(nn$id[[10]], 9:8)

  ## test kNN with query
  x <- data.frame(x = 1:10, row.names = LETTERS[1:10], check.names = FALSE)
  nn <- frNN(x[1:8, , drop=FALSE], x[9:10, , drop = FALSE], eps = 2)

  expect_length(nn$id, 2L)
  expect_identical(nn$id[[1]], 8:7)
  expect_identical(nn$id[[2]], 8L)

  expect_error(frNN(dist(x[1:8, , drop=FALSE]), x[9:10, , drop = FALSE], eps = 2))
})

test_that("frNN sorting and boundary", {
  x <- as.matrix(expand.grid(1:6, 1:6))

  unsorted <- frNN(x, eps = 1.5, sort = FALSE)
  expect_false(unsorted$sort)
  sorted <- frNN(x, eps = 1.5, sort = TRUE)
  expect_true(sorted$sort)
  expect_identical(sort(unsorted), sorted)

  ## ties are broken by id
  for (i in seq_along(sorted$id))
    expect_identical(order(sorted$dist[[i]], sorted$id[[i]]),
      seq_along(sorted$id[[i]]))

  ## decreasing
  dec <- sort(sorted, decreasing = TRUE)
  expect_false(dec$sort)
  for (i in seq_along(dec$dist))
    expect_false(is.unsorted(rev(dec$dist[[i]])))
  expect_identical(sort(dec), sorted)

  ## points exactly at distance eps are neighbors (all search methods)
  y <- cbind(c(0, 1, 3), 0)
  for (s in c("kdtree", "linear", "bdtree", "dist"))
    expect_identical(frNN(y, eps = 1, search = s)$id[[1]], 2L)

  ## isolated points have empty neighborhoods
  expect_identical(frNN(y, eps = .5)$id[[3]], integer(0))
})
