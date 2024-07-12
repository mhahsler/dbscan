test_that("kNN", {
  set.seed(665544)
  n <- 1000
  x <- cbind(
    x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

  ## no duplicates first! All distances should be unique
  x <- x[!duplicated(x),]

  rownames(x) <- paste0("Object_", seq_len(nrow(x)))

  k <- 5L
  nn <- kNN(x, k=k, sort = TRUE)

  ## check dimensions
  expect_identical(nn$k, k)
  expect_identical(dim(nn$dist), c(nrow(x), k))
  expect_identical(dim(nn$id), c(nrow(x), k))

  ## check visually
  #plot(x)
  #points(x[nn$id[1,],], col="red", lwd=5)
  #points(x[nn$id[2,],], col="green", lwd=5)

  ## compare with kNN found using distances
  nn_d <- kNN(dist(x), k, sort = TRUE)

  ## check visually
  #plot(x)
  #points(x[nn_d$id[1,],], col="red", lwd=5)
  #points(x[nn_d$id[2,],], col="green", lwd=5)

  ### will agree since we use sorting
  expect_equal(nn, nn_d)

  ## calculate dist internally
  nn_d2 <- kNN(x, k, search = "dist", sort = TRUE)
  expect_equal(nn, nn_d2)

  ## without sorting
  nn2 <- kNN(x, k=k, sort = FALSE)
  expect_equal(t(apply(nn$id, MARGIN = 1, sort)),
    t(apply(nn2$id, MARGIN = 1, sort)))

  ## search options
  nn_linear <- kNN(x, k=k, search = "linear", sort = TRUE)
  expect_equal(nn, nn_linear)

  ## split options
  for(so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
    nn3 <- kNN(x, k=k, splitRule = so, sort = TRUE)
    expect_equal(nn, nn3)
  }

  ## bucket size
  for (bs in c(5, 10, 15, 100)) {
    nn3 <- kNN(x, k=k, bucketSize = bs, sort = TRUE)
    expect_equal(nn, nn3)
  }

  ## the order is not stable with matching distances which means that the
  ## k-NN are not stable. We add 100 copied points to check if self match
  ## filtering and sort works
  x <- rbind(x, x[sample(seq_len(nrow(x)), 100),])
  rownames(x) <- paste0("Object_", seq_len(nrow(x)))

  k <- 5L
  nn <- kNN(x, k=k, sort = TRUE)

  ## compare with manually found NNs
  nn_d <- kNN(x, k=k, search = "dist", sort = TRUE)

  expect_equal(nn$dist, nn_d$dist)
  ## This is expected to fail: because the ids are not stable for matching distances
  ## expect_equal(nn$id, nn_d$id)
  ## FIXME: write some code to check this!


  ## missing values, but distances are fine
  x_na <- x
  x_na[c(1, 3, 5), 1] <- NA
  expect_error(kNN(x_na, k = 3), regexp = "NA")
  res_d1 <- kNN(x_na, k = 3, search = "dist")
  res_d2 <- kNN(dist(x_na), k = 3)
  expect_equal(res_d1, res_d2)

  ## introduce NAs into dist
  x_na[c(1, 3, 5),] <- NA
  expect_error(kNN(x_na, k = 3), regexp = "NA")
  expect_error(kNN(x_na, k = 3, search = "dist"), regexp = "NA")
  expect_error(kNN(dist(x_na), k = 3), regexp = "NA")

  ## inf
  x_inf <- x
  x_inf[c(1, 3, 5), 2] <- Inf
  kNN(x_inf, k = 3)
  kNN(x_inf, k = 3, search = "dist")
  kNN(dist(x_inf), k = 3)


  ## sort and kNN to reduce k
  nn10 <- kNN(x, k = 10)
  #nn10 <- kNN(x, k = 10, sort = FALSE)
  ## knn now returns sorted lists
  #expect_equal(nn10$sort, FALSE)
  expect_error(kNN(nn10, k = 11))
  nn5 <- kNN(nn10, k = 5)
  expect_true(nn5$sort)
  expect_identical(ncol(nn5$id), 5L)
  expect_identical(ncol(nn5$dist), 5L)

  ## test with simple data
  x <- data.frame(x=1:10, row.names = LETTERS[1:10])
  nn <- kNN(x, k = 5)
  expect_identical(unname(nn$id[1, ]), 2:6)
  expect_identical(unname(nn$id[5, ]), c(4L, 6L, 3L, 7L, 2L))
  expect_identical(unname(nn$id[10, ]), 9:5)

  ## test kNN with query
  x <- data.frame(x=1:10, row.names = LETTERS[1:10])
  nn <- kNN(x[1:8, , drop=FALSE], x[9:10, , drop = FALSE], k = 5)
  expect_identical(nrow(nn$id), 2L)
  expect_identical(unname(nn$id[1, ]), 8:4)
  expect_identical(unname(nn$id[2, ]), 8:4)

  expect_error(kNN(dist(x[1:8, , drop=FALSE]), x[9:10, , drop = FALSE], k = 5))
})
