test_that("sNN", {
  set.seed(665544)
  n <- 1000
  x <- cbind(
    x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

  ## no duplicates first!
  x <- x[!duplicated(x),]

  rownames(x) <- paste0("Object_", seq_len(nrow(x)))

  k <- 5L
  nn <- sNN(x, k=k, sort = TRUE)

  ## check dimensions
  expect_equal(nn$k, k)
  expect_equal(dim(nn$dist), c(nrow(x), k))
  expect_equal(dim(nn$id), c(nrow(x), k))
  expect_type(nn$id, "integer")

  ## check visually
  #plot(x)
  #points(x[nn$id[1,],], col="red", lwd=5)
  #points(x[nn$id[2,],], col="green", lwd=5)

  ## compare with kNN found using distances
  nn_d <- sNN(dist(x), k, sort = TRUE)

  ## check visually
  #plot(x)
  #points(x[nn_d$id[1,],], col="red", lwd=5)
  #points(x[nn_d$id[2,],], col="green", lwd=5)

  ### will aggree minus some tries
  expect_equal(nn, nn_d)

  ## calculate dist internally
  nn_d2 <- sNN(x, k, search = "dist", sort = TRUE)
  expect_equal(nn, nn_d2)

  ## missing values, but distances are fine
  x_na <- x
  x_na[c(1,3,5), 1] <- NA
  expect_error(sNN(x_na, k = 3), regexp = "NA")
  res_d1 <- sNN(x_na, k = 3, search = "dist")
  res_d2 <- sNN(dist(x_na), k = 3)
  expect_equal(res_d1, res_d2)

  ## introduce NAs into dist
  x_na[c(1,3,5),] <- NA
  expect_error(sNN(x_na, k = 3), regexp = "NA")
  expect_error(sNN(x_na, k = 3, search = "dist"), regexp = "NA")
  expect_error(sNN(dist(x_na), k = 3), regexp = "NA")


  ## sort and kNN to reduce k
  nn10 <- sNN(x, k = 10, sort = FALSE)
  expect_false(nn10$sort_shared)
  expect_error(sNN(nn10, k = 11))

  nn5 <- sNN(nn10, k = 5, sort = TRUE)
  nn5_x <- sNN(x, k = 5, sort = TRUE)
  expect_equal(nn5, nn5_x)

  ## test with simple data
  x <- data.frame(x = 1:10, check.names = FALSE)
  nn <- sNN(x, k = 5)

  i <- 1
  j_ind <- 1
  j <- nn$id[i,j_ind]
  intersect(c(i, nn$id[i,]), nn$id[j,])
  nn$shared[i,j_ind]

  # compute the sNN simularity in R
  ss <- matrix(nrow = nrow(x), ncol = nn$k)
  for(i in seq_len(nrow(x)))
    for(j_ind in 1:nn$k)
      ss[i, j_ind] <- length(intersect(c(i, nn$id[i,]), nn$id[nn$id[i,j_ind],]))

  expect_equal(nn$shared, ss)
})

test_that("sNN sorting and arguments", {
  x <- as.matrix(iris[!duplicated(iris[, -5]), -5])
  k <- 6

  nn <- sNN(x, k = k)
  expect_true(nn$sort_shared)
  expect_false(nn$sort)
  for (i in seq_len(nrow(x)))
    expect_identical(order(-nn$shared[i, ], nn$id[i, ]), seq_len(k))

  ## increasing order (ties by id)
  inc <- sort(nn, decreasing = FALSE)
  expect_false(inc$sort_shared)
  for (i in seq_len(nrow(x)))
    expect_identical(order(inc$shared[i, ], inc$id[i, ]), seq_len(k))
  expect_equal(unname(inc$shared), unname(t(apply(nn$shared, 1, sort))))

  ## sorting back gives the original order
  expect_identical(sort(inc), nn)

  ## k is required unless x is a kNN object
  expect_error(sNN(x), "k is required")
  expect_identical(sNN(kNN(x, 4))$k, 4L)

  ## kt can be equal to k, but not larger
  expect_no_error(sNN(x, k = k, kt = k))
  expect_error(sNN(x, k = k, kt = k + 1), "at most k")

  ## thresholded edges are NA
  t3 <- sNN(x, k = k, kt = 3)
  expect_true(all(is.na(t3$id) == is.na(t3$shared)))
  expect_true(all(t3$shared[!is.na(t3$shared)] >= 3))
})
