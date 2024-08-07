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
