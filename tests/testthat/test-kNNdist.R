test_that("kNNdist", {
  set.seed(665544)
  n <- 1000
  x <- cbind(
    x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
    z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
  )

  d <- kNNdist(x, k = 5)
  expect_length(d, n)

  d <- kNNdist(x, k = 5, all = TRUE)
  expect_equal(dim(d), c(n, 5))

  # kNNdistplot(): kNNdist_param goes to kNNdist, ... goes to the plot
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(kNNdistplot(x, k = 4))
  expect_no_error(kNNdistplot(x, minPts = 5, col = "red", main = "test",
    ylab = "custom", kNNdist_param = list(search = "linear")))
  expect_no_error(kNNdistplot(x, k = 1:3, col = 1:3,
    kNNdist_param = list(search = "linear")))
  expect_no_error(kNNdistplot(dist(x[1:50, ]), k = 3))
})
