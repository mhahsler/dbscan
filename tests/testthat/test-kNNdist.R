library("testthat")
library("dbscan")

context("kNNdist")

set.seed(665544)
n <- 1000
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
)


d <- kNNdist(x, k = 5)
expect_equal(length(d), n)

d <- kNNdist(x, k = 5, all = TRUE)
expect_equal(dim(d), c(n, 5))

# does the plot work?
kNNdistplot(x, 5)

