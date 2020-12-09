library("dbscan")
library("testthat")

context("predict")

set.seed(665544)
n <- 100
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
)

res <- dbscan(x, eps = .3, minPts = 3)

# check if points with a little noise are assigned to the same cluster
newdata <- x + rnorm(2*n, 0, .1)
pr <- predict(res, newdata, data = x)

expect_equal(res$cluster, pr)



