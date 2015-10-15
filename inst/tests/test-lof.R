library("dbscan")
library("testthat")


context("LOF")

set.seed(665544)
n <- 600
x <- cbind(
  x=runif(10, 0, 5) + rnorm(n, sd=0.4),
  y=runif(10, 0, 5) + rnorm(n, sd=0.4)
)

### calculate LOF score
system.time(lof <- lof(x, k=4))
expect_identical(length(lof), nrow(x))

system.time(lof_d <- lof(dist(x), k=4))
expect_equal(lof, lof_d)

## compare with lofactor from DMwR
if(requireNamespace("DMwR", quietly = TRUE)) {
  system.time(lof_DMwr <- DMwR::lofactor(x, k=4))

  expect_equal(lof, lof_DMwr)
}
