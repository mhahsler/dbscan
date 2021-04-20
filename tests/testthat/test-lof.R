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
# DMwR is now retired
#if(requireNamespace("DMwR", quietly = TRUE)) {
#  system.time(lof_DMwr <- DMwR::lofactor(x, k=4))
#
#  expect_equal(lof, lof_DMwr)
#}

## missing values, but distances are fine
x_na <- x
x_na[c(1,3,5), 1] <- NA
expect_error(lof(x_na, k=4), regexp = "NA")
res_d1 <- lof(x_na, k =4, search = "dist")
res_d2 <- lof(dist(x_na), k = 4)
expect_equal(res_d1, res_d2)

x_na[c(1,3,5), 2] <- NA
expect_error(lof(x_na, k=4), regexp = "NA")
expect_error(lof(x_na, k=4, search = "dist"),
  regexp = "NA")
expect_error(lof(dist(x_na), k =4), regexp = "NA")

## test with tied distances
x <- rbind(1,2,3,4,5,6,7)
expect_equal(round(lof(x, k = 3), 7),
  c(1.0679012, 1.0679012, 1.0133929, 0.8730159, 1.0133929, 1.0679012, 1.0679012))

expect_equal(round(lof(dist(x), k = 3),7),
  c(1.0679012, 1.0679012, 1.0133929, 0.8730159, 1.0133929, 1.0679012, 1.0679012))

