library("testthat")
library("dbscan")

context("frNN")

set.seed(665544)
n <- 1000
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  z = runif(10, 0, 10) + rnorm(n, sd = 0.2)
)

## no duplicates first!
x <- x[!duplicated(x),]


rownames(x) <- paste("Object_", 1:nrow(x), sep="")

eps <- .5
nn <- dbscan::frNN(x, eps = eps, sort = TRUE)

## check dimensions
expect_equal(nn$eps, eps)
expect_equal(length(nn$dist), nrow(x))
expect_equal(length(nn$id), nrow(x))

expect_equal(sapply(nn$dist, length), sapply(nn$id, length))

## check visually
#plot(x)
#points(x[nn$id[[1]],], col="red", lwd=5)
#points(x[nn$id[[2]],], col="green", lwd=5)
#points(x[1:2,, drop = FALSE], col="blue", pch="+", cex=2)

## compare with manually found NNs
nn_d <- dbscan::frNN(dist(x), eps = eps, sort = TRUE)
expect_equal(nn, nn_d)

nn_d2 <- dbscan::frNN(x, eps = eps, sort = TRUE, search = "dist")
expect_equal(nn, nn_d2)

## without sorting
nn2 <- dbscan::frNN(x, eps = eps, sort = FALSE)
expect_equal(lapply(nn$id, sort),
  lapply(nn2$id, sort))

## search options
nn_linear <- dbscan::frNN(x, eps=eps, search = "linear")
expect_equal(nn, nn_linear)

## split options
for(so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
  nn3 <- dbscan::frNN(x, eps=eps, splitRule = so)
  expect_equal(nn, nn3)
}

## bucket size
for(bs in c(5, 10, 15, 100)) {
  nn3 <- dbscan::frNN(x, eps=eps, bucketSize = bs)
  expect_equal(nn, nn3)
}


## add 100 copied points to check if self match filtering works
x <- rbind(x, x[sample(1:nrow(x), 100),])
rownames(x) <- paste("Object_", 1:nrow(x), sep="")

eps <- .5
nn <- dbscan::frNN(x, eps = eps, sort = TRUE)

## compare with manually found NNs
nn_d <- dbscan::frNN(x, eps = eps, sort = TRUE, search = "dist")

expect_equal(nn, nn_d)
