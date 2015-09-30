library("testthat")
library("dbscan")

context("kNN")

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

k <- 5L
nn <- dbscan::kNN(x, k=k, sort = TRUE)

## check dimensions
expect_identical(nn$k, k)
expect_identical(dim(nn$dist), c(nrow(x), k))
expect_identical(dim(nn$id), c(nrow(x), k))

## check visually
#plot(x)
#points(x[nn$id[1,],], col="red", lwd=5)
#points(x[nn$id[2,],], col="green", lwd=5)

## compare with manually found NNs
d <- as.matrix(dist(x)); diag(d) <- Inf
ids <- t(apply(d, MARGIN = 1, order, decreasing = FALSE))[,1:k]
dimnames(ids) <- list(rownames(x), 1:k)
dists <- t(apply(d, MARGIN = 1, sort, decreasing = FALSE))[,1:k]
dimnames(dists) <- list(rownames(x), 1:k)

## check visually
#plot(x)
#points(x[ids[1,],], col="red", lwd=5)
#points(x[ids[2,],], col="green", lwd=5)

#head(ids)
#head(nn$id)

expect_identical(nn$id, ids)
expect_identical(nn$dist, dists)


## without sorting
nn2 <- dbscan::kNN(x, k=k, sort = FALSE)
expect_identical(t(apply(nn$id, MARGIN = 1, sort)),
  t(apply(nn2$id, MARGIN = 1, sort)))

## search options
nn_linear <- dbscan::kNN(x, k=k, search = "linear")
expect_identical(nn, nn_linear)

## split options
for(so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
  nn3 <- dbscan::kNN(x, k=k, splitRule = so)
  expect_identical(nn, nn3)
}

## bucket size
for(bs in c(5, 10, 15, 100)) {
  nn3 <- dbscan::kNN(x, k=k, bucketSize = bs)
  expect_identical(nn, nn3)
}

## the order is not stable with matching distances which means that the
## k-NN are not stable. add 100 copied points to check if self match
## filtering and sort works
x <- rbind(x, x[sample(1:nrow(x), 100),])
rownames(x) <- paste("Object_", 1:nrow(x), sep="")

k <- 5L
nn <- dbscan::kNN(x, k=k, sort = TRUE)

## compare with manually found NNs
d <- as.matrix(dist(x)); diag(d) <- Inf
ids <- t(apply(d, MARGIN = 1, order, decreasing = FALSE))[,1:k]
dimnames(ids) <- list(rownames(x), 1:k)
dists <- t(apply(d, MARGIN = 1, sort, decreasing = FALSE))[,1:k]
dimnames(dists) <- list(rownames(x), 1:k)

expect_identical(nn$dist, dists)
## This is expected to fail: expect_identical(nn$id, ids)

## only compare the ids which do not have unique distances
ms <- nn$id == ids
cos <- diff(dists)==0; cos[,1] <- FALSE
ms[cos] <- TRUE
expect_false(any(rowSums(ms[,-k]) !=4))

