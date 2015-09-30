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
expect_identical(nn$eps, eps)
expect_identical(length(nn$dist), nrow(x))
expect_identical(length(nn$id), nrow(x))

expect_identical(sapply(nn$dist, length), sapply(nn$id, length))

## check visually
#plot(x)
#points(x[nn$id[[1]],], col="red", lwd=5)
#points(x[nn$id[[2]],], col="green", lwd=5)
#points(x[1:2,, drop = FALSE], col="blue", pch="+", cex=2)

## compare with manually found NNs
d <- as.matrix(dist(x)); diag(d) <- Inf
ids <- apply(d, MARGIN = 1, FUN =
    function(y) {
      o <- order(y, decreasing = FALSE)
      o[y[o] < eps]
    }
)

dists <- lapply(1:nrow(d), FUN =
    function(i) {
      unname(d[i,ids[[i]]])
    }
)
names(dists) <- rownames(x)

## check visually
#plot(x)
#points(x[ids[[1]],], col="red", lwd=5)
#points(x[ids[[2]],], col="green", lwd=5)
#points(x[1:2,, drop = FALSE], col="blue", pch="+", cex=2)

#head(ids)
#head(nn$id)
#head(dists, n=2)
#head(nn$dist, n=2)

expect_identical(nn$id, ids)
expect_identical(nn$dist, dists)


## without sorting
nn2 <- dbscan::frNN(x, eps = eps, sort = FALSE)
expect_identical(lapply(nn$id, sort),
  lapply(nn2$id, sort))

## search options
nn_linear <- dbscan::frNN(x, eps=eps, search = "linear")
expect_identical(nn, nn_linear)

## split options
for(so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
  nn3 <- dbscan::frNN(x, eps=eps, splitRule = so)
  expect_identical(nn, nn3)
}

## bucket size
for(bs in c(5, 10, 15, 100)) {
  nn3 <- dbscan::frNN(x, eps=eps, bucketSize = bs)
  expect_identical(nn, nn3)
}


## add 100 copied points to check if self match filtering works
x <- rbind(x, x[sample(1:nrow(x), 100),])
rownames(x) <- paste("Object_", 1:nrow(x), sep="")

eps <- .5
nn <- dbscan::frNN(x, eps = eps, sort = TRUE)

## compare with manually found NNs
d <- as.matrix(dist(x)); diag(d) <- Inf
ids <- apply(d, MARGIN = 1, FUN =
    function(y) {
      o <- order(y, decreasing = FALSE)
      o[y[o] < eps]
    }
)

dists <- lapply(1:nrow(d), FUN =
    function(i) {
      unname(d[i,ids[[i]]])
    }
)
names(dists) <- rownames(x)

expect_identical(nn$dist, dists)
expect_identical(nn$id, ids)
