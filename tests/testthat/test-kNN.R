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
expect_equal(nn$k, k)
expect_equal(dim(nn$dist), c(nrow(x), k))
expect_equal(dim(nn$id), c(nrow(x), k))

## check visually
#plot(x)
#points(x[nn$id[1,],], col="red", lwd=5)
#points(x[nn$id[2,],], col="green", lwd=5)

## compare with kNN found using distances
nn_d <- kNN(dist(x), k)

## check visually
#plot(x)
#points(x[nn_d$id[1,],], col="red", lwd=5)
#points(x[nn_d$id[2,],], col="green", lwd=5)

### will aggree minus some tries
expect_equal(nn, nn_d)

## calculate dist internally
nn_d2 <- kNN(x, k, search = "dist")
expect_equal(nn, nn_d2)

## without sorting
nn2 <- dbscan::kNN(x, k=k, sort = FALSE)
expect_equal(t(apply(nn$id, MARGIN = 1, sort)),
  t(apply(nn2$id, MARGIN = 1, sort)))

## search options
nn_linear <- dbscan::kNN(x, k=k, search = "linear")
expect_equal(nn, nn_linear)

## split options
for(so in c("STD", "MIDPT", "FAIR", "SL_FAIR")) {
  nn3 <- dbscan::kNN(x, k=k, splitRule = so)
  expect_equal(nn, nn3)
}

## bucket size
for(bs in c(5, 10, 15, 100)) {
  nn3 <- dbscan::kNN(x, k=k, bucketSize = bs)
  expect_equal(nn, nn3)
}

## the order is not stable with matching distances which means that the
## k-NN are not stable. add 100 copied points to check if self match
## filtering and sort works
x <- rbind(x, x[sample(1:nrow(x), 100),])
rownames(x) <- paste("Object_", 1:nrow(x), sep="")

k <- 5L
nn <- dbscan::kNN(x, k=k, sort = TRUE)

## compare with manually found NNs
nn_d <- dbscan::kNN(x, k=k, search = "dist")

expect_equal(nn$dist, nn_d$dist)
## This is expected to fail: expect_equal(nn$id, ids)

## only compare the ids which do not have unique distances
ms <- nn$id == nn_d$id
cos <- diff(nn_d$dist)==0; cos[,1] <- FALSE
ms[cos] <- TRUE
expect_false(any(rowSums(ms[,-k]) !=4))

## missing values, but distances are fine
x_na <- x
x_na[c(1,3,5), 1] <- NA
expect_error(kNN(x_na, k = 3), regexp = "NA")
res_d1 <- kNN(x_na, k = 3, search = "dist")
res_d2 <- kNN(dist(x_na), k = 3)
expect_equal(res_d1, res_d2)

## introduce NAs into dist
x_na[c(1,3,5),] <- NA
expect_error(kNN(x_na, k = 3), regexp = "NA")
expect_error(kNN(x_na, k = 3, search = "dist"), regexp = "NA")
expect_error(kNN(dist(x_na), k = 3), regexp = "NA")
