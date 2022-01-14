library("dbscan")
library("testthat")

context("predict")

set.seed(3)
n <- 100
x_data <- cbind(
  x = runif(5, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(5, 0, 10) + rnorm(n, sd = 0.2)
)

x_noise <- cbind(
  x = runif(n/2, 0, 10),
  y = runif(n/2, 0, 10)
)

x <- rbind(x_data, x_noise)

# check if l points with a little noise are assigned to the same cluster
l <- 20
newdata <- rbind(
  x_data[1:l,] + rnorm(2*l, 0, .05),
  x_noise[1:l,] + rnorm(2*l, 0, .05)
)

idx <- c(1:l, n + (1:l))

#plot(x, col = rep(c("black", "gray"), each = n))
#points(newdata, col = rep(c("red", "gray"), each = l), pch = 16)

# DBSCAN
res <- dbscan(x, eps = .3, minPts = 3)
pr <- predict(res, newdata, data = x)

rbind(true = res$cluster[idx], pred = pr)
expect_equal(res$cluster[idx], pr)
#plot(x, col = ifelse(res$cluster == 0, "gray", res$cluster))
#points(newdata, col = ifelse(pr == 0, "gray", pr), pch = 16)

# OPTICS
res <- optics(x, minPts = 3)
res <- extractDBSCAN(res, eps = .3)
pr <- predict(res, newdata, data = x)

rbind(true = res$cluster[idx], pred = pr)
expect_equal(res$cluster[idx], pr)

# currently no implementation for extractXi

# HDBSCAN (note predict is not perfect for the data.)
res <- hdbscan(x, minPts = 3)
pr <- predict(res, newdata, data = x)

rbind(true = res$cluster[idx], pred = pr)
accuracy <- sum(res$cluster[idx] == pr)/length(pr)
expect_true(accuracy > .9)

# show misclassifications
#plot(x, col = ifelse(res$cluster == 0, "gray", res$cluster))
#points(newdata, col = ifelse(pr == 0, "gray", pr), pch = 16)
#points(newdata[res$cluster[idx] != pr,, drop = FALSE], col = "red", pch = 4, lwd = 2)



