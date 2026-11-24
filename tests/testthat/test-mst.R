test_that("mst_dist", {
  set.seed(1234)
  x <- cbind(rnorm(60), rnorm(60))
  rownames(x) <- paste0("p", 1:60)
  d <- dist(x)

  hc <- mst_dist(d)
  ref <- hclust(d, method = "single")

  expect_s3_class(hc, "hclust")
  expect_equal(hc$height, ref$height)
  expect_equal(hc$labels, ref$labels)
  expect_equal(hc$method, "robust single")
  expect_equal(hc$dist.method, "euclidean")
  expect_false(is.unsorted(hc$height))
  expect_equal(sort(hc$order), 1:60)

  # same hierarchy as hclust (continuous data: no ties)
  for (k in c(1, 2, 5, 20, 60))
    expect_equal(cutree(hc, k), cutree(ref, k))

  # merge matrix passes the base R validity checks used by as.dendrogram
  expect_no_error(dend <- as.dendrogram(hc))
  expect_equal(nobs(dend), 60)
  expect_setequal(labels(dend), rownames(x))
  expect_no_error(plot(hc))

  # merge rows keep the MST order and are not normalized to the hclust
  # convention, but the hierarchy is valid: every observation and every merge
  # step is used exactly once.
  mg <- hc$merge
  expect_equal(sort(-mg[mg < 0]), 1:60)
  expect_equal(sort(mg[mg > 0]), 1:58)
  # clusters are only merged after they were created
  expect_true(all(mg[mg > 0] < row(mg)[mg > 0]))

  # smallest case and invalid input
  expect_setequal(mst_dist(dist(c(0, 3)))$merge, c(-1L, -2L))
  expect_error(mst_dist(dist(1)))
  expect_error(mst_dist(as.matrix(d)))
  d_na <- d
  d_na[3] <- NA
  expect_error(mst_dist(d_na), "finite")
})

test_that("mrdist keeps labels", {
  set.seed(1234)
  x <- cbind(rnorm(30), rnorm(30))
  rownames(x) <- paste0("p", 1:30)

  expect_equal(attr(mrdist(x, minPts = 3), "Labels"), rownames(x))
  expect_equal(attr(mrdist(dist(x), minPts = 3), "Labels"), rownames(x))
  expect_equal(labels(mrdist(x, minPts = 3)), rownames(x))
  expect_equal(mst_dist(mrdist(x, minPts = 3))$labels, rownames(x))

  # no labels stay NULL
  expect_null(attr(mrdist(unname(x), minPts = 3), "Labels"))
})

test_that("extractFOSC cluster_selection_epsilon", {
  data("moons")
  for (eps in c(0, 0.2, 1)) {
    h <- hdbscan(moons, minPts = 5, cluster_selection_epsilon = eps)
    cl <- extractFOSC(mrdist(moons, 5) |> mst_dist(), minPts = 5,
                      cluster_selection_epsilon = eps)$cluster
    # same partition (cluster ids may differ)
    expect_equal(as.integer(factor(h$cluster)) , as.integer(factor(cl)))
    expect_equal(extractFOSC(h, cluster_selection_epsilon = eps)$cluster, cl)
  }
  # default unchanged
  h <- hdbscan(moons, minPts = 5)
  expect_equal(extractFOSC(h)$cluster, extractFOSC(h, cluster_selection_epsilon = 0)$cluster)
  expect_error(extractFOSC(h, cluster_selection_epsilon = -1), "cluster_selection_epsilon")
  expect_error(extractFOSC(h, cluster_selection_epsilon = c(1, 2)), "cluster_selection_epsilon")
})

test_that("C++ mst_mrd matches mst_dist_cpp", {
  set.seed(1234)
  for (n in c(2, 3, 10, 60)) {
    x <- cbind(rnorm(n), rnorm(n), rnorm(n))
    rownames(x) <- paste0("p", seq_len(n))
    cd <- if (n > 3) coredist(x, 3) else runif(n)
    m1 <- dbscan:::mst_dist_cpp(mrdist(x, minPts = 2, coredist = cd), n)
    m2 <- dbscan:::mst_mrd(x, cd)
    expect_equal(m2, m1)
  }

  # ties (duplicated points) give the same tree
  x <- rbind(matrix(0, 5, 2), matrix(1, 5, 2), c(5, 5))
  cd <- rep(0.5, nrow(x))
  expect_equal(dbscan:::mst_mrd(x, cd),
               dbscan:::mst_dist_cpp(mrdist(x, minPts = 2, coredist = cd), nrow(x)))

  expect_error(dbscan:::mst_mrd(x[1, , drop = FALSE], 1), "at least 2")
  expect_error(dbscan:::mst_mrd(x, 1:3), "coredist")
})

test_that("hdbscan with data and dist input gives the same result", {
  data("moons")
  h1 <- hdbscan(moons, minPts = 5)
  h2 <- hdbscan(dist(moons), minPts = 5)
  expect_equal(h1$cluster, h2$cluster)
  expect_equal(h1$hc$merge, h2$hc$merge)
  expect_equal(h1$hc$height, h2$hc$height)
  expect_equal(h1$outlier_scores, h2$outlier_scores)
  expect_s3_class(h1$hc, "hclust")

  x <- as.matrix(moons)
  rownames(x) <- paste0("p", seq_len(nrow(x)))
  expect_equal(hdbscan(x, minPts = 5)$hc$labels, rownames(x))
})

test_that("mst_mrdist", {
  set.seed(1234)
  x <- cbind(rnorm(50), rnorm(50))
  rownames(x) <- paste0("p", 1:50)

  hc <- mst_mrdist(x, minPts = 4)
  ref <- mst_dist(mrdist(x, minPts = 4))
  expect_s3_class(hc, "hclust")
  expect_equal(hc$merge, ref$merge)
  expect_equal(hc$height, ref$height)
  expect_equal(hc$order, ref$order)
  expect_equal(hc$labels, rownames(x))
  expect_equal(hc$minPts, 4L)
  expect_equal(hc$dist.method, "mutual reachability (euclidean)")

  # given core distances and data frames
  cd <- coredist(x, 4)
  expect_equal(mst_mrdist(x, 4, coredist = cd)$merge, hc$merge)
  expect_equal(mst_mrdist(as.data.frame(x), 4)$merge, hc$merge)

  # works with extractFOSC
  expect_equal(extractFOSC(hc)$cluster, extractFOSC(ref, minPts = 4)$cluster)

  expect_error(mst_mrdist(x, minPts = 1))
  expect_error(mst_mrdist(x[1, , drop = FALSE], minPts = 2), "at least 2")
  expect_error(mst_mrdist(x, 4, coredist = 1:3), "coredist")
  x_na <- x; x_na[1, 1] <- NA
  expect_error(mst_mrdist(x_na, 4))
})

test_that("minPts is optional if coredist is given", {
  set.seed(1234)
  x <- cbind(rnorm(40), rnorm(40))
  cd <- coredist(x, 4)

  expect_equal(mrdist(x, coredist = cd), mrdist(x, minPts = 4, coredist = cd))
  expect_equal(mrdist(dist(x), coredist = cd), mrdist(x, minPts = 4, coredist = cd),
               ignore_attr = "call")
  expect_error(mrdist(x), "minPts is required")
  expect_error(mrdist(x, minPts = 1))

  hc <- mst_mrdist(x, coredist = cd)
  expect_equal(hc$merge, mst_mrdist(x, 4, coredist = cd)$merge)
  expect_null(hc$minPts)
  expect_equal(mst_mrdist(x, 4, coredist = cd)$minPts, 4L)
  expect_error(mst_mrdist(x), "minPts is required")
  expect_error(mst_mrdist(x, minPts = 1, coredist = cd))
})

test_that("lower-triangle index does not overflow for long dist objects", {
  # closed form in double arithmetic (exact below 2^53)
  pos0 <- function(n, i, j) {
    n <- as.double(n)
    a <- as.double(pmin(i, j)); b <- as.double(pmax(i, j))
    n * a - (a + 1) * a / 2 + b - a - 1
  }
  n <- 100000L
  i <- c(0L, 5L, 46340L, 46341L, 50000L, 70000L, 99998L)
  j <- c(1L, 99999L, 46341L, 46342L, 99999L, 60000L, 99999L)

  expect_equal(dbscan:::lt_pos_cpp(n, i, j, FALSE), pos0(n, i, j))
  expect_equal(dbscan:::lt_pos_cpp(n, j, i, FALSE), pos0(n, i, j))
  expect_equal(dbscan:::lt_pos_cpp(n, i + 1L, j + 1L, TRUE), pos0(n, i, j))

  # last element of the dist vector
  expect_equal(dbscan:::lt_pos_cpp(n, n - 2L, n - 1L, FALSE), as.double(n) * (n - 1) / 2 - 1)

  # agrees with stats::dist for a small object
  d <- dist(matrix(rnorm(40), ncol = 2))
  m <- as.matrix(d)
  ii <- c(1L, 3L, 7L, 20L); jj <- c(2L, 9L, 20L, 4L)
  expect_equal(unclass(d)[dbscan:::lt_pos_cpp(20L, ii - 1L, jj - 1L, FALSE) + 1],
               m[cbind(ii, jj)])
})

test_that("mst_dist checks the size of the dist object", {
  d <- dist(matrix(rnorm(20), ncol = 2))
  expect_error(dbscan:::mst_dist_cpp(unclass(d), 11), "does not match")
})

test_that("dist_subset works", {
  set.seed(1)
  d <- dist(matrix(rnorm(60), ncol = 2))
  idx <- c(2L, 5L, 11L, 30L)
  expect_equal(as.vector(dbscan:::dist_subset(d, idx)),
               as.vector(as.dist(as.matrix(d)[idx, idx])))
})
