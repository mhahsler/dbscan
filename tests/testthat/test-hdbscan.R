test_that("HDBSCAN", {
  data("iris")

  ## minPts not given
  expect_error(hdbscan(iris))

  ## Expects numerical data; species is factor
  expect_error(dbscan(iris, minPts = 4))

  iris <- as.matrix(iris[,1:4])

  res <- hdbscan(iris, minPts = 4)
  expect_length(res$cluster, nrow(iris))

  ## expected result of table(res$cluster) is:
  expect_identical(table(res$cluster, dnn = NULL),
                    as.table(c("1" = 100L, "2" = 50L)))

  ## compare on moons data
  data("moons")
  res <- hdbscan(moons, minPts = 5)
  expect_length(res$cluster, nrow(moons))

  ## Check hierarchy matches dbscan* at every value
  check <- rep(FALSE, nrow(moons)-1)
  core_dist <- kNNdist(moons, k=5-1)

  ## cutree doesn't distinguish noise as 0, so we make a new method to do it manually
  cut_tree <- function(hcl, eps, core_dist){
    cuts <- unname(cutree(hcl, h=eps))
    cuts[which(core_dist > eps)] <- 0 # Use core distance to distinguish noise
    cuts
  }

  eps_values <- sort(res$hc$height, decreasing = TRUE)+.Machine$double.eps ## Machine eps for consistency between cuts
  for (i in seq_along(eps_values)) {
    cut_cl <- cut_tree(res$hc, eps_values[i], core_dist)
    dbscan_cl <- dbscan(moons, eps = eps_values[i], minPts = 5, borderPoints = FALSE) # DBSCAN* doesn't include border points

    ## Use run length encoding as an ID-independent way to check ordering
    check[i] <- (all.equal(rle(cut_cl)$lengths, rle(dbscan_cl$cluster)$lengths) == "TRUE")
  }

  expect_true(all(check))

  ## Expect generating extra trees doesn't fail
  res <- hdbscan(moons, minPts = 5, gen_hdbscan_tree = TRUE, gen_simplified_tree = TRUE)
  expect_s3_class(res, "hdbscan")

  ## Expect hdbscan tree matches stats:::as.dendrogram version of hclust object
  hc_dend <- as.dendrogram(res$hc)
  expect_s3_class(hc_dend, "dendrogram")
  expect_identical(hc_dend, res$hdbscan_tree)

  ## Expect hdbscan works with non-euclidean distances
  dist_moons <- dist(moons, method = "canberra")
  res <- hdbscan(dist_moons, minPts = 5)
  expect_s3_class(res, "hdbscan")
})

test_that("mrdist", {
  expect_identical(mrdist(cbind(1:10), 2),  mrdist(dist(cbind(1:10)), 2))
  expect_identical(mrdist(cbind(1:11), 3), mrdist(dist(cbind(1:11)), 3))
})

test_that("HDBSCAN(e)", {
  X <- data.frame(
   x = c(
    0.08, 0.46, 0.46, 2.95, 3.50, 1.49, 6.89, 6.87, 0.21, 0.15,
    0.15, 0.39, 0.80, 0.80, 0.37, 3.63, 0.35, 0.30, 0.64, 0.59, 1.20, 1.22,
    1.42, 0.95, 2.70, 6.36, 6.36, 6.36, 6.60, 0.04, 0.71, 0.57, 0.24, 0.24,
    0.04, 0.04, 1.35, 0.82, 1.04, 0.62, 0.26, 5.98, 1.67, 1.67, 0.48, 0.15,
    6.67, 6.67, 1.20, 0.21, 3.99, 0.12, 0.19, 0.15, 6.96, 0.26, 0.08, 0.30,
    1.04, 1.04, 1.04, 0.62, 0.04, 0.04, 0.04, 0.82, 0.82, 1.29, 1.35, 0.46,
    0.46, 0.04, 0.04, 5.98, 5.98, 6.87, 0.37, 6.47, 6.47, 6.47, 6.67, 0.30,
    1.49, 3.21, 3.21, 0.75, 0.75, 0.46, 0.46, 0.46, 0.46, 3.63, 0.39, 3.65,
    4.09, 4.01, 3.36, 1.43, 3.28, 5.94, 6.35, 6.87, 5.60, 5.99, 0.12, 0.00,
    0.32, 0.39, 0.00, 1.63, 1.36, 5.67, 5.60, 5.79, 1.10, 2.99, 0.39, 0.18
    ),
   y = c(
    7.41, 8.01, 8.01, 5.44, 7.11, 7.13, 1.83, 1.83, 8.22, 8.08,
    8.08, 7.20, 7.83, 7.83, 8.29, 5.99, 8.32, 8.22, 7.38, 7.69, 8.22, 7.31,
    8.25, 8.39, 6.34, 0.16, 0.16, 0.16, 1.66, 7.55, 7.90, 8.18, 8.32, 8.32,
    7.97, 7.97, 8.15, 8.43, 7.83, 8.32, 8.29, 1.03, 7.27, 7.27, 8.08, 7.27,
    0.79, 0.79, 8.22, 7.73, 6.62, 7.62, 8.39, 8.36, 1.73, 8.29, 8.04, 8.22,
    7.83, 7.83, 7.83, 8.32, 8.11, 7.69, 7.55, 7.20, 7.20, 8.01, 8.15, 7.55,
    7.55, 7.97, 7.97, 1.03, 1.03, 1.24, 7.20, 0.47, 0.47, 0.47, 0.79, 8.22,
    7.13, 6.48, 6.48, 7.10, 7.10, 8.01, 8.01, 8.01, 8.01, 5.99, 8.04, 5.22,
    5.82, 5.14, 4.81, 7.62, 5.73, 0.55, 1.31, 0.05, 0.95, 1.59, 7.99, 7.48,
    8.38, 7.12, 2.01, 1.40, 0.00, 9.69, 9.47, 9.25, 2.63, 6.89, 0.56, 3.11
   )
  )

  hdbe <- hdbscan(X, minPts = 3, cluster_selection_epsilon = 1)
  #plot(X, col = hdbe$cluster + 1L, main = "HDBSCAN(e)")

  expect_equal(ncluster(hdbe), 5L)
  expect_equal(nnoise(hdbe), 0L)
})

