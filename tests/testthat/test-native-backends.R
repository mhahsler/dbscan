test_that("randomized exact searches agree across native and distance backends", {
  skip_on_cran()
  set.seed(20261027)
  iterations <- as.integer(Sys.getenv("DBSCAN_RANDOMIZED_ITERATIONS", "20"))
  stopifnot(!is.na(iterations), iterations >= 1L)

  same_partition <- function(a, b) {
    expect_identical(a == 0L, b == 0L)
    expect_identical(outer(a, a, "=="), outer(b, b, "=="))
  }

  for (iteration in seq_len(iterations)) {
    n <- sample(12:60, 1)
    dimensions <- sample(1:6, 1)
    x <- matrix(rnorm(n * dimensions), ncol = dimensions)
    distances <- dist(x)
    eps <- runif(1, .2, 2) * sqrt(dimensions)
    min_pts <- sample(2:min(8, n), 1)
    k <- sample(c(1L, min_pts - 1L, n - 1L), 1)
    bucket <- sample(c(1L, 3L, 10L), 1)
    rule <- sample(c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST"), 1)

    reference_knn <- kNN(distances, k = k)
    reference_radius <- frNN(distances, eps = eps)
    reference_lof <- lof(distances, minPts = min_pts)
    reference_density <- pointdensity(distances, eps = eps)
    weights <- runif(n, -.5, 2)
    reference_star <- dbscan(distances, eps, min_pts, borderPoints = FALSE)
    reference_weighted <- dbscan(distances, eps, min_pts, weights = weights,
                                borderPoints = FALSE)
    reference_clusters <- dbscan(distances, eps, min_pts)
    reference_optics <- optics(distances, eps = eps, minPts = min_pts)

    for (search in c("kdtree", "linear", "dist")) {
      nearest <- kNN(x, k, search = search, bucketSize = bucket, splitRule = rule)
      expect_equal(nearest$id, reference_knn$id)
      expect_equal(nearest$dist, reference_knn$dist, tolerance = 1e-12)
      radius <- frNN(x, eps, search = search, bucketSize = bucket, splitRule = rule)
      expect_equal(radius$id, reference_radius$id)
      expect_equal(radius$dist, reference_radius$dist, tolerance = 1e-12)
      expect_equal(lof(x, min_pts, search = search, bucketSize = bucket,
                       splitRule = rule), reference_lof, tolerance = 1e-10)
      expect_equal(pointdensity(x, eps, search = search, bucketSize = bucket,
                                splitRule = rule), reference_density)
      star <- dbscan(x, eps, min_pts, borderPoints = FALSE, search = search,
                     bucketSize = bucket, splitRule = rule)
      same_partition(star$cluster, reference_star$cluster)
      weighted <- dbscan(x, eps, min_pts, weights = weights, borderPoints = FALSE,
                         search = search, bucketSize = bucket, splitRule = rule)
      same_partition(weighted$cluster, reference_weighted$cluster)
      clusters <- dbscan(x, eps, min_pts, search = search, bucketSize = bucket,
                         splitRule = rule)
      # A border point can legitimately belong to either of two clusters.
      expect_identical(clusters$cluster == 0L, reference_clusters$cluster == 0L)
      core <- reference_star$cluster != 0L
      same_partition(clusters$cluster[core], reference_clusters$cluster[core])
      ordering <- optics(x, eps, min_pts, search = search, bucketSize = bucket,
                         splitRule = rule)
      expect_equal(ordering$order, reference_optics$order)
      expect_equal(ordering$coredist, reference_optics$coredist, tolerance = 1e-12)
      expect_equal(ordering$reachdist, reference_optics$reachdist, tolerance = 1e-12)
    }

    # Explicit queries have no self-exclusion, and the distance backend does
    # not accept them. Use cross-distances computed independently by stats.
    query <- matrix(rnorm(4 * dimensions), ncol = dimensions)
    cross <- as.matrix(dist(rbind(query, x)))[seq_len(nrow(query)),
                                             nrow(query) + seq_len(n), drop = FALSE]
    expected_order <- t(apply(cross, 1, order))[, seq_len(k), drop = FALSE]
    expected_distance <- matrix(cross[cbind(rep(seq_len(nrow(query)), each = k),
                                             as.vector(t(expected_order)))],
                                nrow = nrow(query), byrow = TRUE)
    for (search in c("kdtree", "linear")) {
      nearest <- kNN(x, k, query = query, search = search, bucketSize = bucket,
                     splitRule = rule)
      expect_equal(unname(nearest$id), unname(expected_order))
      expect_equal(unname(nearest$dist), expected_distance, tolerance = 1e-12)
      radius <- frNN(x, eps, query = query, search = search, bucketSize = bucket,
                     splitRule = rule)
      for (i in seq_len(nrow(query))) {
        ids <- which(cross[i, ] <= eps)
        ids <- ids[order(cross[i, ids], ids)]
        expect_equal(radius$id[[i]], unname(ids))
        expect_equal(radius$dist[[i]], unname(cross[i, ids]), tolerance = 1e-12)
      }
    }
  }
})

test_that("tied, duplicate, and boundary neighborhoods agree across backends", {
  skip_on_cran()
  set.seed(20261028)
  for (iteration in seq_len(10)) {
    x <- matrix(sample(-2:2, 48, replace = TRUE), ncol = 2)
    x <- rbind(x, x[rep(1L, 5L), , drop = FALSE])
    distances <- dist(x)
    for (k in c(1L, 4L, nrow(x) - 1L)) {
      expected <- kNN(distances, k)
      for (search in c("kdtree", "linear")) {
        nearest <- kNN(x, k, search = search, bucketSize = 1L)
        # A cutoff tie may select different valid IDs, but the distance
        # multiset, unique IDs, and removal of the self-index must agree.
        expect_equal(nearest$dist, expected$dist, tolerance = 1e-12)
        for (i in seq_len(nrow(x))) {
          expect_false(i %in% nearest$id[i, ])
          expect_equal(anyDuplicated(nearest$id[i, ]), 0L)
          expect_equal(unname(as.matrix(distances)[i, nearest$id[i, ]]),
                       unname(nearest$dist[i, ]), tolerance = 1e-12)
        }
      }
    }
    for (eps in c(0, 1, 2)) {
      expected <- frNN(distances, eps)
      for (search in c("kdtree", "linear")) {
        actual <- frNN(x, eps, search = search, bucketSize = 1L)
        expect_equal(actual$id, expected$id)
        expect_equal(actual$dist, expected$dist, tolerance = 1e-12)
      }
    }
    for (search in c("kdtree", "linear"))
      expect_equal(lof(x, 4, search = search, bucketSize = 1L),
                   lof(distances, 4), tolerance = 1e-10)
  }
})

test_that("ANN resources unwind after native errors and subsequent calls work", {
  skip_on_cran()
  x <- matrix(seq_len(60) / 10, ncol = 2)
  for (search in c("kdtree", "linear")) {
    expected <- kNN(x, 3, search = search)
    for (iteration in seq_len(25)) {
      # The length error is thrown in dbscan_int after allocating points and
      # constructing the tree, exercising native exception unwinding.
      expect_error(dbscan(x, eps = .5, minPts = 3, weights = c(1, 2),
                          search = search), "length of weights")
      expect_equal(kNN(x, 3, search = search), expected)
    }
  }
})
