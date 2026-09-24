test_that("integer parameters are scalar, finite, and integer-valued", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
  invalid <- list(
    numeric(),
    c(1, 2),
    NA_real_,
    NaN,
    Inf,
    -Inf,
    1.9,
    0,
    -1,
    "1",
    TRUE
  )

  for (value in invalid) {
    expect_error(kNN(x, k = value), "k must be")
    expect_error(dbscan(x, eps = 1, minPts = value), "minPts must be")
  }

  expect_s3_class(kNN(x, k = 1), "kNN")
  expect_s3_class(dbscan(x, eps = 0, minPts = 1), "dbscan")
  expect_error(optics(x, eps = 1, minPts = 1.5), "minPts must be")
  expect_error(lof(x, minPts = 1.5), "minPts must be")
  expect_error(hdbscan(x, minPts = 1.5), "minPts must be")
  expect_error(sNN(kNN(x, 2), k = 1.5), "k must be")
})

test_that("eps and approx are finite nonnegative scalars", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
  invalid <- list(
    numeric(),
    c(0, 1),
    NA_real_,
    NaN,
    Inf,
    -Inf,
    -1,
    "0",
    FALSE
  )

  for (value in invalid) {
    expect_error(frNN(x, eps = value), "eps must be")
    expect_error(kNN(x, k = 1, approx = value), "approx must be")
  }

  expect_s3_class(frNN(x, eps = 0), "frNN")
  expect_s3_class(kNN(x, k = 1, approx = 0), "kNN")
  expect_error(dbscan(x, eps = Inf, minPts = 2), "eps must be")
  expect_s3_class(optics(x, eps = Inf, minPts = 2), "optics")
  expect_true(is.finite(optics(x, minPts = 2)$eps))

  for (value in list(numeric(), c(0, 1), NA_real_, NaN, -Inf, -1, "0", FALSE)) {
    expect_error(optics(x, eps = value, minPts = 2), "eps must be")
  }

})

test_that("bucketSize is a positive integer scalar", {
  x <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  invalid <- list(
    numeric(),
    c(1, 2),
    NA_real_,
    NaN,
    Inf,
    -Inf,
    0,
    -1,
    1.5,
    "1",
    TRUE
  )

  for (value in invalid)
    expect_error(kNN(x, k = 1, bucketSize = value), "bucketSize must be")

  expect_s3_class(kNN(x, k = 1, bucketSize = 1), "kNN")
})

test_that("all direct neighbor APIs validate search controls", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)

  expect_error(frNN(x, 1, bucketSize = 0), "bucketSize must be")
  expect_error(dbscan(x, 1, minPts = 2, bucketSize = 0), "bucketSize must be")
  expect_error(optics(x, 1, minPts = 2, bucketSize = 0), "bucketSize must be")
  expect_error(lof(x, minPts = 2, bucketSize = 0), "bucketSize must be")
  expect_error(pointdensity(x, 1, bucketSize = 0), "bucketSize must be")

  expect_error(frNN(x, 1, approx = Inf), "approx must be")
  expect_error(dbscan(x, 1, minPts = 2, approx = Inf), "approx must be")
  expect_error(optics(x, 1, minPts = 2, approx = Inf), "approx must be")
  expect_error(lof(x, minPts = 2, approx = Inf), "approx must be")
  expect_error(pointdensity(x, 1, approx = Inf), "approx must be")
})

test_that("matrix inputs are numeric, nonempty, and finite", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)

  expect_error(kNN(matrix(numeric(), nrow = 0, ncol = 2), 1), "0 rows")
  expect_error(kNN(matrix(numeric(), nrow = 2, ncol = 0), 1), "0 columns")
  expect_error(kNN(matrix(letters[1:10], ncol = 2), 1), "numeric matrix")

  for (value in list(NA_real_, NaN, Inf, -Inf)) {
    bad <- x
    bad[1, 1] <- value
    expect_error(kNN(bad, 1), "NA, NaN, or infinite")
  }

  bad_query <- matrix(c(NA_real_, 0), nrow = 1)
  expect_error(kNN(x, 1, query = bad_query), "query cannot contain")

  x_inf <- x
  x_inf[1, 1] <- Inf
  expect_error(frNN(x_inf, 1), "infinite")
  expect_error(dbscan(x_inf, 1, minPts = 2), "infinite")
  expect_error(optics(x_inf, 1, minPts = 2), "infinite")
  expect_error(lof(x_inf, minPts = 2), "infinite")
  expect_error(hdbscan(x_inf, minPts = 2), "infinite")
  expect_error(pointdensity(x_inf, 1), "infinite")
})
