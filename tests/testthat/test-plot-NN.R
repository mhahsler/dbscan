test_that("plot.NN projects data with more than 2 dimensions", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  x <- as.matrix(iris[, -5])
  expect_no_error(plot(kNN(x, 3), x))
  expect_no_error(plot(frNN(x, .5), x))
  expect_no_error(plot(sNN(x, 5), x))
  expect_no_error(plot(kNN(x, 3), as.data.frame(x), xlab = "a", ylab = "b"))
  ## query points
  expect_no_error(plot(kNN(x, 3, query = x[1:2, ]), x))

  ## 2 columns are plotted as is
  expect_no_error(plot(kNN(x[, 1:2], 3), x[, 1:2]))

  ## the plot region covers the principal components (not the first 2 columns)
  plot(kNN(x, 3), x)
  pc <- prcomp(x)$x[, 1:2]
  usr <- par("usr")
  expect_true(usr[1] <= min(pc[, 1]) && usr[2] >= max(pc[, 1]))
  expect_true(usr[3] <= min(pc[, 2]) && usr[4] >= max(pc[, 2]))
  expect_true(usr[1] < min(x[, 1]) - 1)

  expect_error(plot(kNN(x[, 1:2], 3), x[, 1, drop = FALSE]), "at least 2")
})

test_that("projection for plots uses principal component labels", {
  x <- as.matrix(iris[, -5])
  p <- dbscan:::.project_2d(x)
  expect_equal(dim(p$data), c(150L, 2L))
  expect_equal(unname(abs(p$data)), unname(abs(prcomp(x)$x[, 1:2])))
  expect_match(p$xlab, "^PC1 \\([0-9.]+% of variance\\)$")
  expect_match(p$ylab, "^PC2 \\([0-9.]+% of variance\\)$")

  ## user labels are kept; 2 columns are not projected
  p <- dbscan:::.project_2d(x, xlab = "a")
  expect_identical(p$xlab, "a")
  p <- dbscan:::.project_2d(x[, 1:2])
  expect_identical(p$data, x[, 1:2])
  expect_null(p$xlab)
  expect_error(dbscan:::.project_2d(x[, 1, drop = FALSE]), "at least 2")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  cl <- kmeans(x, 3)$cluster
  expect_no_error(hullplot(x, cl))
  expect_no_error(clplot(x, cl, xlab = "a", ylab = "b"))
  expect_no_error(hullplot(x[, 1:2], cl, xlab = "a"))
})
