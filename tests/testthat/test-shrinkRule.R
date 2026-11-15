test_that("shrinkRule gives exact results for all bd-tree shrinking rules", {
  set.seed(20261115)
  x <- matrix(rnorm(200), ncol = 2)
  eps <- .5
  min_pts <- 4

  ref_lof <- lof(x, min_pts, search = "linear")
  ref_density <- pointdensity(x, eps, search = "linear")
  ref_db <- dbscan(x, eps, min_pts, borderPoints = FALSE, search = "linear")
  ref_optics <- optics(x, eps, min_pts, search = "linear")

  for (sr in c("none", "simple", "centroid", "suggest")) {
    expect_equal(lof(x, min_pts, search = "bdtree", shrinkRule = sr), ref_lof,
      tolerance = 1e-10)
    expect_equal(pointdensity(x, eps, search = "bdtree", shrinkRule = sr),
      ref_density)

    db <- dbscan(x, eps, min_pts, borderPoints = FALSE, search = "bdtree",
      shrinkRule = sr)
    expect_identical(db$cluster, ref_db$cluster)

    o <- optics(x, eps, min_pts, search = "bdtree", shrinkRule = sr)
    expect_equal(o$order, ref_optics$order)
    expect_equal(o$reachdist, ref_optics$reachdist, tolerance = 1e-12)

    ## sNN passes the rule on to kNN
    expect_equal(sNN(x, k = 5, search = "bdtree", shrinkRule = sr)$id,
      sNN(x, k = 5, search = "linear")$id)
  }

  for (f in list(
    function() lof(x, min_pts, search = "bdtree", shrinkRule = "bogus"),
    function() pointdensity(x, eps, search = "bdtree", shrinkRule = "bogus"),
    function() dbscan(x, eps, min_pts, search = "bdtree", shrinkRule = "bogus"),
    function() optics(x, eps, min_pts, search = "bdtree", shrinkRule = "bogus")
  ))
    expect_error(f(), "shrinkRule")
})
