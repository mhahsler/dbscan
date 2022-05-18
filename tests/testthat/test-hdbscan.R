library("dbscan")
library("testthat")

context("HDBSCAN")


data("iris")

## minPts not given
expect_error(dbscan::hdbscan(iris))

## Expects numerical data; species is factor
expect_error(dbscan::hdbscan(iris, minPts = 4))

iris <- as.matrix(iris[,1:4])

res <- dbscan::hdbscan(iris, minPts = 4)
expect_identical(length(res$cluster), nrow(iris))

## expected result of table(res$cluster) is:
expect_equivalent(table(res$cluster),
                  as.table(c("1" = 100L, "2" = 50L)))

## compare on moons data
data("moons")
res <- dbscan::hdbscan(moons, minPts = 5)
expect_identical(length(res$cluster), nrow(moons))

## Check hierarchy matches dbscan* at every value
check <- rep(F, nrow(moons)-1)
core_dist <- kNNdist(moons, k=5-1)

## cutree doesn't distinguish noise as 0, so we make a new method to do it manually
cut_tree <- function(hcl, eps, core_dist){
  cuts <- unname(cutree(hcl, h=eps))
  cuts[which(core_dist > eps)] <- 0 # Use core distance to distinguish noise
  cuts
}

eps_values <- sort(res$hc$height, decreasing = T)+.Machine$double.eps ## Machine eps for consistency between cuts
for (i in 1:length(eps_values)) {
  cut_cl <- cut_tree(res$hc, eps_values[i], core_dist)
  dbscan_cl <- dbscan(moons, eps = eps_values[i], minPts = 5, borderPoints = F) # DBSCAN* doesn't include border points

  ## Use run length encoding as an ID-independent way to check ordering
  check[i] <- (all.equal(rle(cut_cl)$lengths, rle(dbscan_cl$cluster)$lengths) == "TRUE")
}

expect_true(all(check))

## Expect generating extra trees doesn't fail
res <- dbscan::hdbscan(moons, minPts = 5, gen_hdbscan_tree = TRUE, gen_simplified_tree = TRUE)
expect_s3_class(res, "hdbscan")

## Expect hdbscan tree matches stats:::as.dendrogram version of hclust object
hc_dend <- as.dendrogram(res$hc)
expect_s3_class(hc_dend, "dendrogram")
expect_equal(hc_dend, res$hdbscan_tree)

## Expect hdbscan works with non-euclidean distances
dist_moons <- dist(moons, method = "canberra")
res <- dbscan::hdbscan(dist_moons, minPts = 5)
expect_s3_class(res, "hdbscan")


## test mrdist()
context("mrdist()")
mrdist(cbind(1:10), 2)
mrdist(cbind(1:11), 2)
mrdist(dist(cbind(1:10)), 2)
mrdist(dist(cbind(1:11)), 2)

