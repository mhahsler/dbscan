#' @title `r packageDescription("dbscan")$Package`: `r packageDescription("dbscan")$Title`
#'
#' @description `r sub("<[^>]+>", "", packageDescription("dbscan")$Description)`
#'
#' @section Key functions:
#' - Clustering: [dbscan()], [hdbscan()], [optics()], [jpclust()], [sNNclust()]
#' - Outliers: [lof()], [glosh()], [pointdensity()]
#' - Nearest Neighbors: [kNN()], [frNN()], [sNN()]
#'
#' @author Michael Hahsler and Matthew Piekenbrock
#'
#' @references
#' Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based Clustering with R. Journal of Statistical Software, 91(1), 1-30. \doi{10.18637/jss.v091.i01}
#'
#' @docType package
#' @name dbscan-package
#'
#' @import Rcpp
#' @importFrom graphics plot points lines text abline polygon par segments
#' @importFrom grDevices palette chull adjustcolor
#' @importFrom stats dist hclust dendrapply as.dendrogram is.leaf prcomp
#' @importFrom utils tail
#'
#' @useDynLib dbscan, .registration=TRUE
NULL
