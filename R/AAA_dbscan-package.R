#' @title `r packageDescription("dbscan")$Package`: `r packageDescription("dbscan")$Title`
#'
#' @description `r packageDescription("dbscan")$Description`
#'
#' @section Key functions:
#' - Clustering: [dbscan()], [hdbscan()], [optics()], [jpclust()], [sNNclust()]
#' - Outliers: [lof()], [glosh()], [pointdensity()]
#' - Nearest Neighbors: [kNN()], [frNN()], [sNN()]
#'
#' @author Michael Hahsler and Matthew Piekenbrock
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
