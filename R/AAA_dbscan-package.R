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
#' @importFrom graphics plot text title
#' @importFrom stats reorder as.dist hclust runif rnorm dist order.dendrogram prcomp
#' @useDynLib dbscan, .registration=TRUE
NULL
