#' Clusterable Data (HDBSCAN Soft Clustering Example)
#'
#' Contains 2309 2-d points with clusters of different shapes, sizes and
#' densities plus background noise.
#' The data is from McInnes, Healy and Astels (2017) and is available as part
#' of the contributed scikit-learn compatible project HDBSCAN.
#'
#' @name clusterable_data
#' @docType data
#' @format A data frame with 2309 observations on the following 2 variables.
#' \describe{
#' \item{X}{a numeric vector}
#' \item{Y}{a numeric vector} }
#' @source Soft clustering explanation in the HDBSCAN documentation:
#' \url{https://hdbscan.readthedocs.io/en/latest/soft_clustering_explanation.html}.
#' Data file:
#' \url{https://github.com/scikit-learn-contrib/hdbscan/blob/master/notebooks/clusterable_data.npy}
#' @references McInnes, Leland, John Healy, and Steve Astels (2017).
#' hdbscan: Hierarchical density based clustering. _Journal of Open Source
#' Software_ 2(11):205. \doi{10.21105/joss.00205}
#' @keywords datasets
#' @examples
#' data(clusterable_data)
#' plot(clusterable_data)
#'
#' cl <- hdbscan(clusterable_data, minPts = 15)
#' cl
#'
#' plot(cl)
#' clplot(clusterable_data, cl)
NULL
