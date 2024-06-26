
#' @importFrom stats nobs
#' @export
nobs.dbscan <- function(object, ...) length(object$cluster)

#' @export
nobs.hdbscan <- function(object, ...) length(object$cluster)

#' @export
nobs.general_clustering <- function(object, ...) length(object$cluster)

