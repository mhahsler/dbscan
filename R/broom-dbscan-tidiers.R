#' Turn a dbscan clustering object into a tidy tibble
#'
#' Provides [tidy()][generics::tidy()], [augment()][generics::augment()], and
#' [glance()][generics::glance()] verbs for clusterings created with algorithms
#' in package `dbscan` to work with [tidymodels](https://www.tidymodels.org/).
#'
#' Methods are available for the clustering objects created by [dbscan()],
#' [hdbscan()], [jpclust()] and [sNNclust()]. Clusterings created by [optics()]
#' are not supported.
#'
#' @details
#' All three verbs return a [tibble][tibble::tibble()]. Cluster 0 represents
#' noise points. Noise points are not assigned to a cluster.
#'
#' ## `tidy()`
#'
#' Returns one row per cluster, including the noise cluster 0, with the columns:
#'
#' * `cluster`: the cluster label as a factor (the levels are `0` to the number
#'   of clusters).
#' * `size`: the number of points in the cluster.
#' * `noise`: `TRUE` for cluster 0 (the noise points).
#' * `cluster_score` (only for [hdbscan()]): the stability score of the cluster
#'   (see the component `cluster_scores` in [hdbscan()]). It is `NA` for the
#'   noise cluster.
#'
#' ## `augment()`
#'
#' Returns the data (or `newdata`) as a tibble with additional columns. The
#' original data used to create the clustering always needs to be passed as
#' `data` because the clustering object does not store it.
#' If `newdata` is specified, then the cluster labels of the new points are
#' predicted with `predict()`. Otherwise, the labels in the clustering are
#' used. The added columns are:
#'
#' * `.cluster`: the cluster label as a factor.
#' * `noise` (only for [dbscan()]): `TRUE` if the point is a noise point.
#' * `.coredist`, `.membership_prob` and `.outlier_scores` (only for
#'   [hdbscan()]): the core distance, the cluster membership probability and
#'   the GLOSH outlier score of each point (see [hdbscan()]). They are `NA` for
#'   points in `newdata`.
#'
#' Augmenting `newdata` is not supported for [jpclust()] and [sNNclust()].
#'
#' ## `glance()`
#'
#' Returns a one-row summary of the clustering with the columns:
#'
#' * `nobs`: the number of points.
#' * `n.clusters`: the number of clusters (without noise).
#' * `nexcluded`: the number of noise points.
#'
#' @param x A clustering object returned by [dbscan()], [hdbscan()],
#'   [jpclust()] or [sNNclust()].
#' @param data The data used to create the clustering.
#' @param newdata New data to predict cluster labels for.
#' @param ... further arguments are ignored without a warning.
#'
#' @returns A [tibble][tibble::tibble()] (see Details for the columns).
#'
#' @name dbscan_tidiers
#' @aliases dbscan_tidiers glance tidy augment
#' @family tidiers
#'
#' @seealso [generics::tidy()], [generics::augment()],
#'  [generics::glance()], [dbscan()]
#'
#' @examplesIf requireNamespace("tibble", quietly = TRUE) && identical(Sys.getenv("NOT_CRAN"), "true")
#'
#' data(iris)
#' x <- scale(iris[, 1:4])
#'
#' ## dbscan
#' db <- dbscan(x, eps = .9, minPts = 5)
#' db
#'
#' # summarize model fit with tidiers
#' tidy(db)
#' glance(db)
#'
#' # augment for this model needs the original data
#' augment(db, x)
#'
#' # to augment new data, the original data is also needed
#' augment(db, x, newdata = x[1:5, ])
#'
#' ## hdbscan
#' hdb <- hdbscan(x, minPts = 5)
#'
#' # summarize model fit with tidiers
#' tidy(hdb)
#' glance(hdb)
#'
#' # augment for this model needs the original data
#' augment(hdb, x)
#'
#' # to augment new data, the original data is also needed
#' augment(hdb, x, newdata = x[1:5, ])
#'
#' ## Jarvis-Patrick clustering
#' cl <- jpclust(x, k = 20, kt = 15)
#'
#' # summarize model fit with tidiers
#' tidy(cl)
#' glance(cl)
#'
#' # augment for this model needs the original data
#' augment(cl, x)
#'
#' ## Shared Nearest Neighbor clustering
#' cl <- sNNclust(x, k = 20, eps = 0.8, minPts = 15)
#'
#' # summarize model fit with tidiers
#' tidy(cl)
#' glance(cl)
#'
#' # augment for this model needs the original data
#' augment(cl, x)
#'
NULL

#' @rdname dbscan_tidiers
#' @importFrom generics tidy
#' @export
generics::tidy


#' @rdname dbscan_tidiers
#' @export
tidy.dbscan <- function(x, ...) {
  n_cl <- max(x$cluster)
  size <- table(factor(x$cluster, levels = 0:n_cl))

  tb <- tibble::tibble(cluster = as.factor(0:n_cl),
         size = as.integer(size))

  tb$noise <- tb$cluster == 0L
  tb
}

#' @rdname dbscan_tidiers
#' @export
tidy.hdbscan <- function(x, ...) {
  n_cl <- max(x$cluster)
  size <- table(factor(x$cluster, levels = 0:n_cl))

  tb <- tibble::tibble(cluster = as.factor(0:n_cl),
         size = as.integer(size))
  tb$cluster_score <- as.numeric(x$cluster_scores[as.character(tb$cluster)])
  tb$noise <- tb$cluster == 0L

  tb
}

#' @rdname dbscan_tidiers
#' @export
tidy.general_clustering <- function(x, ...) {
  n_cl <- max(x$cluster)
  size <- table(factor(x$cluster, levels = 0:n_cl))

  tb <- tibble::tibble(cluster = as.factor(0:n_cl),
         size = as.integer(size))
  tb$noise <- tb$cluster == 0L

  tb
}


## augment

#' @importFrom generics augment
#' @rdname dbscan_tidiers
#' @export
generics::augment


#' @rdname dbscan_tidiers
#' @export
augment.dbscan <- function(x, data = NULL, newdata = NULL, ...) {
  n_cl <- max(x$cluster)

  if (is.null(data) && is.null(newdata))
    stop("Must specify either `data` or `newdata` argument.")

  if (is.null(data) || nrow(data) != length(x$cluster)) {
    stop("The original data needs to be passed as data.")
  }

  if (is.null(newdata)) {
    tb <- tibble::as_tibble(data)
    tb$.cluster <- factor(x$cluster, levels = 0:n_cl)
  } else {
    tb <- tibble::as_tibble(newdata)
    tb$.cluster <- factor(predict(x,
                                  newdata = newdata,
                                  data = data), levels = 0:n_cl)
  }

  tb$noise <- tb$.cluster == 0L

  tb
}

#' @rdname dbscan_tidiers
#' @export
augment.hdbscan <- function(x, data = NULL, newdata = NULL, ...) {
  n_cl <- max(x$cluster)

  if (is.null(data) || nrow(data) != length(x$cluster)) {
    stop("The original data needs to be passed as data.")
  }

  if (is.null(newdata)) {
    tb <- tibble::as_tibble(data)
    tb$.cluster <- factor(x$cluster, levels = 0:n_cl)
    tb$.coredist <- x$coredist
    tb$.membership_prob <- x$membership_prob
    tb$.outlier_scores <- x$outlier_scores
  } else {
    tb <- tibble::as_tibble(newdata)
    tb$.cluster <- factor(
        predict(x, newdata = newdata, data = data), levels = 0:n_cl)
    tb$.coredist <- NA_real_
    tb$.membership_prob <- NA_real_
    tb$.outlier_scores <- NA_real_
  }

  tb
}

#' @rdname dbscan_tidiers
#' @export
augment.general_clustering <- function(x, data = NULL, newdata = NULL, ...) {
  n_cl <- max(x$cluster)

  if (is.null(data) || nrow(data) != length(x$cluster)) {
    stop("The original data needs to be passed as data.")
  }

  if (is.null(newdata)) {
    tb <- tibble::as_tibble(data)
    tb$.cluster <- factor(x$cluster, levels = 0:n_cl)
  } else {
    stop("augmenting new data is not supported.")
  }

  tb
}



## glance
#' @importFrom generics glance
#' @rdname dbscan_tidiers
#' @export
generics::glance


#' @rdname dbscan_tidiers
#' @export
glance.dbscan <- function(x, ...) {
  tibble::tibble(
    nobs = length(x$cluster),
    n.clusters = length(table(x$cluster[x$cluster != 0L])),
    nexcluded = sum(x$cluster == 0L)
  )
}

#' @rdname dbscan_tidiers
#' @export
glance.hdbscan <- function(x, ...) {
  tibble::tibble(
    nobs = length(x$cluster),
    n.clusters = length(table(x$cluster[x$cluster != 0L])),
    nexcluded = sum(x$cluster == 0L)
  )
}

#' @rdname dbscan_tidiers
#' @export
glance.general_clustering <- function(x, ...) {
  tibble::tibble(
    nobs = length(x$cluster),
    n.clusters = length(table(x$cluster[x$cluster != 0L])),
    nexcluded = sum(x$cluster == 0L)
  )
}

