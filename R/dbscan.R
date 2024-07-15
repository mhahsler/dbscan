#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#' Density-based Spatial Clustering of Applications with Noise (DBSCAN)
#'
#' Fast reimplementation of the DBSCAN (Density-based spatial clustering of
#' applications with noise) clustering algorithm using a kd-tree.
#'
#' The
#' implementation is significantly faster and can work with larger data sets
#' than [fpc::dbscan()] in \pkg{fpc}. Use `dbscan::dbscan()` (with specifying the package) to
#' call this implementation when you also load package \pkg{fpc}.
#'
#' **The algorithm**
#'
#' This implementation of DBSCAN follows the original
#' algorithm as described by Ester et al (1996). DBSCAN performs the following steps:
#'
#' 1. Estimate the density
#'   around each data point by counting the number of points in a user-specified
#'   eps-neighborhood and applies a used-specified minPts thresholds to identify
#'      - core points (points with more than minPts points in their neighborhood),
#'      - border points (non-core points with a core point in their neighborhood) and
#'      - noise points (all other points).
#' 2. Core points form the backbone of clusters by joining them into
#'   a cluster if they are density-reachable from each other (i.e., there is a chain of core
#'   points where one falls inside the eps-neighborhood of the next).
#' 3. Border points are assigned to clusters. The algorithm needs parameters
#'   `eps` (the radius of the epsilon neighborhood) and `minPts` (the
#'   density threshold).
#'
#' Border points are arbitrarily assigned to clusters in the original
#' algorithm. DBSCAN* (see Campello et al 2013) treats all border points as
#' noise points. This is implemented with `borderPoints = FALSE`.
#'
#' **Specifying the data**
#'
#' If `x` is a matrix or a data.frame, then fast fixed-radius nearest
#' neighbor computation using a kd-tree is performed using Euclidean distance.
#' See [frNN()] for more information on the parameters related to
#' nearest neighbor search. **Note** that only numerical values are allowed in `x`.
#'
#' Any precomputed distance matrix (dist object) can be specified as `x`.
#' You may run into memory issues since distance matrices are large.
#'
#' A precomputed frNN object can be supplied as `x`. In this case
#' `eps` does not need to be specified. This option us useful for large
#' data sets, where a sparse distance matrix is available. See
#' [frNN()] how to create frNN objects.
#'
#' **Setting parameters for DBSCAN**
#'
#' The parameters `minPts` and `eps` define the minimum density required
#' in the area around core points which form the backbone of clusters.
#' `minPts` is the number of points
#' required in the neighborhood around the point defined by the parameter `eps`
#' (i.e., the radius around the point). Both parameters
#' depend on each other and changing one typically requires changing
#' the other one as well. The parameters also depend on the size of the data set with
#' larger datasets requiring a larger `minPts` or a smaller `eps`.
#'
#' * `minPts:` The original
#' DBSCAN paper (Ester et al, 1996) suggests to start by setting \eqn{\text{minPts} \ge d + 1},
#' the data dimensionality plus one or higher with a minimum of 3. Larger values
#' are preferable since increasing the parameter suppresses more noise in the data
#' by requiring more points to form clusters.
#' Sander et al (1998) uses in the examples two times the data dimensionality.
#' Note that setting \eqn{\text{minPts} \le 2} is equivalent to hierarchical clustering
#' with the single link metric and the dendrogram cut at height `eps`.
#'
#' * `eps:` A suitable neighborhood size
#' parameter `eps` given a fixed value for `minPts` can be found
#' visually by inspecting the [kNNdistplot()] of the data using
#' \eqn{k = \text{minPts} - 1} (`minPts` includes the point itself, while the
#' k-nearest neighbors distance does not). The k-nearest neighbor distance plot
#' sorts all data points by their k-nearest neighbor distance. A sudden
#' increase of the kNN distance (a knee) indicates that the points to the right
#' are most likely outliers. Choose `eps` for DBSCAN where the knee is.
#'
#' **Predict cluster memberships**
#'
#' [predict()] can be used to predict cluster memberships for new data
#' points. A point is considered a member of a cluster if it is within the eps
#' neighborhood of a core point of the cluster. Points
#' which cannot be assigned to a cluster will be reported as
#' noise points (i.e., cluster ID 0).
#' **Important note:** `predict()` currently can only use Euclidean distance to determine
#' the neighborhood of core points. If `dbscan()` was called using distances other than Euclidean,
#' then the neighborhood calculation will not be correct and only approximated by Euclidean
#' distances. If the data contain factor columns (e.g., using Gower's distance), then
#' the factors in `data` and `query` first need to be converted to numeric to use the
#' Euclidean approximation.
#'
#'
#' @aliases dbscan DBSCAN print.dbscan_fast
#' @family clustering functions
#'
#' @param x a data matrix, a data.frame, a [dist] object or a [frNN] object with
#' fixed-radius nearest neighbors.
#' @param eps size (radius) of the epsilon neighborhood. Can be omitted if
#' `x` is a frNN object.
#' @param minPts number of minimum points required in the eps neighborhood for
#' core points (including the point itself).
#' @param weights numeric; weights for the data points. Only needed to perform
#' weighted clustering.
#' @param borderPoints logical; should border points be assigned to clusters.
#' The default is `TRUE` for regular DBSCAN. If `FALSE` then border
#' points are considered noise (see DBSCAN* in Campello et al, 2013).
#' @param ...  additional arguments are passed on to the fixed-radius nearest
#' neighbor search algorithm. See [frNN()] for details on how to
#' control the search strategy.
#'
#' @return `dbscan()` returns an object of class `dbscan_fast` with the following components:
#'
#' \item{eps }{ value of the `eps` parameter.}
#' \item{minPts }{ value of the `minPts` parameter.}
#' \item{metric }{ used distance metric.}
#' \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points.}
#'
#' `is.corepoint()` returns a logical vector indicating for each data point if it is a
#'   core point.
#'
#' @author Michael Hahsler
#' @references Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast
#' Density-Based Clustering with R.  _Journal of Statistical Software,_
#' 91(1), 1-30.
#' \doi{10.18637/jss.v091.i01}
#'
#' Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu (1996). A
#' Density-Based Algorithm for Discovering Clusters in Large Spatial Databases
#' with Noise. Institute for Computer Science, University of Munich.
#' _Proceedings of 2nd International Conference on Knowledge Discovery and
#' Data Mining (KDD-96),_ 226-231.
#' \url{https://dl.acm.org/doi/10.5555/3001460.3001507}
#'
#' Campello, R. J. G. B.; Moulavi, D.; Sander, J. (2013). Density-Based
#' Clustering Based on Hierarchical Density Estimates. Proceedings of the
#' 17th Pacific-Asia Conference on Knowledge Discovery in Databases, PAKDD
#' 2013, _Lecture Notes in Computer Science_ 7819, p. 160.
#' \doi{10.1007/978-3-642-37456-2_14}
#'
#' Sander, J., Ester, M., Kriegel, HP. et al. (1998). Density-Based
#' Clustering in Spatial Databases: The Algorithm GDBSCAN and Its Applications.
#' _Data Mining and Knowledge Discovery_ 2, 169-194.
#' \doi{10.1023/A:1009745219419}
#'
#' @keywords model clustering
#' @examples
#' ## Example 1: use dbscan on the iris data set
#' data(iris)
#' iris <- as.matrix(iris[, 1:4])
#'
#' ## Find suitable DBSCAN parameters:
#' ## 1. We use minPts = dim + 1 = 5 for iris. A larger value can also be used.
#' ## 2. We inspect the k-NN distance plot for k = minPts - 1 = 4
#' kNNdistplot(iris, minPts = 5)
#'
#' ## Noise seems to start around a 4-NN distance of .7
#' abline(h=.7, col = "red", lty = 2)
#'
#' ## Cluster with the chosen parameters
#' res <- dbscan(iris, eps = .7, minPts = 5)
#' res
#'
#' pairs(iris, col = res$cluster + 1L)
#'
#' ## Use a precomputed frNN object
#' fr <- frNN(iris, eps = .7)
#' dbscan(fr, minPts = 5)
#'
#' ## Example 2: use data from fpc
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
#'   y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
#'   )
#'
#' res <- dbscan(x, eps = .3, minPts = 3)
#' res
#'
#' ## plot clusters and add noise (cluster 0) as crosses.
#' plot(x, col = res$cluster)
#' points(x[res$cluster == 0, ], pch = 3, col = "grey")
#'
#' hullplot(x, res)
#'
#' ## Predict cluster membership for new data points
#' ## (Note: 0 means it is predicted as noise)
#' newdata <- x[1:5,] + rnorm(10, 0, .3)
#' hullplot(x, res)
#' points(newdata, pch = 3 , col = "red", lwd = 3)
#' text(newdata, pos = 1)
#'
#' pred_label <- predict(res, newdata, data = x)
#' pred_label
#' points(newdata, col = pred_label + 1L,  cex = 2, lwd = 2)
#'
#' ## Compare speed against fpc version (if microbenchmark is installed)
#' ## Note: we use dbscan::dbscan to make sure that we do now run the
#' ## implementation in fpc.
#' \dontrun{
#' if (requireNamespace("fpc", quietly = TRUE) &&
#'     requireNamespace("microbenchmark", quietly = TRUE)) {
#'   t_dbscan <- microbenchmark::microbenchmark(
#'     dbscan::dbscan(x, .3, 3), times = 10, unit = "ms")
#'   t_dbscan_linear <- microbenchmark::microbenchmark(
#'     dbscan::dbscan(x, .3, 3, search = "linear"), times = 10, unit = "ms")
#'   t_dbscan_dist <- microbenchmark::microbenchmark(
#'     dbscan::dbscan(x, .3, 3, search = "dist"), times = 10, unit = "ms")
#'   t_fpc <- microbenchmark::microbenchmark(
#'     fpc::dbscan(x, .3, 3), times = 10, unit = "ms")
#'
#'   r <- rbind(t_fpc, t_dbscan_dist, t_dbscan_linear, t_dbscan)
#'   r
#'
#'   boxplot(r,
#'     names = c('fpc', 'dbscan (dist)', 'dbscan (linear)', 'dbscan (kdtree)'),
#'     main = "Runtime comparison in ms")
#'
#'   ## speedup of the kd-tree-based version compared to the fpc implementation
#'   median(t_fpc$time) / median(t_dbscan$time)
#' }}
#'
#' ## Example 3: manually create a frNN object for dbscan (dbscan only needs ids and eps)
#' nn <- structure(list(ids = list(c(2,3), c(1,3), c(1,2,3), c(3,5), c(4,5)), eps = 1),
#'   class =  c("NN", "frNN"))
#' nn
#' dbscan(nn, minPts = 2)
#'
#' @export
dbscan <-
  function(x,
    eps,
    minPts = 5,
    weights = NULL,
    borderPoints = TRUE,
    ...) {
    if (inherits(x, "frNN") && missing(eps))
      eps <- x$eps

    if (inherits(x, "dist")) {
      .check_dist(x)
      dist_method <- attr(x, "method")
    } else
      dist_method <- "euclidean"

    if (is.null(dist_method))
      dist_method <- "unknown"

    ### extra contains settings for frNN
    ### search = "kdtree", bucketSize = 10, splitRule = "suggest", approx = 0
    ### also check for MinPts for fpc compatibility (does not work for
    ### search method dist)
    extra <- list(...)
    args <-
      c("MinPts", "search", "bucketSize", "splitRule", "approx")
    m <- pmatch(names(extra), args)
    if (anyNA(m))
      stop("Unknown parameter: ",
        toString(names(extra)[is.na(m)]))
    names(extra) <- args[m]

    if (!is.null(extra$MinPts)) {
      warning("converting argument MinPts (fpc) to minPts (dbscan)!")
      minPts <- extra$MinPts
    }
    extra$MinPts <- NULL

    search <- extra$search %||% "kdtree"
    splitRule <- extra$splitRule %||% "suggest"
    search <- .parse_search(search)
    splitRule <- .parse_splitRule(splitRule)

    bucketSize <- if (is.null(extra$bucketSize))
      10L
    else
      as.integer(extra$bucketSize)

    approx <-
      if (is.null(extra$approx))
        0L
    else
      as.integer(extra$approx)

    ### do dist search
    if (search == 3) {
      if (!inherits(x, "dist"))
        if (.matrixlike(x))
          x <- dist(x)
      else
        stop("x needs to be a matrix to calculate distances")
    }

    ## for dist we provide the R code with a frNN list and no x
    frNN <- list()
    if (inherits(x, "dist")) {
      frNN <- frNN(x, eps, ...)$id
      x <- matrix(0.0, nrow = 0, ncol = 0)
    } else if (inherits(x, "frNN")) {
      if (x$eps != eps) {
        eps <- x$eps
        warning("Using the eps of ",
          eps,
          " provided in the fixed-radius NN object.")
      }
      frNN <- x$id
      x <- matrix(0.0, nrow = 0, ncol = 0)

    } else{
      if (!.matrixlike(x))
        stop("x needs to be a matrix or data.frame.")
      ## make sure x is numeric
      x <- as.matrix(x)
      if (storage.mode(x) == "integer")
        storage.mode(x) <- "double"
      if (storage.mode(x) != "double")
        stop("all data in x has to be numeric.")
    }

    if (length(frNN) == 0 && anyNA(x))
      stop("data/distances cannot contain NAs for dbscan (with kd-tree)!")

    ## add self match and use C numbering if frNN is used
    if (length(frNN) > 0)
      frNN <-
      lapply(
        seq_along(frNN),
        FUN = function(i)
          c(i - 1L, frNN[[i]] - 1L)
      )

    if (length(minPts) != 1 ||
        !is.finite(minPts) ||
        minPts < 0)
      stop("minPts need to be a single integer >=0.")

    if (is.null(eps) ||
        is.na(eps) || eps < 0)
      stop("eps needs to be >=0.")

    ret <- dbscan_int(
      x,
      as.double(eps),
      as.integer(minPts),
      as.double(weights),
      as.integer(borderPoints),
      as.integer(search),
      as.integer(bucketSize),
      as.integer(splitRule),
      as.double(approx),
      frNN
    )

    structure(
      list(
        cluster = ret,
        eps = eps,
        minPts = minPts,
        metric = dist_method,
        borderPoints = borderPoints
      ),
      class = c("dbscan_fast", "dbscan")
    )
  }

#' @export
print.dbscan_fast <- function(x, ...) {
  cl <- unique(x$cluster)
  cl <- length(cl[cl != 0L])

  writeLines(c(
    paste0("DBSCAN clustering for ", length(x$cluster), " objects."),
    paste0("Parameters: eps = ", x$eps, ", minPts = ", x$minPts),
    paste0(
      "Using ",
      x$metric,
      " distances and borderpoints = ",
      x$borderPoints
    ),
    paste0(
      "The clustering contains ",
      cl,
      " cluster(s) and ",
      sum(x$cluster == 0L),
      " noise points."
    )
  ))

  print(table(x$cluster))
  cat("\n")

  writeLines(strwrap(paste0(
    "Available fields: ",
    toString(names(x))
  ), exdent = 18))
}

#' @rdname dbscan
#' @export
is.corepoint <- function(x, eps, minPts = 5, ...)
  lengths(frNN(x, eps = eps, ...)$id) >= (minPts - 1)
