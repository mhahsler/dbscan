#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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

#' Ordering Points to Identify the Clustering Structure (OPTICS)
#'
#' Implementation of the OPTICS (Ordering points to identify the clustering
#' structure) point ordering algorithm using a kd-tree.
#'
#' ## The Algorithm
#'
#' This implementation of OPTICS implements the original
#' algorithm as described by Ankerst et al (1999). OPTICS is an ordering
#' algorithm with methods to extract a clustering from the ordering.
#' While using similar concepts as DBSCAN, `minPts` in OPTICS has a different
#' effect than in DBSCAN. Since it is also used to calculate the reachability
#' distance, larger values will make the reachability distance plot smoother.
#' The parameter `eps` is optional and defaults to `Inf`.
#' It represents
#' an upper limit for the neighborhood size used to reduce
#' computational complexity which is helpful for large data sets.
#'
#' OPTICS linearly orders the data points such that points which are spatially
#' closest become neighbors in the ordering. The closest analog to this
#' ordering is dendrogram in single-link hierarchical clustering. The algorithm
#' also calculates the reachability distance for each point.
#' `plot()` (see [reachability_plot])
#' produces a reachability plot which shows each points reachability distance
#' between two consecutive points
#' where the points are sorted by OPTICS. Valleys represent clusters (the
#' deeper the valley, the more dense the cluster) and high points indicate
#' points between clusters.
#'
#' ## Specifying the Data
#'
#' If `x` is specified as a data matrix, then Euclidean distances and fast
#' nearest neighbor lookup using a kd-tree are used. See [kNN()] for
#' details on the parameters for the kd-tree.
#'
#' ## Extracting a Clustering
#'
#' Several methods to extract a clustering from the order returned by OPTICS are
#' implemented:
#'
#' * `extractDBSCAN()` extracts a clustering from an OPTICS ordering that is
#'   similar to what DBSCAN would produce with an eps set to `eps_cl` (see
#'   Ankerst et al, 1999). The only difference to a DBSCAN clustering is that
#'   OPTICS is not able to assign some border points and reports them instead as
#'   noise.
#'
#' * `extractXi()` extract clusters hierarchically specified in Ankerst et al
#'   (1999) based on the steepness of the reachability plot. One interpretation
#'   of the `xi` parameter is that it classifies clusters by change in
#'   relative cluster density. The used algorithm was originally contributed by
#'   the ELKI framework and is explained in Schubert et al (2018), but contains a
#'   set of fixes.
#'
#' ## Predict Cluster Memberships
#'
#' `predict()` requires an extracted DBSCAN clustering with `extractDBSCAN()` and then
#' uses predict for `dbscan()`.
#'
#' @aliases optics OPTICS
#' @family clustering functions
#'
#' @param x a data matrix or a [dist] object.
#' @param eps maximum epsilon neighborhood size only used for performance.
#' When set to `Inf`, then the actual maximal needed
#' radius is estimated from `minPts` and the data.
#' The upper limit can be further reduced to improve
#' performance. If set
#' too low then many reachability values will erroneously become `Inf`
#' shown as dashed lines in the reachability plot. `eps` should be increased.
#' @param minPts the parameter is used to identify dense neighborhoods and the
#' reachability distance is calculated as the distance to the minPts nearest
#' neighbor. Controls the smoothness of the reachability distribution. Default
#' is 5 points.
#' @param eps_cl Threshold to identify clusters (`eps_cl <= eps`).
#' @param xi Steepness threshold to identify clusters hierarchically using the
#' Xi method.
#' @param object an object of class `optics`.
#' @param minimum logical, representing whether or not to extract the minimal
#' (non-overlapping) clusters in the Xi clustering algorithm.
#' @param correctPredecessors logical, correct a common artifact by pruning
#' the steep up area for points that have predecessors not in the
#' cluster--found by the ELKI framework, see details below.
#' @param ...  additional arguments are passed on to fixed-radius nearest
#' neighbor search algorithm. See [frNN()] for details on how to
#' control the search strategy.
#' @param cluster,predecessor plot clusters and predecessors.
#'
#' @return An object of class `optics` with components:
#' \item{eps }{ value of `eps` parameter. }
#' \item{minPts }{ value of `minPts` parameter. }
#' \item{order }{ optics order for the data points in `x`. }
#' \item{reachdist }{ [reachability] distance for each data point in `x`. }
#' \item{coredist }{ core distance for each data point in `x`. }
#'
#' For `extractDBSCAN()`, in addition the following
#' components are available:
#' \item{eps_cl }{ the value of the `eps_cl` parameter. }
#' \item{cluster }{ assigned cluster labels in the order of the data points in `x`. }
#'
#' For `extractXi()`, in addition the following components
#' are available:
#' \item{xi}{ Steepness threshold `xi`. }
#' \item{cluster }{ assigned cluster labels in the order of the data points in `x`.}
#' \item{clusters_xi }{ data.frame containing the start and end of each cluster
#' found in the OPTICS ordering. }
#'
#' @author Michael Hahsler and Matthew Piekenbrock
#' @seealso Density [reachability].
#'
#' @references Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg
#' Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure.
#' _ACM SIGMOD international conference on Management of data._ ACM Press. pp.
#' \doi{10.1145/304181.304187}
#'
#' Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast Density-Based
#' Clustering with R.  _Journal of Statistical Software_, 91(1), 1-30.
#' \doi{10.18637/jss.v091.i01}
#'
#' Erich Schubert, Michael Gertz (2018). Improving the Cluster Structure
#' Extracted from OPTICS Plots. In _Lernen, Wissen, Daten, Analysen (LWDA 2018),_
#' pp. 318-329.
#' @keywords model clustering
#' @examples
#' set.seed(2)
#' n <- 400
#'
#' x <- cbind(
#'   x = runif(4, 0, 1) + rnorm(n, sd = 0.1),
#'   y = runif(4, 0, 1) + rnorm(n, sd = 0.1)
#'   )
#'
#' plot(x, col=rep(1:4, times = 100))
#'
#' ### run OPTICS (Note: we use the default eps calculation)
#' res <- optics(x, minPts = 10)
#' res
#'
#' ### get order
#' res$order
#'
#' ### plot produces a reachability plot
#' plot(res)
#'
#' ### plot the order of points in the reachability plot
#' plot(x, col = "grey")
#' polygon(x[res$order, ])
#'
#' ### extract a DBSCAN clustering by cutting the reachability plot at eps_cl
#' res <- extractDBSCAN(res, eps_cl = .065)
#' res
#'
#' plot(res)  ## black is noise
#' hullplot(x, res)
#'
#' ### re-cut at a higher eps threshold
#' res <- extractDBSCAN(res, eps_cl = .07)
#' res
#' plot(res)
#' hullplot(x, res)
#'
#' ### extract hierarchical clustering of varying density using the Xi method
#' res <- extractXi(res, xi = 0.01)
#' res
#'
#' plot(res)
#' hullplot(x, res)
#'
#' # Xi cluster structure
#' res$clusters_xi
#'
#' ### use OPTICS on a precomputed distance matrix
#' d <- dist(x)
#' res <- optics(d, minPts = 10)
#' plot(res)
#' @export
optics <- function(x, eps = Inf, minPts = 5, ...) {

  minPts <- .validate_integer_scalar(minPts, "minPts")
  eps <- .validate_nonnegative_scalar(eps, "eps", allow_infinite = TRUE)

  ### For infinity we use eps from minPts which gives the same result
  if (is.infinite(eps))
    eps <- max(kNNdist(x, k =  minPts))

  ### extra contains settings for frNN
  ### search = "kdtree", bucketSize = 10, splitRule = "suggest", approx = 0
  extra <- list(...)
  args <- c("search", "bucketSize", "splitRule", "approx")
  m <- pmatch(names(extra), args)
  if (anyNA(m))
    stop("Unknown parameter: ",
      toString(names(extra)[is.na(m)]))
  names(extra) <- args[m]

  search <- .parse_search(extra$search %||% "kdtree")
  splitRule <- .parse_splitRule(extra$splitRule %||% "suggest")
  bucketSize <- .validate_bucket_size(extra$bucketSize %||% 10L)
  approx <- .validate_nonnegative_scalar(extra$approx %||% 0, "approx")

  ### dist search
  if (search == 3L && !inherits(x, "dist")) {
    if (.matrixlike(x))
      x <- dist(x)
    else
      stop("x needs to be a matrix to calculate distances")
  }

  ## for dist we provide the R code with a frNN list and no x
  frNN <- list()
  if (inherits(x, "dist")) {
    frNN <- frNN(x, eps, ...)
    ## add self match and use C numbering
    frNN$id <- lapply(
      seq_along(frNN$id),
      FUN = function(i)
        c(i - 1L, frNN$id[[i]] - 1L)
    )
    frNN$dist <- lapply(
      seq_along(frNN$dist),
      FUN = function(i)
        c(0, frNN$dist[[i]]) ^ 2
    )

    x <- matrix()
    storage.mode(x) <- "double"

  } else{
    x <- .as_finite_numeric_matrix(x)
  }

  ret <-
    optics_int(
      as.matrix(x),
      as.double(eps),
      as.integer(minPts),
      as.integer(search),
      as.integer(bucketSize),
      as.integer(splitRule),
      as.double(approx),
      frNN
    )

  ret$minPts <- as.integer(minPts)
  ret$eps <- as.double(eps)
  ret$eps_cl <- NA_real_
  ret$xi <- NA_real_
  class(ret) <- "optics"

  ret
}

#' @rdname optics
#' @export
print.optics <- function(x, ...) {
  writeLines(c(
    paste0(
      "OPTICS ordering/clustering for ",
      length(x$order),
      " objects."
    ),
    paste0(
      "Parameters: ",
      "minPts = ",
      x$minPts,
      ", eps = ",
      x$eps,
      ", eps_cl = ",
      x$eps_cl,
      ", xi = ",
      x$xi
    )
  ))

  if (!is.null(x$cluster)) {

    if (is.na(x$xi)) {
      writeLines(paste0(
        "The clustering contains ",
        ncluster(x),
        " cluster(s) and ",
        nnoise(x),
        " noise points."
      ))

      print(table(x$cluster))
    } else {
      writeLines(
        paste0(
          "The clustering contains ",
          nrow(x$clusters_xi),
          " cluster(s) and ",
          nnoise(x),
          " noise points."
        )
      )
    }
    cat("\n")
  }
  writeLines(strwrap(paste0(
    "Available fields: ",
    toString(names(x))
  ), exdent = 18))
}

#' @rdname optics
#' @export
plot.optics <-
  function(x,
    cluster = TRUE,
    predecessor = FALSE,
    ...) {
    # OPTICS cluster extraction methods
    if (inherits(x$cluster, "xics") ||
        all(c("start", "end", "cluster_id") %in% names(x$clusters_xi))) {
      # Sort clusters by size
      hclusters <-
        x$clusters_xi[order(x$clusters_xi$end - x$clusters_xi$start), ]

      # .1 means to leave 15% for the cluster lines
      def.par <- par(no.readonly = TRUE)
      par(mar = c(2, 4, 4, 2) + 0.1, omd = c(0, 1, .15, 1))

      # Need to know how to spread out lines
      y_max <- max(x$reachdist[!is.infinite(x$reachdist)])
      y_increments <- (y_max / 0.85 * .15) / (nrow(hclusters) + 1L)

      # Get top level cluster labels
      # top_level <- extractClusterLabels(x$clusters_xi, x$order)
      plot(
        as.reachability(x),
        col = x$cluster[x$order] + 1L,
        xlab = NA,
        xaxt = 'n',
        yaxs = "i",
        ylim = c(0, y_max),
        ...
      )

      # Lines beneath plotting region indicating Xi clusters
      i <- seq_len(nrow(hclusters))
      segments(
        x0 = hclusters$start[i],
        y0 = -(y_increments * i),
        x1 = hclusters$end[i],
        col = hclusters$cluster_id[i] + 1L,
        lwd = 2,
        xpd = NA
      )
      ## Restore previous settings
      par(def.par)
    } else if (is.numeric(x$cluster) &&
        !is.null(x$eps_cl)) {
      # Works for integers too
      ## extractDBSCAN clustering
      plot(as.reachability(x), col = x$cluster[x$order] + 1L, ...)
      lines(
        x = c(0, length(x$cluster)),
        y = c(x$eps_cl, x$eps_cl),
        col = "black",
        lty = 2
      )
    } else {
      # Regular reachability plot
      plot(as.reachability(x), ...)
    }
  }

# Simple conversion between OPTICS objects and reachability objects
#' @rdname optics
#' @export
as.reachability.optics <- function(object, ...) {
  structure(list(reachdist = object$reachdist, order = object$order),
    class = "reachability")
}

# Conversion between OPTICS objects and dendrograms
#' @rdname optics
#' @export
as.dendrogram.optics <- function(object, ...) {
  if (object$minPts > length(object$order)) {
    stop("'minPts' should be less or equal to the points in the dataset.")
  }
  if (sum(is.infinite(object$reachdist)) > 1)
    stop(
      "Eps value is not large enough to capture the complete hierarchical structure of the dataset. Please use a large eps value (such as Inf)."
    )
  as.dendrogram(as.reachability(object))
}

#' @rdname optics
#' @export
extractDBSCAN <- function(object, eps_cl) {
  if (!inherits(object, "optics"))
    stop("extractDBSCAN only accepts objects resulting from dbscan::optics!")

  reachdist <- object$reachdist[object$order]
  coredist <- object$coredist[object$order]
  n <- length(object$order)
  cluster <- integer(n)

  clusterid <- 0L         ### 0 is noise
  for (i in seq_len(n)) {
    if (reachdist[i] > eps_cl) {
      if (coredist[i] <= eps_cl) {
        clusterid <- clusterid + 1L
        cluster[i] <- clusterid
      } else{
        cluster[i] <- 0L  ### noise
      }
    } else{
      cluster[i] <- clusterid
    }
  }

  object$eps_cl <- eps_cl
  object$xi <- NA_real_
  ### fix the order so cluster is in the same order as the rows in x
  cluster[object$order] <- cluster
  object$cluster <- cluster

  object
}

### Extract clusters (minimum == T extracts clusters that do not contain other clusters) from a given ordering of points
extractClusterLabels <- function(cl, order, minimum = FALSE) {
  ## Add cluster_id to clusters
  if (!all(c("start", "end") %in% names(cl)))
    stop("extractClusterLabels expects start and end references")
  if (!"cluster_id" %in% names(cl))
    cl <- cbind(cl, cluster_id = seq_len(nrow(cl)))

  ## Sort cl based on minimum parameter / cluster size
  if (!"cluster_size" %in% names(cl))
    cl <- cbind(cl, list(cluster_size = (cl$end - cl$start)))
  cl <-
    if (minimum) {
      cl[order(cl$cluster_size), ]
    } else {
      cl[order(-cl$cluster_size), ]
    }

  ## Fill in the [cluster] vector with cluster IDs
  clusters <- integer(length(order))
  for (cid in cl$cluster_id) {
    cluster <- cl[cl$cluster_id == cid, ]
    if (minimum) {
      if (all(clusters[cluster$start:cluster$end] == 0)) {
        clusters[cluster$start:cluster$end] <- cid
      }
    } else
      clusters[cluster$start:cluster$end] <- cid
  }

  # Fix the ordering
  clusters[order] <- clusters
  return(clusters)
}
