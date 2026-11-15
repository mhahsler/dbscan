#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2017 Michael Hahsler

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

# number of shared nearest neighbors including the point itself.


#' Find Shared Nearest Neighbors
#'
#' Calculates the number of shared nearest neighbors
#' and creates a shared nearest neighbors graph.
#'
#' The number of shared nearest neighbors of two points p and q is the size of
#' the intersection of the kNN neighborhoods of the two points.
#' Note that each point is considered to be part
#' of its own kNN neighborhood.
#' The range for the shared nearest neighbors is
#' \eqn{[0, k]}. The result is a n-by-k matrix called `shared`.
#' Each row is a point and the columns are the point's k nearest neighbors.
#' The value `shared[i, j]` is the count of the shared neighbors of point `i` and
#' its j-th neighbor.
#'
#' The shared nearest neighbor graph connects a point with all its nearest neighbors
#' if they have at least one shared neighbor. The number of shared neighbors can be used
#' as an edge weight.
#' Jarvis and Patrick (1973) use a slightly
#' modified (see parameter `jp`) shared nearest neighbor graph for
#' clustering.
#'
#' The `k` nearest neighbors are found with [kNN()]. The search parameters are
#' passed on to it (see the section on the search parameters below).
#'
#' ## Using a kNN object as `x`
#'
#' If `x` is a `kNN` object, then `k` can be omitted and `x$k` is used. If `k`
#' is smaller than `x$k`, then the first `k` neighbors (by distance) are used.
#' An error is raised if `x` contains fewer than `k` neighbors. The search
#' parameters are ignored. A `query` is not supported.
#'
#' ## Thresholding and sorting
#'
#' * **Threshold `kt`:** If `kt` is specified, then only edges with at least
#'   `kt` shared neighbors are kept (`kt` needs to be at most `k`). The entries
#'   for removed edges are set to `NA` in `id`, `dist` and `shared`.
#'   [adjacencylist()] skips them.
#' * **Sorting:** The neighbors are found without sorting. With `sort = TRUE`
#'   (default), the neighbors of each point (the columns of `id`, `dist` and
#'   `shared`) are then reordered by decreasing number of shared neighbors
#'   (ties are broken by id). The columns then no longer follow the distance
#'   rank, and `$sort` is `FALSE`. A `sNN` object can be sorted later using
#'   `sort()`.
#'
#' @aliases sNN snn
#' @family NN functions
#'
#' @templateVar type sNN
#' @template NN-search-section
#' @template NN-search-params
#'
#' @param x a data matrix, a [dist] object or a [kNN] object.
#' @param k number of neighbors to consider to calculate the shared nearest
#' neighbors.
#' @param kt minimum threshold on the number of shared nearest neighbors to
#' build the shared nearest neighbor graph. Edges are only preserved if
#' `kt` or more neighbors are shared.
#' @param jp In regular sNN graphs, a point and one of its neighbors can have
#' shared neighbors even if the point is not one of the neighbor's `k` nearest
#' neighbors.
#' Jarvis and Patrick (1973) require the neighbor relation to be mutual: the
#' point also needs to be among the neighbor's `k` nearest neighbors, otherwise
#' the count is zeroed out. `TRUE` uses this behavior.
#' @param sort sort the neighbors of each point by the number of shared nearest
#' neighbors (decreasing, ties are broken by id)? Note that this
#' is expensive and `sort = FALSE` is much faster. sNN objects can be
#' sorted later using `sort()`.
#' @param decreasing for `sort()`, sort by decreasing number of shared
#' neighbors (the default)? Use `FALSE` for increasing order.
#' @param ... further arguments (ignored).
#' @return An object of class `sNN` (subclass of [kNN] and [NN]) containing a list
#' with the following components (the matrices have one row per point and `k`
#' columns; edges removed by `kt` are `NA`):
#' \item{id }{a matrix with the ids (row numbers in `x`) of the neighbors. }
#' \item{dist}{a matrix with the distances to the neighbors. }
#' \item{shared }{a matrix with the number of shared nearest neighbors. }
#' \item{k }{the number `k` used. }
#' \item{sort }{logical; are the neighbors sorted by increasing distance? This
#'   is `FALSE` if `sort = TRUE` was used, since the neighbors are then sorted by
#'   the number of shared neighbors. }
#' \item{sort_shared }{logical; are the neighbors sorted by decreasing number of
#'   shared neighbors? }
#' \item{kt }{the threshold `kt` used (`NULL` if none). }
#' \item{metric }{the used distance metric. }
#'
#' @author Michael Hahsler
#' @references R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a
#' Similarity Measure Based on Shared Near Neighbors. _IEEE Trans. Comput._
#' 22, 11 (November 1973), 1025-1034.
#' \doi{10.1109/T-C.1973.223640}
#'
#' David M. Mount and Sunil Arya (2010). ANN: A Library for
#' Approximate Nearest Neighbor Searching,
#' \url{http://www.cs.umd.edu/~mount/ANN/}.
#' @seealso [sNNclust()] and [jpclust()], which cluster using the shared nearest
#' neighbor graph.
#' @keywords model
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # finding kNN and add the number of shared nearest neighbors.
#' k <- 5
#' nn <- sNN(x, k = k)
#' nn
#'
#' # shared nearest neighbor distribution
#' table(as.vector(nn$shared))
#'
#' # explore number of shared points for the k-neighborhood of point 10
#' i <- 10
#' nn$shared[i,]
#'
#' plot(nn, x)
#'
#' # apply a threshold to create a sNN graph with edges
#' # if more than 3 neighbors are shared.
#' nn_3 <- sNN(nn, kt = 3)
#' plot(nn_3, x)
#'
#' # get an adjacency list for the shared nearest neighbor graph
#' adjacencylist(nn_3)
#'
#' # Jarvis-Patrick graph (requires mutual nearest neighbors)
#' nn_jp <- sNN(x, k = k, jp = TRUE)
#' table(as.vector(nn_jp$shared))
#' @export
sNN <- function(x,
  k,
  kt = NULL,
  jp = FALSE,
  sort = TRUE,
  search = "kdtree",
  bucketSize = 10,
  splitRule = "suggest",
  shrinkRule = "suggest",
  approx = 0) {
  if (missing(k)) {
    if (!inherits(x, "kNN"))
      stop("k is required unless x is a kNN object.")
    k <- x$k
  }
  k <- .validate_integer_scalar(k, "k")

  if (inherits(x, "kNN")) {
    if (k != x$k) {
      if (ncol(x$id) < k)
        stop("kNN object does not contain enough neighbors!")
      if (!x$sort)
        x <- sort.kNN(x)
      x$id <- x$id[, 1:k, drop = FALSE]
      x$dist <- x$dist[, 1:k, drop = FALSE]
      x$k <- k
    }

  } else
    x <-
      kNN(
        x,
        k,
        sort = FALSE,
        search = search,
        bucketSize = bucketSize,
        splitRule = splitRule,
        shrinkRule = shrinkRule,
        approx = approx
      )

  x$shared <- SNN_sim_int(x$id, as.logical(jp[1]))
  x$sort_shared <- FALSE

  class(x) <- c("sNN", "kNN", "NN")

  if (sort)
    x <- sort.sNN(x)

  x$kt <- kt

  if (!is.null(kt)) {
    if (kt > k)
      stop("kt needs to be at most k.")
    rem <- x$shared < kt
    x$id[rem] <- NA
    x$dist[rem] <- NA
    x$shared[rem] <- NA
  }

  x
}

#' @rdname sNN
#' @export
sort.sNN <- function(x, decreasing = TRUE, ...) {
  ## x$sort_shared = TRUE means: sorted by decreasing number of shared neighbors
  if (isTRUE(x$sort_shared) && decreasing)
    return(x)
  if (is.null(x$shared))
    stop("Unable to sort. Number of shared neighbors is missing.")
  if (ncol(x$id) < 2) {
    x$sort <- TRUE
    x$sort_shared <- TRUE
    return(x)
  }

  ## sort by number of shared points and break ties by id (increasing)
  k <- ncol(x$shared)
  sgn <- if (decreasing) -1L else 1L
  o <- vapply(
    seq_len(nrow(x$shared)),
    function(i) order(sgn * x$shared[i, ], x$id[i, ]),
    integer(k)
  )
  for (i in seq_len(ncol(o))) {
    x$shared[i, ] <- x$shared[i, ][o[, i]]
    x$dist[i, ] <- x$dist[i, ][o[, i]]
    x$id[i, ] <- x$id[i, ][o[, i]]
  }

  x$sort <- FALSE
  x$sort_shared <- decreasing

  x
}

#' @rdname sNN
#' @export
print.sNN <- function(x, ...) {
  cat(
    "shared-nearest neighbors for ",
    nrow(x$id),
    " objects (k=",
    x$k,
    ", kt=",
    x$kt %||% "NULL",
    ").",
    "\n",
    sep = ""
  )
  cat("Available fields: ", toString(names(x)), "\n", sep = "")
}
