#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2017 Michael Hahsler

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

# number of shared nearest neighbors including the point itself.


#' Find Shared Nearest Neighbors
#'
#' Calculates the number of shared nearest neighbors
#' and creates a shared nearest neighbors graph.
#'
#' The number of shared nearest neighbors of two points p and q is the
#' intersection of the kNN neighborhood of two points.
#' Note: that each point is considered to be part
#' of its own kNN neighborhood.
#' The range for the shared nearest neighbors is
#' \eqn{[0, k]}. The result is a n-by-k matrix called `shared`.
#' Each row is a point and the columns are the point's k nearest neighbors.
#' The value is the count of the shared neighbors.
#'
#' The shared nearest neighbor graph connects a point with all its nearest neighbors
#' if they have at least one shared neighbor. The number of shared neighbors can be used
#' as an edge weight.
#' Javis and Patrick (1973) use a slightly
#' modified (see parameter `jp`) shared nearest neighbor graph for
#' clustering.
#'
#' @aliases sNN snn
#' @family NN functions
#'
#' @param x a data matrix, a [dist] object or a [kNN] object.
#' @param k number of neighbors to consider to calculate the shared nearest
#' neighbors.
#' @param kt minimum threshold on the number of shared nearest neighbors to
#' build the shared nearest neighbor graph. Edges are only preserved if
#' `kt` or more neighbors are shared.
#' @param jp In regular sNN graphs, two points that are not neighbors
#' can have shared neighbors.
#' Javis and Patrick (1973) requires the two points to be neighbors, otherwise
#' the count is zeroed out. `TRUE` uses this behavior.
#' @param search nearest neighbor search strategy (one of `"kdtree"`, `"linear"` or
#' `"dist"`).
#' @param sort sort by the number of shared nearest neighbors? Note that this
#' is expensive and `sort = FALSE` is much faster. sNN objects can be
#' sorted using `sort()`.
#' @param bucketSize max size of the kd-tree leafs.
#' @param splitRule rule to split the kd-tree. One of `"STD"`, `"MIDPT"`, `"FAIR"`,
#' `"SL_MIDPT"`, `"SL_FAIR"` or `"SUGGEST"` (SL stands for sliding). `"SUGGEST"` uses
#' ANNs best guess.
#' @param approx use approximate nearest neighbors. All NN up to a distance of
#' a factor of `(1 + approx) eps` may be used. Some actual NN may be omitted
#' leading to spurious clusters and noise points.  However, the algorithm will
#' enjoy a significant speedup.
#' @param decreasing logical; sort in decreasing order?
#' @param ... additional parameters are passed on.
#' @return An object of class `sNN` (subclass of [kNN] and [NN]) containing a list
#' with the following components:
#' \item{id }{a matrix with ids. }
#' \item{dist}{a matrix with the distances. }
#' \item{shared }{a matrix with the number of shared nearest neighbors. }
#' \item{k }{number of `k` used. }
#' \item{metric }{the used distance metric. }
#'
#' @author Michael Hahsler
#' @references R. A. Jarvis and E. A. Patrick. 1973. Clustering Using a
#' Similarity Measure Based on Shared Near Neighbors. _IEEE Trans. Comput._
#' 22, 11 (November 1973), 1025-1034.
#' \doi{10.1109/T-C.1973.223640}
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
#' @export
sNN <- function(x,
  k,
  kt = NULL,
  jp = FALSE,
  sort = TRUE,
  search = "kdtree",
  bucketSize = 10,
  splitRule = "suggest",
  approx = 0) {
  if (missing(k))
    k <- x$k

  if (inherits(x, "kNN")) {
    if (k != x$k) {
      if (ncol(x$id) < k)
        stop("kNN object does not contain enough neighbors!")
      if (!x$sort)
        x <- sort.kNN(x)
      x$id <- x$id[, 1:k]
      x$dist <- x$dist[, 1:k]
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
      stop("kt needs to be less than k.")
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
  if (isTRUE(x$sort_shared))
    return(x)
  if (is.null(x$shared))
    stop("Unable to sort. Number of shared neighbors is missing.")
  if (ncol(x$id) < 2) {
    x$sort <- TRUE
    x$sort_shared <- TRUE
    return(x)
  }

  ## sort first by number of shared points (decreasing) and break ties by id (increasing)
  k <- ncol(x$shared)
  o <- vapply(
    seq_len(nrow(x$shared)),
    function(i) order(k - x$shared[i, ], x$id[i, ], decreasing = !decreasing),
    integer(k)
  )
  for (i in seq_len(ncol(o))) {
    x$shared[i, ] <- x$shared[i, ][o[, i]]
    x$dist[i, ] <- x$dist[i, ][o[, i]]
    x$id[i, ] <- x$id[i, ][o[, i]]
  }

  x$sort <- FALSE
  x$sort_shared <- TRUE

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
    ifelse(is.null(x$kt), "NULL", x$kt),
    ").",
    "\n",
    sep = ""
  )
  cat("Available fields: ", toString(names(x)), "\n", sep = "")
}
