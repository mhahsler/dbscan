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


#' Shared Nearest Neighbors
#'
#' Calculates the number of shared nearest neighbors, the shared nearest
#' neighbor similarity and creates a shared nearest neighbors graph.
#'
#' The number of shared nearest neighbors is the intersection of the kNN
#' neighborhood of two points. Note: that each point is considered to be part
#' of its own kNN neighborhood. The range for the shared nearest neighbors is
#' \eqn{[0, k]}.
#'
#' Javis and Patrick (1973) use the shared nearest neighbor graph for
#' clustering. They only count shared neighbors between points that are in each
#' other's kNN neighborhood.
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
#' @param jp use the definition by Javis and Patrick (1973), where shared
#' neighbors are only counted between points that are in each other's
#' neighborhood, otherwise 0 is returned. If `FALSE`, then the number of shared
#' neighbors is returned, even if the points are not neighbors.
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
#' @return An object of class `sNN` (subclass of [kNN] and [NN]) containing a list
#' with the following components:
#' \item{id }{a matrix with ids. }
#' \item{dist}{a matrix with the distances. }
#' \item{shared }{a matrix with the number of shared nearest neighbors. }
#' \item{k }{number of `k` used. }
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
#' # explore neighborhood of point 10
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
#' @export sNN
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


sort.sNN <- function(x, decreasing = TRUE, ...) {
  if (!is.null(x$sort_shared) && x$sort_shared)
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
  o <- sapply(
    1:nrow(x$shared),
    FUN =
      function(i)
        order(k - x$shared[i, ], x$id[i, ], decreasing = !decreasing)
  )
  for (i in 1:ncol(o)) {
    x$shared[i, ] <- x$shared[i, ][o[, i]]
    x$dist[i, ] <- x$dist[i, ][o[, i]]
    x$id[i, ] <- x$id[i, ][o[, i]]
  }

  x$sort <- FALSE
  x$sort_shared <- TRUE

  x
}

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
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}
