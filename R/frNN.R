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


#' Find the Fixed Radius Nearest Neighbors
#'
#' This function finds the fixed radius nearest neighbors
#' (including distances) fast.
#'
#' If `x` is specified as a data matrix, then Euclidean distances and fast
#' nearest neighbor lookup using a k-d tree (default), a box-decomposition tree
#' or linear search are used (see the section on the search parameters).
#' All points within a distance of `eps` of a point (including points at
#' exactly distance `eps`) are its neighbors.
#'
#' ## Sorting and self-matches
#'
#' * **Sorting:** With `sort = TRUE` (default), the neighbors of each point are
#'   ordered by increasing distance and ties are broken by id. With
#'   `sort = FALSE`, the neighbors are returned in no particular order, which
#'   avoids the cost of sorting. A `frNN` object can be sorted later using `sort()`.
#' * **Self-matches:** If no `query` is specified, then self-matches are
#'   removed. If `query` is specified, then a query point that is also in `x`
#'   is returned as its own neighbor (at distance 0).
#'
#' ## Using a frNN object as `x`
#'
#' If `x` is a `frNN` object, then the neighborhoods are reduced to the
#' smaller radius `eps`. An error is raised if `eps` is larger than the radius
#' stored in `x`. All other arguments are ignored.
#'
#' ## Input requirements
#'
#' * `x` can be a numeric matrix, a data.frame with only numeric columns, a
#'   [dist] object or a `frNN` object.
#' * `query` needs to be a numeric matrix or data.frame with the same number of
#'   columns as `x`, and can only be used if `x` is a data matrix.
#'
#' ## Memory
#'
#' The size of the result grows with the number of neighbors of all points, which
#' can approach \eqn{O(n^2)} for a large `eps`. Choose `eps` carefully for large
#' data sets.
#'
#' ## Creating a frNN object from scratch
#'
#' A `frNN` object is a list with class `c("frNN", "NN")` and the components
#' described in the value section. You need to supply at least `id` and `eps`.
#' `dist` is needed for sorting and plotting, and `sort` should be set to
#' `FALSE` unless the neighbors are already sorted. `metric` is informational.
#'
#' @aliases frNN frnn print.frnn
#' @family NN functions
#'
#' @templateVar type frNN
#' @template NN-search-section
#' @template NN-search-params
#'
#' @param x a data matrix, a dist object or a frNN object.
#' @param eps radius of the neighborhood (a finite, non-negative number).
#' Points within this distance are neighbors.
#' @param query a data matrix with the points to query. If query is not
#' specified, the NN for all the points in `x` is returned. If query is
#' specified then `x` needs to be a data matrix.
#' @param sort sort the neighbors of each point by distance (ties are broken by
#' id)? This is expensive and can be done later using `sort()`.
#' @param decreasing for `sort()`, sort the neighbors by decreasing distance
#' instead of increasing distance?
#' @param ... further arguments (ignored).
#'
#' @returns
#'
#' `frNN()` returns an object of class [frNN] (subclass of
#' [NN]) containing a list with the following components:
#' \item{id }{a list with one integer vector per point in `x` (or `query`).
#' Each vector contains the ids (row numbers in `x`) of the fixed radius nearest
#' neighbors. A point without neighbors has an empty vector (`integer(0)`).
#' The list is named with the row names if available. }
#' \item{dist }{a list with the distances (same structure as
#' `id`). }
#' \item{eps }{ neighborhood radius `eps` that was used. }
#' \item{sort }{logical; are the neighbors of each point sorted by increasing
#' distance (with ties broken by id)? }
#' \item{metric }{ used distance metric. This is `"euclidean"` for data
#' matrices. For a [dist] object, it is the method stored in the object (which
#' can be `NULL`). }
#'
#' `adjacencylist()` returns a list with one entry per data point in `x`. Each entry
#' contains the id of the nearest neighbors.
#'
#' @author Michael Hahsler
#'
#' @references David M. Mount and Sunil Arya (2010). ANN: A Library for
#' Approximate Nearest Neighbor Searching,
#' \url{http://www.cs.umd.edu/~mount/ANN/}.
#' @seealso [dbscan()], which accepts a precomputed `frNN` object.
#' @keywords model
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # Example 1: Find fixed radius nearest neighbors for each point
#' nn <- frNN(x, eps = .5)
#' nn
#'
#' # Number of neighbors
#' hist(lengths(adjacencylist(nn)),
#'   xlab = "Number of neighbors", main = "Number of Neighbors",
#'   sub = paste("Neighborhood size eps =", nn$eps))
#'
#' # Explore neighbors of point i = 10
#' i <- 10
#' nn$id[[i]]
#' nn$dist[[i]]
#' plot(x, col = ifelse(seq_len(nrow(iris)) %in% nn$id[[i]], "red", "black"))
#'
#' # get an adjacency list
#' head(adjacencylist(nn))
#'
#' # plot the fixed radius neighbors (and then reduced to a radius of .3)
#' plot(nn, x)
#' plot(frNN(nn, eps = .3), x)
#'
#' ## Example 2: find fixed-radius NN for query points
#' q <- x[c(1,100),]
#' nn <- frNN(x, eps = .5, query = q)
#'
#' plot(nn, x[, 1:2], col = "grey")
#' points(q, pch = 3, lwd = 2)
#'
#' ## Example 3: a precomputed frNN object can be used for DBSCAN
#' nn <- frNN(x, eps = .5)
#' db <- dbscan(nn, minPts = 5)
#' db
#'
#' clplot(x, cl = db)
#' @export frNN
frNN <-
  function(x,
    eps,
    query = NULL,
    sort = TRUE,
    search = "kdtree",
    bucketSize = 10,
    splitRule = "suggest",
    shrinkRule = "suggest",
    approx = 0) {
    eps <- .validate_nonnegative_scalar(eps, "eps")

    if (inherits(x, "frNN")) {
      if (x$eps < eps)
        stop("frNN in x has not a sufficient eps radius.")

      for (i in seq_along(x$dist)) {
        take <- x$dist[[i]] <= eps
        x$dist[[i]] <- x$dist[[i]][take]
        x$id[[i]] <- x$id[[i]][take]
      }
      x$eps <- eps

      return(x)
    }

    search <- .parse_search(search)
    splitRule <- .parse_splitRule(splitRule)
    shrinkRule <- .parse_shrinkRule(shrinkRule)
    bucketSize <- .validate_bucket_size(bucketSize)
    approx <- .validate_nonnegative_scalar(approx, "approx")

    ### dist search
    if (search == 3 && !inherits(x, "dist")) {
      if (.matrixlike(x))
        x <- dist(x)
      else
        stop("x needs to be a matrix to calculate distances")
    }

    ### get kNN from a dist object in R
    if (inherits(x, "dist")) {
      if (!is.null(query))
        stop("query can only be used if x contains the data.")

      if (anyNA(x))
        stop("data/distances cannot contain NAs for frNN (with k-d tree)!")

      return(dist_to_frNN(x, eps = eps, sort = sort))
    }

    ## make sure x is numeric
    x <- .as_finite_numeric_matrix(x)

    if (!is.null(query)) {
      query <- .as_finite_numeric_matrix(query, "query")
      if (ncol(x) != ncol(query))
        stop("x and query need to have the same number of columns!")
    }

    ## returns NO self matches
    if (!is.null(query)) {
      ret <-
        frNN_query_int(
          as.matrix(x),
          as.matrix(query),
          as.double(eps),
          as.integer(search),
          as.integer(bucketSize),
          as.integer(splitRule),
          as.integer(shrinkRule),
          as.double(approx)
        )
      names(ret$dist) <- rownames(query)
      names(ret$id) <- rownames(query)
      ret$metric <- "euclidean"
    } else {
      ret <- frNN_int(
        as.matrix(x),
        as.double(eps),
        as.integer(search),
        as.integer(bucketSize),
        as.integer(splitRule),
        as.integer(shrinkRule),
        as.double(approx)
      )
      names(ret$dist) <- rownames(x)
      names(ret$id) <- rownames(x)
      ret$metric <- "euclidean"
    }

    ret$eps <- eps
    ret$sort <- FALSE
    class(ret) <- c("frNN", "NN")

    if (sort)
      ret <- sort.frNN(ret)

    ret
  }

# extract a row from a distance matrix without doubling space requirements
dist_row <- function(x, i, self_val = 0) {
  n <- attr(x, "Size")

  i <- rep(i, times = n)
  j <- seq_len(n)
  swap_idx <- i > j
  tmp <- i[swap_idx]
  i[swap_idx] <- j[swap_idx]
  j[swap_idx] <- tmp

  diag_idx <- i == j
  idx <- n * (i - 1) - i * (i - 1) / 2 + j - i
  idx[diag_idx] <- NA

  val <- x[idx]
  val[diag_idx] <- self_val
  val
}

dist_to_frNN <- function(x, eps, sort = FALSE) {
  .check_dist(x)

  n <- attr(x, "Size")

  id <- list()
  d <- list()

  for (i in seq_len(n)) {
    ### Inf -> no self-matches
    y <- dist_row(x, i, self_val = Inf)
    o <- which(y <= eps)
    id[[i]] <- o
    d[[i]] <- y[o]
  }
  names(id) <- labels(x)
  names(d) <- labels(x)

  ret <-
    structure(list(
      dist = d,
      id = id,
      eps = eps,
      metric = attr(x, "method"),
      sort = FALSE
    ),
      class = c("frNN", "NN"))

  if (sort)
    ret <- sort.frNN(ret)

  return(ret)
}

#' @rdname frNN
#' @export
sort.frNN <- function(x, decreasing = FALSE, ...) {
  ## x$sort = TRUE means: sorted by increasing distance (ties by id)
  if (isTRUE(x$sort) && !decreasing)
    return(x)
  if (is.null(x$dist))
    stop("Unable to sort. Distances are missing.")

  ## FIXME: This is slow do this in C++
  n <- names(x$id)

  o <- lapply(
    seq_along(x$dist),
    FUN =
      function(i)
        order(x$dist[[i]], x$id[[i]], decreasing = decreasing)
  )
  x$dist <-
    lapply(
      seq_along(o),
      FUN = function(p)
        x$dist[[p]][o[[p]]]
    )
  x$id <- lapply(
    seq_along(o),
    FUN = function(p)
      x$id[[p]][o[[p]]]
  )

  names(x$dist) <- n
  names(x$id) <- n

  x$sort <- !decreasing

  x
}

#' @rdname frNN
#' @export
adjacencylist.frNN <- function(x, ...)
  x$id

#' @rdname frNN
#' @export
print.frNN <- function(x, ...) {
  cat(
    "fixed radius nearest neighbors for ",
    length(x$id),
    " objects (eps=",
    x$eps,
    ").",
    "\n",
    sep = ""
  )

  cat("Distance metric:", x$metric, "\n")
  cat("\nAvailable fields: ", toString(names(x)), "\n", sep = "")
}
