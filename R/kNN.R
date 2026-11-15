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


#' Find the k Nearest Neighbors
#'
#' This function finds all k nearest neighbors in a data
#' matrix (including distances).
#'
#' To speed up k nearest neighbors, a k-d tree is used
#' for data matrices by default. Box-decomposition trees, linear search,
#' and precomputing distance matrices are also available.
#'
#' ## Using a kNN object as `x`
#'
#' If `x` is a `kNN` object, then the first `k` neighbors are returned
#' (sorted by distance). All other arguments are ignored. An error is
#' raised if `x` contains fewer than `k` neighbors.
#'
#' ## Sorting, ties and self-matches
#'
#' * **Sorting:** With `sort = TRUE` (default), the neighbors of each point are
#'   ordered by distance and ties are broken by id. With `sort = FALSE`, the
#'   neighbors are returned in the order found by the search (the tree-based
#'   searches return them by increasing distance, but tied neighbors are
#'   not ordered by id), which avoids the cost of sorting.
#'   A `kNN` object can be sorted later using `sort()`.
#' * **Ties:** If the kth and the (k+1)th nearest neighbor are tied, then the
#'   neighbor found first is returned and the other one is ignored.
#' * **Self-matches:** If no `query` is specified, then self-matches are
#'   removed. If `query` is specified, then a query point that is also in `x`
#'   is returned as its own nearest neighbor (at distance 0).
#'
#' ## Input requirements
#'
#' * `x` can be a numeric matrix, a data.frame with only numeric columns, a
#'   [dist] object or a `kNN` object.
#' * `query` needs to be a numeric matrix or
#'   data.frame with the same number of columns as `x`, and can only be used
#'   if `x` is a data matrix.
#' * `k` needs to be a whole number with `1 <= k < nrow(x)`.
#'
#' ## Non-Euclidean distances
#'
#' To use non-Euclidean distances, a distance matrix needs to be precomputed
#' by the user and supplied as `x`. Note that this leads to restrictions
#' of the number of observations that can be used since the distance
#' matrix grows with \eqn{O(n^2)}.
#'
#' @aliases kNN knn
#' @family NN functions
#'
#' @templateVar type kNN
#' @template NN-search-section
#' @template NN-search-params
#'
#' @param x a numeric data matrix (or data.frame), a [dist] object or a [kNN]
#' object.
#' @param k number of neighbors to find.
#' @param query a data matrix with the points to query. If query is not
#' specified, the NN for all the points in `x` is returned. If query is
#' specified then `x` needs to be a data matrix.
#' @param sort sort the neighbors of each point by distance (ties are broken by
#' id)? `sort = FALSE` skips this step and returns the neighbors in the order
#' found by the search, which is faster. kNN objects can be sorted later using
#' `sort()`.
#' @param decreasing for `sort()`, sort the neighbors by decreasing distance
#' instead of increasing distance?
#' @param ... further arguments (ignored).
#'
#' @return An object of class `kNN` (subclass of [NN]) containing a
#' list with the following components:
#' \item{dist }{a matrix with distances (one row per point in `x` or `query`,
#'   one column per neighbor). }
#' \item{id }{a matrix with the row numbers in `x` of the neighbors. }
#' \item{k }{number `k` used. }
#' \item{sort }{logical; are the neighbors in each row sorted by increasing
#'   distance (with ties broken by id)? }
#' \item{metric }{ used distance metric. This is `"euclidean"` for data
#'   matrices. For a [dist] object, it is the method stored in the object (which
#'   can be `NULL`). }
#'
#' @author Michael Hahsler
#' @references David M. Mount and Sunil Arya (2010). ANN: A Library for
#' Approximate Nearest Neighbor Searching,
#' \url{https://www.cs.umd.edu/~mount/ANN/}.
#' @keywords model
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # Example 1: finding kNN for all points in a data matrix (using a k-d tree)
#' nn <- kNN(x, k = 5)
#' nn
#'
#' # explore neighborhood of point 10
#' i <- 10
#' nn$id[i,]
#' plot(x, col = ifelse(seq_len(nrow(iris)) %in% nn$id[i,], "red", "black"))
#'
#' # visualize the 5 nearest neighbors
#' plot(nn, x)
#'
#' # visualize a reduced 2-NN graph
#' plot(kNN(nn, k = 2), x)
#'
#' # Example 2: find kNN for query points
#' q <- x[c(1,100),]
#' nn <- kNN(x, k = 10, query = q)
#'
#' # we plot the first two features for iris
#' plot(nn, x[, 1:2], col = "grey")
#' points(q, pch = 3, lwd = 2)
#'
#' # Example 3: find kNN using a precomputed (non-Euclidean) distance matrix
#' d <- dist(x, method = "manhattan")
#' nn <- kNN(d, k = 1)
#' plot(nn, x)
#' @export
kNN <-
  function(x,
    k,
    query = NULL,
    sort = TRUE,
    search = "kdtree",
    bucketSize = 10,
    splitRule = "suggest",
    shrinkRule = "suggest",
    approx = 0) {
    k <- .validate_integer_scalar(k, "k")

    if (inherits(x, "kNN")) {
      if (x$k < k)
        stop("kNN in x has not enough nearest neighbors.")
      if (!x$sort)
        x <- sort(x)
      x$id <- x$id[, 1:k, drop = FALSE]
      if (!is.null(x$dist))
        x$dist <- x$dist[, 1:k, drop = FALSE]
      if (!is.null(x$shared))
        x$shared <- x$shared[, 1:k, drop = FALSE]
      x$k <- k
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

    ### get kNN from a dist object
    if (inherits(x, "dist")) {
      if (!is.null(query))
        stop("query can only be used if x contains a data matrix.")

      if (anyNA(x))
        stop("distances cannot be NAs for kNN!")

      return(dist_to_kNN(x, k = k))
    }

    ## make sure x is numeric
    x <- .as_finite_numeric_matrix(x)

    if (!is.null(query)) {
      query <- .as_finite_numeric_matrix(query, "query")
      if (ncol(x) != ncol(query))
        stop("x and query need to have the same number of columns!")
    }

    if (k >= nrow(x))
      stop("Not enough neighbors in data set!")

    ## returns NO self matches
    if (!is.null(query)) {
      ret <- kNN_query_int(
        as.matrix(x),
        as.matrix(query),
        as.integer(k),
        as.integer(search),
        as.integer(bucketSize),
        as.integer(splitRule),
        as.integer(shrinkRule),
        as.double(approx)
      )
      dimnames(ret$dist) <- list(rownames(query), 1:k)
      dimnames(ret$id) <- list(rownames(query), 1:k)
    } else {
      ret <- kNN_int(
        as.matrix(x),
        as.integer(k),
        as.integer(search),
        as.integer(bucketSize),
        as.integer(splitRule),
        as.integer(shrinkRule),
        as.double(approx)
      )
      dimnames(ret$dist) <- list(rownames(x), 1:k)
      dimnames(ret$id) <- list(rownames(x), 1:k)
    }

    class(ret) <- c("kNN", "NN")

    ### ANN already returns them sorted (by dist but not by ID)
    if (sort)
      ret <- sort(ret)

    ret$metric <- "euclidean"

    ret
  }

# make sure we have a lower-triangle representation w/o diagonal
.check_dist <- function(x) {
  if (!inherits(x, "dist"))
    stop("x needs to be a dist object")

  # cluster::dissimilarity does not have Diag or Upper attributes, but is a lower triangle
  # representation
  if (inherits(x, "dissimilarity"))
    return(TRUE)

  # check that dist objects have diag = FALSE, upper = FALSE
  if (attr(x, "Diag") || attr(x, "Upper"))
    stop("x needs to be a dist object with attributes Diag and Upper set to FALSE. Use as.dist(x, diag = FALSE, upper = FALSE) first.")
  }

dist_to_kNN <- function(x, k) {
  .check_dist(x)

  n <- attr(x, "Size")

  id <- structure(integer(n * k), dim = c(n, k))
  d <- matrix(NA_real_, nrow = n, ncol = k)

  for (i in seq_len(n)) {
    ### Inf -> no self-matches
    y <- dist_row(x, i, self_val = Inf)
    o <- order(y, decreasing = FALSE)
    o <- o[seq_len(k)]
    id[i, ] <- o
    d[i, ] <- y[o]
  }
  dimnames(id) <- list(labels(x), seq_len(k))
  dimnames(d) <- list(labels(x), seq_len(k))

  ret <-
    structure(list(
      dist = d,
      id = id,
      k = k,
      sort = TRUE,
      metric = attr(x, "method")
    ),
      class = c("kNN", "NN"))

  return(ret)
}

#' @rdname kNN
#' @export
sort.kNN <- function(x, decreasing = FALSE, ...) {
  ## x$sort = TRUE means: sorted by increasing distance (ties by id)
  if (isTRUE(x$sort) && !decreasing)
    return(x)
  if (is.null(x$dist))
    stop("Unable to sort. Distances are missing.")
  if (ncol(x$id) < 2) {
    x$sort <- TRUE
    return(x)
  }

  ## sort first by dist and break ties using id
  o <- vapply(
    seq_len(nrow(x$dist)),
    function(i) order(x$dist[i, ], x$id[i, ], decreasing = decreasing),
    integer(ncol(x$id))
  )
  for (i in seq_len(ncol(o))) {
    x$dist[i, ] <- x$dist[i, ][o[, i]]
    x$id[i, ] <- x$id[i, ][o[, i]]
  }
  x$sort <- !decreasing

  x
}

#' @rdname kNN
#' @export
adjacencylist.kNN <- function(x, ...)
  lapply(
    seq_len(nrow(x$id)),
    FUN = function(i) {
      ## filter NAs
      tmp <- x$id[i, ]
      tmp[!is.na(tmp)]
    }
  )

#' @rdname kNN
#' @export
print.kNN <- function(x, ...) {
  cat("k-nearest neighbors for ",
    nrow(x$id),
    " objects (k=",
    x$k,
    ").",
    "\n",
    sep = "")
  cat("Distance metric:", x$metric, "\n")
  cat("\nAvailable fields: ", toString(names(x)), "\n", sep = "")
}

# Convert names to integers for C++
.parse_search <- function(search) {
  search <- pmatch(tolower(search), .ANNsearchTypes)
  if (is.na(search))
    stop("Unknown NN search method! Available for search: ",
         paste(sQuote(.ANNsearchTypes), collapse = ", "))
  search
}

.parse_shrinkRule <- function(shrinkRule) {
  shrinkRule <- pmatch(tolower(shrinkRule), .ANNshrinkRule) - 1L
  if (length(shrinkRule) != 1L || is.na(shrinkRule))
    stop("Unknown shrinkRule! Available: ",
         paste(sQuote(.ANNshrinkRule), collapse = ", "))
  shrinkRule
}

.parse_splitRule <- function(splitRule) {
  splitRule <- pmatch(tolower(splitRule), .ANNsplitRule) - 1L
  if (is.na(splitRule))
    stop("Unknown splitRule! Available: ",
         paste(sQuote(.ANNsplitRule), collapse = ", "))
  splitRule
}
