#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

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


#' Find the k Nearest Neighbors
#'
#' This function uses a kd-tree to find all k nearest neighbors in a data
#' matrix (including distances) fast.
#'
#' **Ties:** If the kth and the (k+1)th nearest neighbor are tied, then the
#' neighbor found first is returned and the other one is ignored.
#'
#' **Self-matches:** If no query is specified, then self-matches are
#' removed.
#'
#' Details on the search parameters:
#'
#' * `search` controls if
#' a kd-tree or linear search (both implemented in the ANN library; see Mount
#' and Arya, 2010). Note, that these implementations cannot handle NAs.
#' `search = "dist"` precomputes Euclidean distances first using R. NAs are
#' handled, but the resulting distance matrix cannot contain NAs. To use other
#' distance measures, a precomputed distance matrix can be provided as `x`
#' (`search` is ignored).
#'
#' * `bucketSize` and `splitRule` influence how the kd-tree is
#' built. `approx` uses the approximate nearest neighbor search
#' implemented in ANN. All nearest neighbors up to a distance of
#' `eps / (1 + approx)` will be considered and all with a distance
#' greater than `eps` will not be considered. The other points might be
#' considered. Note that this results in some actual nearest neighbors being
#' omitted leading to spurious clusters and noise points. However, the
#' algorithm will enjoy a significant speedup. For more details see Mount and
#' Arya (2010).
#'
#' @aliases kNN knn
#' @family NN functions
#'
#' @param x a data matrix, a [dist] object or a [kNN] object.
#' @param k number of neighbors to find.
#' @param query a data matrix with the points to query. If query is not
#' specified, the NN for all the points in `x` is returned. If query is
#' specified then `x` needs to be a data matrix.
#' @param search nearest neighbor search strategy (one of `"kdtree"`, `"linear"` or
#' `"dist"`).
#' @param sort sort the neighbors by distance? Note that some search methods
#' already sort the results. Sorting is expensive and `sort = FALSE` may
#' be much faster for some search methods. kNN objects can be sorted using
#' `sort()`.
#' @param bucketSize max size of the kd-tree leafs.
#' @param splitRule rule to split the kd-tree. One of `"STD"`, `"MIDPT"`, `"FAIR"`,
#' `"SL_MIDPT"`, `"SL_FAIR"` or `"SUGGEST"` (SL stands for sliding). `"SUGGEST"` uses
#' ANNs best guess.
#' @param approx use approximate nearest neighbors. All NN up to a distance of
#' a factor of `1 + approx` eps may be used. Some actual NN may be omitted
#' leading to spurious clusters and noise points.  However, the algorithm will
#' enjoy a significant speedup.
#' @param decreasing sort in decreasing order?
#' @param ... further arguments
#'
#' @return An object of class `kNN` (subclass of [NN]) containing a
#' list with the following components:
#' \item{dist }{a matrix with distances. }
#' \item{id }{a matrix with `ids`. }
#' \item{k }{number `k` used. }
#' \item{metric }{ used distance metric. }
#'
#' @author Michael Hahsler
#' @references David M. Mount and Sunil Arya (2010). ANN: A Library for
#' Approximate Nearest Neighbor Searching,
#' \url{http://www.cs.umd.edu/~mount/ANN/}.
#' @keywords model
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # Example 1: finding kNN for all points in a data matrix (using a kd-tree)
#' nn <- kNN(x, k = 5)
#' nn
#'
#' # explore neighborhood of point 10
#' i <- 10
#' nn$id[i,]
#' plot(x, col = ifelse(1:nrow(iris) %in% nn$id[i,], "red", "black"))
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
#' plot(nn, x, col = "grey")
#' points(q, pch = 3, lwd = 2)
#'
#' # Example 3: find kNN using distances
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
    approx = 0) {
    if (inherits(x, "kNN")) {
      if (x$k < k)
        stop("kNN in x has not enough nearest neighbors.")
      if (!x$sort)
        x <- sort(x)
      x$id <- x$id[, 1:k]
      if (!is.null(x$dist))
        x$dist <- x$dist[, 1:k]
      if (!is.null(x$shared))
        x$dist <- x$shared[, 1:k]
      x$k <- k
      return(x)
    }

    search <- .parse_search(search)
    splitRule <- .parse_splitRule(splitRule)

    k <- as.integer(k)
    if (k < 1)
      stop("Illegal k: needs to be k>=1!")

    ### dist search
    if (search == 3) {
      if (!inherits(x, "dist"))
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
    if (!.matrixlike(x))
      stop("x needs to be a matrix to calculate distances")
    x <- as.matrix(x)
    if (storage.mode(x) == "integer")
      storage.mode(x) <- "double"
    if (storage.mode(x) != "double")
      stop("x has to be a numeric matrix.")

    if (!is.null(query)) {
      query <- as.matrix(query)
      if (storage.mode(query) == "integer")
        storage.mode(query) <- "double"
      if (storage.mode(query) != "double")
        stop("query has to be NULL or a numeric matrix.")
      if (ncol(x) != ncol(query))
        stop("x and query need to have the same number of columns!")
    }

    if (k >= nrow(x))
      stop("Not enough neighbors in data set!")


    if (anyNA(x))
      stop("data/distances cannot contain NAs for kNN (with kd-tree)!")

    ## returns NO self matches
    if (!is.null(query)) {
      ret <- kNN_query_int(
        as.matrix(x),
        as.matrix(query),
        as.integer(k),
        as.integer(search),
        as.integer(bucketSize),
        as.integer(splitRule),
        as.double(approx)
      )
      dimnames(ret$dist) <- list(rownames(query), 1:k)
      dimnames(ret$id) <- list(rownames(query), 1:k)
    } else{
      ret <- kNN_int(
        as.matrix(x),
        as.integer(k),
        as.integer(search),
        as.integer(bucketSize),
        as.integer(splitRule),
        as.double(approx)
      )
      dimnames(ret$dist) <- list(rownames(x), 1:k)
      dimnames(ret$id) <- list(rownames(x), 1:k)
    }

    class(ret) <- c("kNN", "NN")

    ### ANN already returns them sorted (by dist but not by ID)
    if (sort)
      ret <- sort(ret)

    ret$metric = "euclidean"

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
  if(attr(x, "Diag") || attr(x, "Upper"))
    stop("x needs to be a dist object with attributes Diag and Upper set to FALSE. Use as.dist(x, diag = FALSE, upper = FALSE) fist.")
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
  if (isTRUE(x$sort))
    return(x)
  if (is.null(x$dist))
    stop("Unable to sort. Distances are missing.")
  if (ncol(x$id) < 2) {
    x$sort <- TRUE
    return(x)
  }

  ## sort first by dist and break ties using id
  o <- sapply(
    1:nrow(x$dist),
    FUN =
      function(i)
        order(x$dist[i,], x$id[i,], decreasing = decreasing)
  )
  for (i in 1:ncol(o)) {
    x$dist[i,] <- x$dist[i,][o[, i]]
    x$id[i,] <- x$id[i,][o[, i]]
  }
  x$sort <- TRUE

  x
}

#' @rdname kNN
#' @export
adjacencylist.kNN <- function(x, ...)
  lapply(
    seq(nrow(x$id)),
    FUN = function(i) {
      ## filter NAs
      tmp <- x$id[i,]
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
  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if (is.na(search))
    stop("Unknown NN search type!")
  search
}

.parse_splitRule <- function(splitRule) {
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule) - 1L
  if (is.na(splitRule))
    stop("Unknown splitRule!")
  splitRule
}
