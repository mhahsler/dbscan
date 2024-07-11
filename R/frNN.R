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


#' Find the Fixed Radius Nearest Neighbors
#'
#' This function uses a kd-tree to find the fixed radius nearest neighbors
#' (including distances) fast.
#'
#' If `x` is specified as a data matrix, then Euclidean distances an fast
#' nearest neighbor lookup using a kd-tree are used.
#'
#' To create a frNN object from scratch, you need to supply at least the
#' elements `id` with a list of integer vectors with the nearest neighbor
#' ids for each point and `eps` (see below).
#'
#' **Self-matches:** Self-matches are not returned!
#'
#' @aliases frNN frnn print.frnn
#' @family NN functions
#'
#' @param x a data matrix, a dist object or a frNN object.
#' @param eps neighbors radius.
#' @param query a data matrix with the points to query. If query is not
#' specified, the NN for all the points in `x` is returned. If query is
#' specified then `x` needs to be a data matrix.
#' @param sort sort the neighbors by distance? This is expensive and can be
#' done later using `sort()`.
#' @param search nearest neighbor search strategy (one of `"kdtree"`, `"linear"` or
#' `"dist"`).
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
#' @returns
#'
#' `frNN()` returns an object of class [frNN] (subclass of
#' [NN]) containing a list with the following components:
#' \item{id }{a list of
#' integer vectors. Each vector contains the ids of the fixed radius nearest
#' neighbors. }
#' \item{dist }{a list with distances (same structure as
#' `id`). }
#' \item{eps }{ neighborhood radius `eps` that was used. }
#' \item{metric }{ used distance metric. }
#'
#' `adjacencylist()` returns a list with one entry per data point in `x`. Each entry
#' contains the id of the nearest neighbors.
#'
#' @author Michael Hahsler
#'
#' @references David M. Mount and Sunil Arya (2010). ANN: A Library for
#' Approximate Nearest Neighbor Searching,
#' \url{http://www.cs.umd.edu/~mount/ANN/}.
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
#'   xlab = "k", main="Number of Neighbors",
#'   sub = paste("Neighborhood size eps =", nn$eps))
#'
#' # Explore neighbors of point i = 10
#' i <- 10
#' nn$id[[i]]
#' nn$dist[[i]]
#' plot(x, col = ifelse(1:nrow(iris) %in% nn$id[[i]], "red", "black"))
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
#' plot(nn, x, col = "grey")
#' points(q, pch = 3, lwd = 2)
#' @export frNN
frNN <-
  function(x,
    eps,
    query = NULL,
    sort = TRUE,
    search = "kdtree",
    bucketSize = 10,
    splitRule = "suggest",
    approx = 0) {
    if (is.null(eps) ||
        is.na(eps) || eps < 0)
      stop("eps needs to be >=0.")

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

    ### dist search
    if (search == 3) {
      if (!inherits(x, "dist"))
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
        stop("data/distances cannot contain NAs for frNN (with kd-tree)!")

      return(dist_to_frNN(x, eps = eps, sort = sort))
    }

    ## make sure x is numeric
    if (!.matrixlike(x))
      stop("x needs to be a matrix or a data.frame.")
    x <- as.matrix(x)
    if (storage.mode(x) == "integer")
      storage.mode(x) <- "double"
    if (storage.mode(x) != "double")
      stop("all data in x has to be numeric.")

    if (!is.null(query)) {
      if (!.matrixlike(query))
        stop("query needs to be a matrix or a data.frame.")
      query <- as.matrix(query)
      if (storage.mode(query) == "integer")
        storage.mode(query) <- "double"
      if (storage.mode(query) != "double")
        stop("query has to be NULL or a numeric matrix or data.frame.")
      if (ncol(x) != ncol(query))
        stop("x and query need to have the same number of columns!")
    }

    if (anyNA(x))
      stop("data/distances cannot contain NAs for frNN (with kd-tree)!")

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
  if (isTRUE(x$sort))
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

  x$sort <- TRUE

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
