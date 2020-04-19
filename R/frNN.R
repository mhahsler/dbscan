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

frNN <- function(x, eps, query = NULL, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0) {

  if(is.null(eps) || is.na(eps) || eps < 0) stop("eps needs to be >=0.")

  if(inherits(x, "frNN")) {
    if(x$eps < eps) stop("frNN in x has not a sufficient eps radius.")

    for(i in 1:length(x$dist)) {
      take <- x$dist[[i]] <= eps
      x$dist[[i]] <- x$dist[[i]][take]
      x$id[[i]] <- x$id[[i]][take]
    }
    x$eps <- eps

    return(x)
  }

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if(is.na(search)) stop("Unknown NN search type!")

  ### dist search
  if(search == 3) {
    if(!inherits(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  ### get kNN from a dist object in R
  if(inherits(x, "dist")) {

    if(!is.null(query)) stop("query can only be used if x contains the data.")

    if(any(is.na(x))) stop("data/distances cannot contain NAs for frNN (with kd-tree)!")

    x <- as.matrix(x)
    diag(x) <- Inf

    id <- lapply(1:nrow(x), FUN = function(i) {
          y <- x[i, ]
          o <- order(y, decreasing = FALSE)
          o[y[o] <= eps]

      }
    )
    names(id) <- rownames(x)

    d <- lapply(1:nrow(x), FUN = function(i) {
          unname(x[i,id[[i]]])
        }
    )
    names(d) <- rownames(x)

    ret <- structure(list(dist = d, id = id, eps = eps, sort = TRUE),
      class = c("frNN", "NN"))

    return(ret)
  }

  ## make sure x is numeric
  if(!.matrixlike(x)) stop("x needs to be a matrix to calculate distances")
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")

  if(!is.null(query)) {
    query <- as.matrix(query)
    if(storage.mode(query) == "integer") storage.mode(query) <- "double"
    if(storage.mode(query) != "double") stop("query has to be NULL or a numeric matrix.")
    if(ncol(x) != ncol(query)) stop("x and query need to have the same number of columns!")
  }

  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  if(any(is.na(x))) stop("data/distances cannot contain NAs for frNN (with kd-tree)!")

  ## returns NO self matches
  if(!is.null(query)) {
    ret <- frNN_query_int(as.matrix(x), as.matrix(query), as.double(eps),
      as.integer(search), as.integer(bucketSize),
      as.integer(splitRule), as.double(approx))
    names(ret$dist) <- rownames(query)
    names(ret$id) <- rownames(query)
  }else{
    ret <- frNN_int(as.matrix(x), as.double(eps),
      as.integer(search), as.integer(bucketSize),
      as.integer(splitRule), as.double(approx))
    names(ret$dist) <- rownames(x)
    names(ret$id) <- rownames(x)
  }

  ret$eps <- eps
  ret$sort <- FALSE
  class(ret) <- c("frNN", "NN")

  if(sort) ret <- sort.frNN(ret)

  ret
}

sort.frNN <- function(x, decreasing = FALSE, ...) {
  if(!is.null(x$sort) && x$sort) return(x)
  if(is.null(x$dist)) stop("Unable to sort. Distances are missing.")

  ## FIXME: This is slow do this in C++
  n <- names(x$id)

  o <- lapply(1:length(x$dist), FUN =
      function(i) order(x$dist[[i]], x$id[[i]], decreasing=decreasing))
  x$dist <- lapply(1:length(o), FUN = function(p) x$dist[[p]][o[[p]]])
  x$id <- lapply(1:length(o), FUN = function(p) x$id[[p]][o[[p]]])

  names(x$dist) <- n
  names(x$id) <- n

  x$sort <- TRUE

  x
}

print.frNN <- function(x, ...) {
  cat("fixed radius nearest neighbors for ", length(x$id),
    " objects (eps=", x$eps,").", "\n", sep = "")
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}
