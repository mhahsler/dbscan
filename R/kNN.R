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


# Convert names to integers for C++
.parse_search <- function(search) {
  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if(is.na(search)) stop("Unknown NN search type!")
  search
}

.parse_splitRule <- function(splitRule) {
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")
  splitRule
}

kNN <- function(x, k, query = NULL, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0) {

  if(inherits(x, "kNN")) {
    if(x$k < k) stop("kNN in x has not enough nearest neighbors.")
    if(!x$sort) x <- sort(x)
    x$id <- x$id[,1:k]
    if(!is.null(x$dist)) x$dist <- x$dist[,1:k]
    if(!is.null(x$shared)) x$dist <- x$shared[,1:k]
    x$k <- k
    return(x)
  }

  search <- .parse_search(search)
  splitRule <- .parse_splitRule(splitRule)

  k <- as.integer(k)
  if(k < 1) stop("Illegal k: needs to be k>=1!")

  ### dist search
  if(search == 3) {
    if(!inherits(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  ### get kNN from a dist object
  if(inherits(x, "dist")) {

    if(!is.null(query)) stop("query can only be used if x contains the data.")

    if(any(is.na(x))) stop("distances cannot be NAs for kNN!")

    x <- as.matrix(x)
    diag(x) <- NA

    if(k >= nrow(x)) stop("Not enough neighbors in data set!")

    o <- t(apply(x, 1, order, decreasing = FALSE))
    id <- o[,1:k, drop = FALSE]
    dimnames(id) <- list(rownames(x), 1:k)
    d <- t(sapply(1:nrow(id), FUN = function(i) x[i, id[i,], drop = FALSE]))
    if(k==1) d <- t(d) ### sapply drops an array with a single column!
    dimnames(d) <- list(rownames(x), 1:k)

    return(structure(list(dist = d, id = id, k = k, sort = TRUE),
      class = c("kNN", "NN")))
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

  if(k >= nrow(x)) stop("Not enough neighbors in data set!")


  if(any(is.na(x))) stop("data/distances cannot contain NAs for kNN (with kd-tree)!")

  ## returns NO self matches
  if(!is.null(query)) {
    ret <- kNN_query_int(as.matrix(x), as.matrix(query), as.integer(k),
      as.integer(search), as.integer(bucketSize),
      as.integer(splitRule), as.double(approx))
    dimnames(ret$dist) <- list(rownames(query), 1:k)
    dimnames(ret$id) <- list(rownames(query), 1:k)
  }else{
    ret <- kNN_int(as.matrix(x), as.integer(k),
      as.integer(search), as.integer(bucketSize),
      as.integer(splitRule), as.double(approx))
    dimnames(ret$dist) <- list(rownames(x), 1:k)
    dimnames(ret$id) <- list(rownames(x), 1:k)
  }

  class(ret) <- c("kNN", "NN")

  ### ANN already returns them sorted (by dist but not by ID)
  if(sort) ret <- sort(ret)

  ret
}


### FIXME: This is expensive! We should do this in C++
sort.kNN <- function(x, decreasing = FALSE, ...) {
  if(!is.null(x$sort) && x$sort) return(x)
  if(is.null(x$dist)) stop("Unable to sort. Distances are missing.")
  if(ncol(x$id)<2) {
    x$sort <- TRUE
    return(x)
  }

  ## sort first by dist and break ties using id
  o <- sapply(1:nrow(x$dist), FUN =
      function(i) order(x$dist[i,], x$id[i,], decreasing=decreasing))
  for(i in 1:ncol(o)) {
    x$dist[i,] <- x$dist[i,][o[,i]]
    x$id[i,] <- x$id[i,][o[,i]]
  }
  x$sort <- TRUE

  x
}

print.kNN <- function(x, ...) {
  cat("k-nearest neighbors for ", nrow(x$id), " objects (k=", x$k,").",
    "\n", sep = "")
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}

