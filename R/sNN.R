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
sNN <- function(x, k, kt = NULL, jp = FALSE, sort = TRUE,
  search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0){

  if(missing(k)) k <- x$k

  if(inherits(x, "kNN")) {
    if(k != x$k) {
      if(ncol(x$id) < k) stop("kNN object does not contain enough neighbors!")
      if(!x$sort) x <- sort.kNN(x)
      x$id <- x$id[,1:k]
      x$dist <- x$dist[,1:k]
      x$k <- k
    }

  } else x <- kNN(x, k, sort = FALSE, search = search, bucketSize = bucketSize,
    splitRule = splitRule, approx = approx)

  x$shared <- SNN_sim_int(x$id, as.logical(jp[1]))
  x$sort_shared <- FALSE

  class(x) <- c("sNN", "kNN", "NN")

  if(sort) x <- sort.sNN(x)

  x$kt <- kt

  if(!is.null(kt)) {
    if(kt > k) stop("kt needs to be less than k.")
    rem <- x$shared < kt
    x$id[rem] <- NA
    x$dist[rem] <- NA
    x$shared[rem] <- NA
  }

  x
}


sort.sNN <- function(x, decreasing = TRUE, ...) {
  if(!is.null(x$sort_shared) && x$sort_shared) return(x)
  if(is.null(x$shared)) stop("Unable to sort. Number of shared neighbors is missing.")
  if(ncol(x$id)<2) {
    x$sort <- TRUE
    x$sort_shared <- TRUE
    return(x)
  }

  ## sort first by number of shared points (decreasing) and break ties by id (increasing)
  k <- ncol(x$shared)
  o <- sapply(1:nrow(x$shared), FUN =
      function(i) order(k-x$shared[i,], x$id[i,], decreasing=!decreasing))
  for(i in 1:ncol(o)) {
    x$shared[i,] <- x$shared[i,][o[,i]]
    x$dist[i,] <- x$dist[i,][o[,i]]
    x$id[i,] <- x$id[i,][o[,i]]
  }

  x$sort <- FALSE
  x$sort_shared <- TRUE

  x
}

print.sNN <- function(x, ...) {
  cat("shared-nearest neighbors for ", nrow(x$id), " objects (k=", x$k,
    ", kt=", ifelse(is.null(x$kt), "NULL", x$kt), ").", "\n", sep = "")
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}
