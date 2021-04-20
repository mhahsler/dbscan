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

lof <- function(x, k, ...) {

  extra <- list(...)
  args <- c("search", "bucketSize", "splitRule", "approx")
  m <- pmatch(names(extra), args)
  if(any(is.na(m))) stop("Unknown parameter: ",
    paste(names(extra)[is.na(m)], collapse = ", "))
  names(extra) <- args[m]

  search <- if(is.null(extra$search)) "kdtree" else extra$search
  search <- .parse_search(search)
  splitRule <- if(is.null(extra$splitRule)) "suggest" else extra$splitRule
  splitRule <- .parse_splitRule(splitRule)
  bucketSize <- if(is.null(extra$bucketSize)) 10L else
    as.integer(extra$bucketSize)
  approx <- if(is.null(extra$approx)) 0L else as.integer(extra$approx)

  ### dist search
  if(search == 3) {
    if(!inherits(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  # get n
  if(inherits(x, "dist")) n <- attr(x, "Size")
  else n <- nrow(x)

  if(is.null(n)) stop("x needs to be a matrix or a dist object!")

  if(k<1 || k>=n)
    stop("k has to be larger than 1 and smaller than the number of points")


  ### get LOF from a dist object
  if(inherits(x, "dist")) {

    if(any(is.na(x))) stop("NAs not allowed in dist for LOF!")

    # find k-NN distance, ids and distances
    x <- as.matrix(x)
    diag(x) <- Inf ### no self-matches
    o <- t(apply(x, 1, order, decreasing = FALSE))
    k_dist <- x[cbind(o[,k], seq_len(n))]
    ids <- lapply(seq_len(n), FUN = function(i) which(x[i, ] <= k_dist[i]))
    dist <- lapply(seq_len(n), FUN = function(i) x[i, x[i, ] <= k_dist[i]])

    ret <- list(dist = dist, ids = ids, k_dist = k_dist)

  }else{

    if(any(is.na(x))) stop("NAs not allowed for LOF using kdtree!")
    ### Use kd-tree

    ret <- lof_kNN(as.matrix(x), as.integer(k),
      as.integer(search), as.integer(bucketSize),
      as.integer(splitRule), as.double(approx))
  }

  # calculate local reachability density (LRD)
  # reachability-distance_k(A,B) = max{k-distance(B), d(A,B)}
  # lrdk(A) = 1/(sum_B \in N_k(A) reachability-distance_k(A, B) / |N_k(A)|)
  lrd <- numeric(n)
  for(A in seq_len(n)) {
    Bs <- ret$ids[[A]]
    lrd[A] <- 1/(sum(pmax.int(ret$k_dist[Bs], ret$dist[[A]])) / length(Bs))
  }

  # calculate local outlier factor (LOF)
  # LOF_k(A) = sum_B \in N_k(A) lrd_k(B)/(|N_k(A)| lrdk(A))
  lof <- numeric(n)
  for (A in seq_len(n)) {
    Bs <- ret$ids[[A]]
    lof[A] <- sum(lrd[Bs])/ length(Bs) / lrd[A]
  }

  # with more than k duplicates lrd can become infinity
  # we define them not to be outliers
  lof[is.nan(lof)] <- 1

  lof
}
