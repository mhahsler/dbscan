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


#' Local Outlier Factor Score
#'
#' Calculate the Local Outlier Factor (LOF) score for each data point using a
#' kd-tree to speed up kNN search.
#'
#' LOF compares the local readability density (lrd) of an point to the lrd of
#' its neighbors. A LOF score of approximately 1 indicates that the lrd around
#' the point is comparable to the lrd of its neighbors and that the point is
#' not an outlier. Points that have a substantially lower lrd than their
#' neighbors are considered outliers and produce scores significantly larger
#' than 1.
#'
#' If a data matrix is specified, then Euclidean distances and fast nearest
#' neighbor search using a kd-tree is used.
#'
#' **Note on duplicate points:** If there are more than `minPts`
#' duplicates of a point in the data, then LOF the local readability distance
#' will be 0 resulting in an undefined LOF score of 0/0. We set LOF in this
#' case to 1 since there is already enough density from the points in the same
#' location to make them not outliers. The original paper by Breunig et al
#' (2000) assumes that the points are real duplicates and suggests to remove
#' the duplicates before computing LOF. If duplicate points are removed first,
#' then this LOF implementation in \pkg{dbscan} behaves like the one described
#' by Breunig et al.
#'
#' @aliases lof LOF
#' @family Outlier Detection Functions
#'
#' @param x a data matrix or a [dist] object.
#' @param minPts number of nearest neighbors used in defining the local
#' neighborhood of a point (includes the point itself).
#' @param ... further arguments are passed on to [kNN()].
#' Note: `sort` cannot be specified here since `lof()`
#' uses always `sort = TRUE`.
#'
#' @return A numeric vector of length `ncol(x)` containing LOF values for
#' all data points.
#'
#' @author Michael Hahsler
#' @references Breunig, M., Kriegel, H., Ng, R., and Sander, J. (2000). LOF:
#' identifying density-based local outliers. In _ACM Int. Conf. on
#' Management of Data,_ pages 93-104.
#' \doi{10.1145/335191.335388}
#' @keywords model
#' @examples
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
#'   y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
#'   )
#'
#' ### calculate LOF score with a neighborhood of 3 points
#' lof <- lof(x, minPts = 3)
#'
#' ### distribution of outlier factors
#' summary(lof)
#' hist(lof, breaks = 10, main = "LOF (minPts = 3)")
#'
#' ### plot sorted lof. Looks like outliers start arounf a LOF of 2.
#' plot(sort(lof), type = "l",  main = "LOF (minPts = 3)",
#'   xlab = "Points sorted by LOF", ylab = "LOF")
#'
#' ### point size is proportional to LOF and mark points with a LOF > 2
#' plot(x, pch = ".", main = "LOF (minPts = 3)", asp = 1)
#' points(x, cex = (lof - 1) * 2, pch = 1, col = "red")
#' text(x[lof > 2,], labels = round(lof, 1)[lof > 2], pos = 3)
#' @export
lof <- function(x, minPts = 5, ...) {
  ### parse extra parameters
  extra <- list(...)

  # check for deprecated k
  if (!is.null(extra[["k"]])) {
    minPts <- extra[["k"]] + 1
    extra[["k"]] <- NULL
    warning("lof: k is now deprecated. use minPts = ", minPts, " instead .")
  }

  args <- c("search", "bucketSize", "splitRule", "approx")
  m <- pmatch(names(extra), args)
  if (anyNA(m))
    stop("Unknown parameter: ",
      toString(names(extra)[is.na(m)]))
  names(extra) <- args[m]

  search <- if (is.null(extra$search))
    "kdtree"
  else
    extra$search
  search <- .parse_search(search)
  splitRule <-
    if (is.null(extra$splitRule))
      "suggest"
  else
    extra$splitRule
  splitRule <- .parse_splitRule(splitRule)
  bucketSize <- if (is.null(extra$bucketSize))
    10L
  else
    as.integer(extra$bucketSize)
  approx <- if (is.null(extra$approx))
    0
  else
    as.double(extra$approx)

  ### precompute distance matrix for dist search
  if (search == 3) {
    if (!inherits(x, "dist"))
      if (.matrixlike(x))
        x <- dist(x)
    else
      stop("x needs to be a matrix to calculate distances")
  }

  # get and check n
  if (inherits(x, "dist"))
    n <- attr(x, "Size")
  else
    n <- nrow(x)
  if (is.null(n))
    stop("x needs to be a matrix or a dist object!")
  if (minPts < 2 || minPts > n)
    stop("minPts has to be at least 2 and not larger than the number of points")


  ### get LOF from a dist object
  if (inherits(x, "dist")) {
    if (anyNA(x))
      stop("NAs not allowed in dist for LOF!")

    # find k-NN distance, ids and distances
    x <- as.matrix(x)
    diag(x) <- Inf ### no self-matches
    o <- t(apply(x, 1, order, decreasing = FALSE))
    k_dist <- x[cbind(o[, minPts - 1], seq_len(n))]
    ids <-
      lapply(
        seq_len(n),
        FUN = function(i)
          which(x[i,] <= k_dist[i])
      )
    dist <-
      lapply(
        seq_len(n),
        FUN = function(i)
          x[i, x[i,] <= k_dist[i]]
      )

    ret <- list(k_dist = k_dist,
      ids = ids,
      dist = dist)

  } else{
    ### Use kd-tree

    if (anyNA(x))
      stop("NAs not allowed for LOF using kdtree!")

    ret <- lof_kNN(
      as.matrix(x),
      as.integer(minPts),
      as.integer(search),
      as.integer(bucketSize),
      as.integer(splitRule),
      as.double(approx)
    )
  }

  # calculate local reachability density (LRD)
  # reachability-distance_k(A,B) = max{k-distance(B), d(A,B)}
  # lrdk(A) = 1/(sum_B \in N_k(A) reachability-distance_k(A, B) / |N_k(A)|)
  lrd <- numeric(n)
  for (A in seq_len(n)) {
    Bs <- ret$ids[[A]]
    lrd[A] <-
      1 / (sum(pmax.int(ret$k_dist[Bs], ret$dist[[A]])) / length(Bs))
  }

  # calculate local outlier factor (LOF)
  # LOF_k(A) = sum_B \in N_k(A) lrd_k(B)/(|N_k(A)| lrdk(A))
  lof <- numeric(n)
  for (A in seq_len(n)) {
    Bs <- ret$ids[[A]]
    lof[A] <- sum(lrd[Bs]) / length(Bs) / lrd[A]
  }

  # with more than k duplicates lrd can become infinity
  # we define them not to be outliers
  lof[is.nan(lof)] <- 1

  lof
}
