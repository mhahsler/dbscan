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

#' Calculate and Plot k-Nearest Neighbor Distances
#'
#' Fast calculation of the k-nearest neighbor distances for a dataset
#' represented as a matrix of points. The kNN distance is defined as the
#' distance from a point to its k nearest neighbor. The kNN distance plot
#' displays the kNN distance of all points sorted from smallest to largest. The
#' plot can be used to help find suitable parameter values for the for DBSCAN.
#'
#' @family Outlier Detection Functions
#' @family NN functions
#'
#' @param x the data set as a matrix of points (Euclidean distance is used) or
#' a precalculated [dist] object.
#' @param k number of nearest neighbors used for the distance calculation.
#' @param all should a matrix with the distances to all k nearest neighbors be
#' returned?
#' @param ... further arguments (e.g., kd-tree related parameters) are passed
#' on to [kNN()].
#'
#' @return `kNNdist()` returns a numeric vector with the distance to its k
#' nearest neighbor. If `all = TRUE` then a matrix with k columns
#' containing the distances to all 1st, 2nd, ..., kth nearest neighbors is
#' returned instead.
#'
#' @author Michael Hahsler
#' @keywords model plot
#' @examples
#' data(iris)
#' iris <- as.matrix(iris[, 1:4])
#'
#' ## Find the 4-NN distance for each observation (see ?kNN
#' ## for different search strategies)
#' kNNdist(iris, k = 4)
#'
#' ## Get a matrix with distances to the 1st, 2nd, ..., 4th NN.
#' kNNdist(iris, k = 4, all = TRUE)
#'
#' ## Produce a k-NN distance plot to determine a suitable eps for
#' ## DBSCAN with MinPts = 5. Use k = 4 (= MinPts -1).
#' ## The knee is visible around a distance of .7
#' kNNdistplot(iris, k = 4)
#'
#' cl <- dbscan(iris, eps = .7, minPts = 5)
#' pairs(iris, col = cl$cluster + 1L)
#' ## Note: black points are noise points
#' @export kNNdist
kNNdist <- function(x, k, all = FALSE, ...) {
  kNNd <- dbscan::kNN(x, k, sort = TRUE, ...)$dist
  if (!all)
    kNNd <- kNNd[, k]
  kNNd
}

#' @rdname kNNdist
kNNdistplot <- function(x, k, ...) {
  kNNdist <- sort(kNNdist(x, k , ...))
  plot(
    sort(kNNdist),
    type = "l",
    ylab = paste(k, "-NN distance", sep = ""),
    xlab = "Points (sample) sorted by distance"
  )
}
