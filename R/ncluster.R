#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler, Matt Piekenbrock

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

#' Number of Clusters, Noise Points, and Observations
#'
#' Extract the number of clusters or the number of noise points for
#' a clustering. This function works with any clustering result that
#' contains a list element named `cluster` with a clustering vector. In
#' addition, `nobs` (see [stats::nobs()]) is also available to retrieve
#' the number of clustered points.
#'
#' @name ncluster
#' @aliases ncluster nnoise nobs
#' @family clustering functions
#'
#' @param object a clustering result object containing a `cluster` element.
#' @param ...  additional arguments are unused.
#'
#' @return returns the number if clusters or noise points.
#' @examples
#' data(iris)
#' iris <- as.matrix(iris[, 1:4])
#'
#' res <- dbscan(iris, eps = .7, minPts = 5)
#' res
#'
#' ncluster(res)
#' nnoise(res)
#' nobs(res)
#'
#' # the functions also work with kmeans and other clustering algorithms.
#' cl <- kmeans(iris, centers = 3)
#' ncluster(cl)
#' nnoise(cl)
#' nobs(res)
#' @export
ncluster <- function(object, ...) {
  UseMethod("ncluster")
}

#' @export
ncluster.default <- function(object, ...) {
  if (!is.list(object) || !is.numeric(object$cluster))
    stop("ncluster() requires a clustering object with a cluster component containing the cluster labels.")

  length(setdiff(unique(object$cluster), 0L))
}

#' @rdname ncluster
#' @export
nnoise <- function(object, ...) {
  UseMethod("nnoise")
}

#' @export
nnoise.default <- function(object, ...) {
  if (!is.list(object) || !is.numeric(object$cluster))
    stop("ncluster() requires a clustering object with a cluster component containing the cluster labels.")

  sum(object$cluster == 0L)
}
