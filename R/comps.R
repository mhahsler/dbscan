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

#' Find Connected Components in a Nearest-neighbor Graph
#'
#' Generic function and methods to find connected components in nearest neighbor graphs.
#'
#' Note that for kNN graphs, one point may be in the kNN of the other but nor vice versa.
#' `mutual = TRUE` requires that both points are in each other's kNN.
#'
#' @family NN functions
#' @aliases components
#'
#' @param x the [NN] object representing the graph or a [dist] object
#' @param eps threshold on the distance
#' @param mutual for a pair of points, do both have to be in each other's neighborhood?
#' @param ... further arguments are currently unused.
#'
#' @return an integer vector with component assignments.
#'
#' @author Michael Hahsler
#' @keywords model
#' @examples
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
#'   y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
#'   )
#' plot(x, pch = 16)
#'
#' # Connected components on a graph where each pair of points
#' # with a distance less or equal to eps are connected
#' d <- dist(x)
#' components <- comps(d, eps = .8)
#' plot(x, col = components, pch = 16)
#'
#' # Connected components in a fixed radius nearest neighbor graph
#' # Gives the same result as the threshold on the distances above
#' frnn <- frNN(x, eps = .8)
#' components <- comps(frnn)
#' plot(frnn, data = x, col = components)
#'
#' # Connected components on a k nearest neighbors graph
#' knn <- kNN(x, 3)
#' components <- comps(knn, mutual = FALSE)
#' plot(knn, data = x, col = components)
#'
#' components <- comps(knn, mutual = TRUE)
#' plot(knn, data = x, col = components)
#'
#' # Connected components in a shared nearest neighbor graph
#' snn <- sNN(x, k = 10, kt = 5)
#' components <- comps(snn)
#' plot(snn, data = x, col = components)
#' @export
comps <- function(x, ...) UseMethod("comps", x)

#' @rdname comps
#' @export
comps.dist <- function(x, eps, ...)
  stats::cutree(stats::hclust(x, method = "single"), h = eps)

#' @rdname comps
#' @export
comps.kNN <- function(x, mutual = FALSE, ...)
  as.integer(factor(comps_kNN(x$id, as.logical(mutual))))

# sNN and frNN are symmetric so no need for mutual
#' @rdname comps
#' @export
comps.sNN <- function(x, ...) comps.kNN(x, mutual = FALSE)

#' @rdname comps
#' @export
comps.frNN <- function(x, ...) comps_frNN(x$id, mutual = FALSE)
