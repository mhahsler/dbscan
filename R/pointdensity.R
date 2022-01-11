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

#' Calculate Local Density at Each Data Point
#'
#' Calculate the local density at each data point as either the number of
#' points in the eps-neighborhood (as used in `dbscan()`) or perform kernel density
#' estimation (KDE) using a uniform kernel. The function uses a kd-tree for fast
#' fixed-radius nearest neighbor search.
#'
#' `dbscan()` estimates the density around a point as the number of points in the
#' eps-neighborhood of the point (including the query point itself).
#' Kernel density estimation (KDE) using a uniform kernel, which is just this point
#' count in the eps-neighborhood divided by \eqn{(2\,eps\,n)}{(2 eps n)}, where
#' \eqn{n} is the number of points in `x`.
#'
#' Points with low local density often indicate noise (see e.g., Wishart (1969)
#' and Hartigan (1975)).
#'
#' @aliases pointdensity density
#' @family Outlier Detection Functions
#'
#' @param x a data matrix.
#' @param eps radius of the eps-neighborhood, i.e., bandwidth of the uniform
#' kernel).
#' @param type `"frequency"` or `"density"`. should the raw count of
#' points inside the eps-neighborhood or the kde be returned.
#' @param search,bucketSize,splitRule,approx algorithmic parameters for
#' [frNN()].
#'
#' @return A vector of the same length as data points (rows) in `x` with
#' the count or density values for each data point.
#'
#' @author Michael Hahsler
#' @seealso [frNN()], [stats::density()].
#' @references Wishart, D. (1969), Mode Analysis: A Generalization of Nearest
#' Neighbor which Reduces Chaining Effects, in _Numerical Taxonomy,_ Ed., A.J.
#' Cole, Academic Press, 282-311.
#'
#' John A. Hartigan (1975), _Clustering Algorithms,_ John Wiley & Sons, Inc.,
#' New York, NY, USA.
#' @keywords model
#' @examples
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
#'   y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
#'   )
#' plot(x)
#'
#' ### calculate density
#' d <- pointdensity(x, eps = .5, type = "density")
#'
#' ### density distribution
#' summary(d)
#' hist(d, breaks = 10)
#'
#' ### plot with point size is proportional to Density
#' plot(x, pch = 19, main = "Density (eps = .5)", cex = d*5)
#'
#' ### Wishart (1969) single link clustering after removing low-density noise
#' # 1. remove noise with low density
#' f <- pointdensity(x, eps = .5, type = "frequency")
#' x_nonoise <- x[f >= 5,]
#'
#' # 2. use single-linkage on the non-noise points
#' hc <- hclust(dist(x_nonoise), method = "single")
#' plot(x, pch = 19, cex = .5)
#' points(x_nonoise, pch = 19, col= cutree(hc, k = 4) + 1L)
#' @export pointdensity
pointdensity <- function(x,
  eps,
  type = "frequency",
  search = "kdtree",
  bucketSize = 10,
  splitRule = "suggest",
  approx = 0) {
  type <- match.arg(type, choices = c("frequency", "density"))

  search <- .parse_search(search)
  splitRule <- .parse_splitRule(splitRule)

  d <- dbscan_density_int(
    as.matrix(x),
    as.double(eps),
    as.integer(search),
    as.integer(bucketSize),
    as.integer(splitRule),
    as.double(approx)
  )

  if (type == "density")
    d <- d / (2 * eps * nrow(x))

  d
}

#gof <- function(x, eps, ...) {
#  d <- pointdensity(x, eps, ...)
#  1/(d/mean(d))
#}
