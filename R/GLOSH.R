#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler, Matthew Piekenbrock

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

#' Global-Local Outlier Score from Hierarchies
#'
#' Calculate the Global-Local Outlier Score from Hierarchies (GLOSH) score for
#' each data point using a kd-tree to speed up kNN search.
#'
#' GLOSH compares the density of a point to densities of any points associated
#' within current and child clusters (if any). Points that have a substantially
#' lower density than the density mode (cluster) they most associate with are
#' considered outliers. GLOSH is computed from a hierarchy a clusters.
#'
#' Specifically, consider a point \emph{x} and a density or distance threshold
#' \emph{lambda}. GLOSH is calculated by taking 1 minus the ratio of how long
#' any of the child clusters of the cluster \emph{x} belongs to "survives"
#' changes in \emph{lambda} to the highest \emph{lambda} threshold of x, above
#' which x becomes a noise point.
#'
#' Scores close to 1 indicate outliers. For more details on the motivation for
#' this calculation, see Campello et al (2015).
#'
#' @aliases glosh GLOSH
#' @family Outlier Detection Functions
#'
#' @param x an [hclust] object, data matrix, or [dist] object.
#' @param k size of the neighborhood.
#' @param ... further arguments are passed on to [kNN()].
#' @return A numeric vector of length equal to the size of the original data
#' set containing GLOSH values for all data points.
#' @author Matt Piekenbrock
#'
#' @references Campello, Ricardo JGB, Davoud Moulavi, Arthur Zimek, and Joerg
#' Sander. Hierarchical density estimates for data clustering, visualization,
#' and outlier detection. _ACM Transactions on Knowledge Discovery from Data
#' (TKDD)_ 10, no. 1 (2015).
#' \doi{10.1145/2733381}
#' @keywords model
#' @examples
#' set.seed(665544)
#' n <- 100
#' x <- cbind(
#'   x=runif(10, 0, 5) + rnorm(n, sd = 0.4),
#'   y=runif(10, 0, 5) + rnorm(n, sd = 0.4)
#'   )
#'
#' ### calculate GLOSH score
#' glosh <- glosh(x, k = 3)
#'
#' ### distribution of outlier scores
#' summary(glosh)
#' hist(glosh, breaks = 10)
#'
#' ### simple function to plot point size is proportional to GLOSH score
#' plot_glosh <- function(x, glosh){
#'   plot(x, pch = ".", main = "GLOSH (k = 3)")
#'   points(x, cex = glosh*3, pch = 1, col = "red")
#'   text(x[glosh > 0.80, ], labels = round(glosh, 3)[glosh > 0.80], pos = 3)
#' }
#' plot_glosh(x, glosh)
#'
#' ### GLOSH with any hierarchy
#' x_dist <- dist(x)
#' x_sl <- hclust(x_dist, method = "single")
#' x_upgma <- hclust(x_dist, method = "average")
#' x_ward <- hclust(x_dist, method = "ward.D2")
#'
#' ## Compare what different linkage criterion consider as outliers
#' glosh_sl <- glosh(x_sl, k = 3)
#' plot_glosh(x, glosh_sl)
#'
#' glosh_upgma <- glosh(x_upgma, k = 3)
#' plot_glosh(x, glosh_upgma)
#'
#' glosh_ward <- glosh(x_ward, k = 3)
#' plot_glosh(x, glosh_ward)
#'
#' ## GLOSH is automatically computed with HDBSCAN
#' all(hdbscan(x, minPts = 3)$outlier_scores == glosh(x, k = 3))
#' @export
glosh <- function(x, k = 4, ...) {
  if (inherits(x, "data.frame"))
    x <- as.matrix(x)

  # get n
  if (inherits(x, "dist") || inherits(x, "matrix")) {
    if (inherits(x, "dist"))
      n <- attr(x, "Size")
    else
      n <- nrow(x)
    # get k nearest neighbors + distances
    d <- kNN(x, k - 1, ...)
    x_dist <-
      if (inherits(x, "dist"))
        x
    else
      dist(x, method = "euclidean") # copy since mrd changes by reference!
    mrd <- mrd(x_dist, d$dist[, k - 1])

    # need to assemble hclust object manually
    mst <- prims(mrd, n)
    hc <- hclustMergeOrder(mst, order(mst[, 3]))
  } else if (inherits(x, "hclust")) {
    hc <- x
    n <- nrow(hc$merge) + 1
  }
  else
    stop("x needs to be a matrix, dist, or hclust object!")

  if (k < 2 || k >= n)
    stop("k has to be larger than 1 and smaller than the number of points")

  res <- computeStability(hc, k, compute_glosh = TRUE)

  # return
  attr(res, "glosh")
}
