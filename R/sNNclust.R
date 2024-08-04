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


#' Shared Nearest Neighbor Clustering
#'
#' Implements the shared nearest neighbor clustering algorithm by Ertoz,
#' Steinbach and Kumar (2003).
#'
#' **Algorithm:**
#'
#' 1. Constructs a shared nearest neighbor graph for a given k. The edge
#' weights are the number of shared k nearest neighbors (in the range of
#' \eqn{[0, k]}).
#'
#' 2. Find each points SNN density, i.e., the number of points which have a
#' similarity of `eps` or greater.
#'
#' 3. Find the core points, i.e., all points that have an SNN density greater
#' than `MinPts`.
#'
#' 4. Form clusters from the core points and assign border points (i.e.,
#' non-core points which share at least `eps` neighbors with a core point).
#'
#' Note that steps 2-4 are equivalent to the DBSCAN algorithm (see [dbscan()])
#' and that `eps` has a different meaning than for DBSCAN. Here it is
#' a threshold on the number of shared neighbors (see [sNN()])
#' which defines a similarity.
#'
#' @aliases sNNclust snnclust
#' @family clustering functions
#'
#' @param x a data matrix/data.frame (Euclidean distance is used), a
#' precomputed [dist] object or a kNN object created with [kNN()].
#' @param k Neighborhood size for nearest neighbor sparsification to create the
#' shared NN graph.
#' @param eps Two objects are only reachable from each other if they share at
#' least `eps` nearest neighbors. Note: this is different from the `eps` in DBSCAN!
#' @param minPts minimum number of points that share at least `eps`
#' nearest neighbors for a point to be considered a core points.
#' @param borderPoints should border points be assigned to clusters like in
#' [DBSCAN]?
#' @param ...  additional arguments are passed on to the k nearest neighbor
#' search algorithm. See [kNN()] for details on how to control the
#' search strategy.
#'
#' @return A object of class `general_clustering` with the following
#' components:
#' \item{cluster }{A integer vector with cluster assignments. Zero
#' indicates noise points.}
#' \item{type }{ name of used clustering algorithm.}
#' \item{param }{ list of used clustering parameters. }
#'
#' @author Michael Hahsler
#'
#' @references Levent Ertoz, Michael Steinbach, Vipin Kumar, Finding Clusters
#' of Different Sizes, Shapes, and Densities in Noisy, High Dimensional Data,
#' _SIAM International Conference on Data Mining,_ 2003, 47-59.
#' \doi{10.1137/1.9781611972733.5}
#' @keywords model clustering
#' @examples
#' data("DS3")
#'
#' # Out of k = 20 NN 7 (eps) have to be shared to create a link in the sNN graph.
#' # A point needs a least 16 (minPts) links in the sNN graph to be a core point.
#' # Noise points have cluster id 0 and are shown in black.
#' cl <- sNNclust(DS3, k = 20, eps = 7, minPts = 16)
#' cl
#'
#' clplot(DS3, cl)
#'
#' @export
sNNclust <- function(x, k, eps, minPts, borderPoints = TRUE, ...) {
  nn <- sNN(x, k = k, jp = TRUE, ...)

  # convert into a frNN object which already enforces eps
  nn_list <- lapply(seq_len(nrow(nn$id)),
    FUN = function(i) unname(nn$id[i, nn$shared[i, ] >= eps]))
  snn <- structure(list(id = nn_list, eps = eps, metric = nn$metric),
    class = c("NN", "frNN"))

  # run dbscan
  cl <- dbscan(snn, minPts = minPts, borderPoints = borderPoints)

  structure(list(cluster = cl$cluster,
    type = "SharedNN clustering",
    param = list(k = k, eps = eps, minPts = minPts, borderPoints = borderPoints),
    metric = cl$metric),
    class = "general_clustering")
}
