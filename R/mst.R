#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#' Minimum Spanning Tree Hierarchies
#'
#' Computes the minimum spanning tree (MST) of a graph with Prim's algorithm
#' and returns it as a single-linkage hierarchy ([hclust] object).
#'
#' * `mst_dist()` computes the MST of the complete graph defined by a [dist] object.
#'
#' * `mst_mrdist()` computes the MST of the graph of mutual reachability
#'   distances (see [mrdist()]) directly from a data matrix and `minPts`. The
#'   core distances are calculated first (see [coredist()]). The mutual
#'   reachability distances are computed on the fly, so **no complete distance
#'   matrix is created**. This saves memory (\eqn{O(n d)} instead of
#'   \eqn{O(n^2)}) and is used by [hdbscan()] for data matrices. The result is
#'   the same as `mst_dist(mrdist(x, minPts))`. Only Euclidean distances are
#'   supported.
#'
#' MSTs are not unique if there are ties in the distances. Many ties are very
#' likely for mutual reachability distances constructed as
#' `mrd(a, b) = max(coredist(a), coredist(b), dist(a, b))`.
#' Ties in this implementation are broken in order of the rows in the
#' input data which is consistent with many otehr implementations.
#'
#' @family HDBSCAN functions
#'
#' @param x for `mst_dist()` a [dist] object with finite distances for at least two
#'   objects. For `mst_mrdist()` a data matrix (or an object that can be
#'   converted to a numeric matrix) with finite values.
#' @param minPts integer; minimum number of points used to calculate the core
#'   distances (see [coredist()]).
#' @param coredist optional vector with the core distance of each point. If
#'   `NULL`, then it is calculated with [coredist()].
#'
#' @returns an [hclust] object representing the single-linkage hierarchy
#'  of the MST. Note: the rows of the `merge` matrix keep
#'  the order in which the MST was built and are **not** normalized to the
#'  convention documented in [hclust] (singletons before clusters, singletons in
#'  increasing index order, clusters in increasing step order).
#'  The raw order is kept so that cluster ids stay consistent with [hdbscan()].
#'  Code that relies on this convention (e.g., other packages reading `$merge`
#'  directly) may fail or give unexpected results. The `method` component is
#'  `"robust single"`. `dist.method` is taken from the `"method"` attribute of
#'  `x` for `mst_dist()` (it is omitted if `x` has none) and is
#'  `"mutual reachability (euclidean)"` for `mst_mrdist()`.
#'  `mst_mrdist()` also stores `minPts` in the result.
#'
#' @seealso [hdbscan()], [mrdist()], [extractFOSC()]
#' @examples
#' data("moons")
#'
#' ## MST of a distance matrix
#' hc <- mst_dist(dist(moons))
#' hc
#'
#' ## MST of the mutual reachability graph without creating the distance matrix.
#' hc_mr <- mst_mrdist(moons, minPts = 5)
#'
#' ## same as
#' hc_mr2 <- mst_dist(mrdist(moons, minPts = 5))
#' all.equal(hc_mr$merge, hc_mr2$merge)
#'
#' ## extract clusters (this is HDBSCAN)
#' cl <- extractFOSC(hc_mr)
#' hullplot(moons, cl$cluster, main = "HDBSCAN")
#' @export
mst_dist <- function(x) {
  .check_dist(x)
  n <- attr(x, "Size")
  if (n < 2L)
    stop("x needs to contain at least 2 objects.")
  if (!all(is.finite(x)))
    stop("x needs to contain only finite distances.")

  m <- mst_dist_cpp(x, n)
  hc <- hclustMergeOrder(m, order(m[, 3]))

  hc$labels <- attr(x, "Labels")
  hc$call <- match.call()
  hc$dist.method <- attr(x, "method")
  return(hc)
}

#' @rdname mst_dist
#' @export
mst_mrdist <- function(x, minPts, coredist = NULL) {
  minPts <- .validate_integer_scalar(minPts, "minPts", min = 2)
  x <- .as_finite_numeric_matrix(x)
  if (nrow(x) < 2L)
    stop("x needs to contain at least 2 objects.")

  if (is.null(coredist))
    coredist <- coredist(x, minPts)
  else if (length(coredist) != nrow(x))
    stop("coredist needs to have one value for each row of x.")

  m <- mst_mrd(x, coredist)
  hc <- hclustMergeOrder(m, order(m[, 3]))

  hc$labels <- rownames(x)
  hc$call <- match.call()
  hc$dist.method <- "mutual reachability (euclidean)"
  hc$minPts <- minPts
  return(hc)
}
