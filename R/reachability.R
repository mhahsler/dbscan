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


#' Density Reachability Structures for OPTICS
#'
#' Class `reachability` provides general functions for representing various
#' hierarchical representations as reachability plots, as originally defined
#' by Ankerst et al (1999) for [OPTICS]. Methods include fast implementations of the
#' conversion algorithms introduced by Sanders et al (2003) to convert between
#' a [dendrogram] and a reachability object.
#'
#' Dendrograms are a popular visualization tool for representing hierarchical
#' relationships. In agglomerative clustering, dendrograms can be constructed
#' using a variety of linkage criterion (such as single or complete linkage),
#' many of which are frequently used to
#'
#' 1.  visualize the density-based
#' relationships in the data or
#' 2.  extract cluster labels from the data the
#' dendrogram represents.
#'
#' The original ordering algorithm OPTICS as described by Ankerst et al (1999)
#' introduced the notion of 2-dimensional representation of so-called
#' _density-reachability_ that was shown to be useful for data visualization.
#' This representation was shown to essentially convey the same information as
#' the more traditional dendrogram structure by Sanders et al (2003).
#'
#' Different hierarchical representations, such as dendrograms or reachability
#' plots, may be preferable depending on the context. In smaller datasets,
#' cluster memberships may be more easily identifiable through a dendrogram
#' representation, particularly is the user is already familiar with tree-like
#' representations. For larger datasets however, a reachability plot may be
#' preferred for visualizing macro-level density relationships.
#'
#' The central idea behind a reachability plot is that the ordering in which
#' points are plotted identifies underlying hierarchical density
#' representation. OPTICS linearly orders the data points such that points
#' which are spatially closest become neighbors in the ordering. Valleys
#' represent clusters, which can be represented hierarchically. Although the
#' ordering is crucial to the structure of the reachability plot, its important
#' to note that OPTICS, like DBSCAN, is not entirely deterministic and, just
#' like the dendrogram, isomorphisms may exist
#'
#' A variety of cluster extraction methods have been proposed using
#' reachability plots. Because both cluster extraction depend directly on the
#' ordering OPTICS produces, they are part of the optics interface.
#' Nonetheless, reachability plots can be created directly from other types of
#' linkage trees, and vice versa.
#'
#' @name reachability
#' @aliases reachability print.reachability
#'
#' @param object any object that can be coerced to class
#' `reachability`, such as an object of class [optics] or [stats::dendrogram].
#' @param x object of class `reachability`.
#' @param order_labels whether to plot text labels for each points reachability
#' distance.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main Title of the plot.
#' @param ...  graphical parameters are passed on to `plot()`,
#'   or arguments for other methods.
#'
#' @return An object of class `reachability` with components:
#' \item{order }{order to use for the data points in `x`. }
#' \item{reachdist }{reachability distance for each data point in `x`. }
#'
#' @author Matthew Piekenbrock
#' @seealso [optics()], [stats::hclust()].
#' @references Ankerst, M., M. M. Breunig, H.-P. Kriegel, J. Sander (1999).
#' OPTICS: Ordering Points To Identify the Clustering Structure. _ACM
#' SIGMOD international conference on Management of data._ ACM Press. pp.
#' 49--60.
#'
#' Sander, J., X. Qin, Z. Lu, N. Niu, and A. Kovarsky (2003). Automatic
#' extraction of clusters from hierarchical clustering representations.
#' _Pacific-Asia Conference on Knowledge Discovery and Data Mining._
#' Springer Berlin Heidelberg.
#' @keywords model clustering hierarchical clustering
#' @examples
#' set.seed(2)
#' n <- 20
#'
#' x <- cbind(
#'   x = runif(4, 0, 1) + rnorm(n, sd = 0.1),
#'   y = runif(4, 0, 1) + rnorm(n, sd = 0.1)
#' )
#'
#' plot(x, xlim = range(x), ylim = c(min(x) - sd(x), max(x) + sd(x)), pch = 20)
#' text(x = x, labels = 1:nrow(x), pos = 3)
#'
#' ### run OPTICS
#' res <- optics(x, eps = 10,  minPts = 2)
#' res
#'
#' ### plot produces a reachability plot
#' plot(res)
#'
#' ### Extract reachability components from OPTICS
#' reach <- as.reachability(res)
#' reach
#'
#' ### plot still produces a reachability plot; points ids
#' ### (rows in the original data) can be displayed with order_labels = TRUE
#' plot(reach, order_labels = TRUE)
#'
#' ### Reachability objects can be directly converted to dendrograms
#' dend <- as.dendrogram(reach)
#' dend
#' plot(dend)
#'
#' ### A dendrogram can be converted back into a reachability object
#' plot(as.reachability(dend))
NULL


print.reachability <- function(x, ...) {
  avg_reach <- mean(x$reachdist[which(x$reachdist != Inf)], na.rm = T)
  cat(
    "Reachability collection for ",
    length(x$order),
    " objects.\n",
    "Avg minimum reachability distance: ",
    avg_reach,
    "\n",
    "Available Fields: order, reachdist",
    sep = ""
  )
}

#' @rdname reachability
as.reachability.dendrogram <- function(object, ...) {
  if (!inherits(object, "dendrogram"))
    stop("The as.reachability method requires a dendrogram object.")
  # Rcpp doesn't seem to import attributes well for vectors
  fix_x <- dendrapply(object, function(leaf) {
    new_leaf <-
      as.list(leaf)
    attributes(new_leaf) <- attributes(leaf)
    new_leaf
  })
  res <- dendrogram_to_reach(fix_x)
  # Refix the ordering
  res$reachdist <- res$reachdist[order(res$order)]

  return(res)
}

#' @rdname reachability
plot.reachability <- function(x,
  order_labels = FALSE,
  xlab = "Order",
  ylab = "Reachability dist.",
  main = "Reachability Plot",
  ...) {
  if (is.null(x$order) ||
      is.null(x$reachdist))
    stop("reachability objects need 'reachdist' and 'order' fields")
  plot(
    x$reachdist[x$order],
    xlab = xlab,
    ylab = ylab,
    main = main,
    type = "h",
    ...
  )
  lines(x = c(1, 1),
    y = c(0, max(x$reachdist[x$reachdist != Inf])),
    lty = 3)
  if (order_labels) {
    text(
      x = 1:length(x$order),
      y = x$reachdist[x$order],
      labels = x$order,
      pos = 3
    )
  }
}

# calculate midpoints for dendrogram
# from stats, but not exported
# see stats:::midcache.dendrogram

.midcache.dendrogram <- function (x, type = "hclust", quiet = FALSE)
{
  type <- match.arg(type)
  stopifnot(inherits(x, "dendrogram"))
  verbose <- getOption("verbose", 0) >= 2
  setmid <- function(d, type) {
    depth <- 0L
    kk <- integer()
    jj <- integer()
    dd <- list()
    repeat {
      if (!is.leaf(d)) {
        k <- length(d)
        if (k < 1)
          stop("dendrogram node with non-positive #{branches}")
        depth <- depth + 1L
        if (verbose)
          cat(sprintf(" depth(+)=%4d, k=%d\n", depth,
            k))
        kk[depth] <- k
        if (storage.mode(jj) != storage.mode(kk))
          storage.mode(jj) <- storage.mode(kk)
        dd[[depth]] <- d
        d <- d[[jj[depth] <- 1L]]
        next
      }
      while (depth) {
        k <- kk[depth]
        j <- jj[depth]
        r <- dd[[depth]]
        r[[j]] <- unclass(d)
        if (j < k)
          break
        depth <- depth - 1L
        if (verbose)
          cat(sprintf(" depth(-)=%4d, k=%d\n", depth,
            k))
        midS <- sum(vapply(r, .midDend, 0))
        if (!quiet && type == "hclust" && k != 2)
          warning("midcache() of non-binary dendrograms only partly implemented")
        attr(r, "midpoint") <- (.memberDend(r[[1L]]) +
            midS) / 2
        d <- r
      }
      if (!depth)
        break
      dd[[depth]] <- r
      d <- r[[jj[depth] <- j + 1L]]
    }
    d
  }
  setmid(x, type = type)
}

.midDend <- function (x)
{
  if (is.null(mp <- attr(x, "midpoint")))
    0
  else
    mp
}

.memberDend <- function (x)
{
  r <- attr(x, "x.member")
  if (is.null(r)) {
    r <- attr(x, "members")
    if (is.null(r))
      r <- 1L
  }
  r
}
