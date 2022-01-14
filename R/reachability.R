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


#' Reachability Plots
#'
#' Reachability plots can be used to show hierarchical relationships between data points.
#' The idea was originally introduced by Ankerst et al (1999) for [OPTICS]. Later,
#' Sanders et al (2003) showed that the visualization is useful for other hierarchical
#' structures and introduced an algorithm to convert [dendrogram] representation to
#' reachability plots.
#'
#' A reachability plot displays the points as vertical bars, were the height is the
#' reachability distance between two consecutive points.
#' The central idea behind reachability plots is that the ordering in which
#' points are plotted identifies underlying hierarchical density
#' representation as mountains and valleys of high and low reachability distance.
#' The original ordering algorithm OPTICS as described by Ankerst et al (1999)
#' introduced the notion of reachability plots.
#'
#' OPTICS linearly orders the data points such that points
#' which are spatially closest become neighbors in the ordering. Valleys
#' represent clusters, which can be represented hierarchically. Although the
#' ordering is crucial to the structure of the reachability plot, its important
#' to note that OPTICS, like DBSCAN, is not entirely deterministic and, just
#' like the dendrogram, isomorphisms may exist
#'
#' Reachability plots were shown to essentially convey the same information as
#' the more traditional dendrogram structure by Sanders et al (2003). An dendrograms
#' can be converted into reachability plots.
#'
#' Different hierarchical representations, such as dendrograms or reachability
#' plots, may be preferable depending on the context. In smaller datasets,
#' cluster memberships may be more easily identifiable through a dendrogram
#' representation, particularly is the user is already familiar with tree-like
#' representations. For larger datasets however, a reachability plot may be
#' preferred for visualizing macro-level density relationships.
#'
#' A variety of cluster extraction methods have been proposed using
#' reachability plots. Because both cluster extraction depend directly on the
#' ordering OPTICS produces, they are part of the [optics()] interface.
#' Nonetheless, reachability plots can be created directly from other types of
#' linkage trees, and vice versa.
#'
#' _Note:_ The reachability distance for the first point is by definition not defined
#' (it has no preceeding point).
#' Also, the reachability distances can be undefined when a point does not have enough
#' neighbors in the epsilon neighborhood. We represent these undefined cases as `Inf`
#' and represent them in the plot as a dashed line.
#'
#' @name reachability_plot
#' @aliases reachability reachability_plot print.reachability
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
#' ### plot produces a reachability plot.
#' plot(res)
#'
#' ### Manually extract reachability components from OPTICS
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
    "Reachability plot collection for ",
    length(x$order),
    " objects.\n",
    "Avg minimum reachability distance: ",
    avg_reach,
    "\n",
    "Available Fields: order, reachdist",
    sep = ""
  )
}

#' @rdname reachability_plot
plot.reachability <- function(x,
  order_labels = FALSE,
  xlab = "Order",
  ylab = "Reachability dist.",
  main = "Reachability Plot",
  ...) {
  if (is.null(x$order) ||
      is.null(x$reachdist))
    stop("reachability objects need 'reachdist' and 'order' fields")
  reachdist <- x$reachdist[x$order]

  plot(
    reachdist,
    xlab = xlab,
    ylab = ylab,
    main = main,
    type = "h",
    ...
  )
  abline(v = which(is.infinite(reachdist)),
    lty = 3)
  if (order_labels) {
    text(
      x = 1:length(x$order),
      y = reachdist,
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

#' @rdname reachability_plot
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
