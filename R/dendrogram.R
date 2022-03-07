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

#' Coersions to Dendrogram
#'
#' Coerce [hclust], [hdbscan] or [reachability] objects to [dendrogram].
#'
#' The coersion from hclust is a faster reimplementation of [stats::as.dendrogram()].
#'
#' The coersion from [hdbscan] builds the non-simplified HDBSCAN hierarchy as a
#' dendrogram object.
#'
#' @name dendrogram
#' @aliases dendrogram
#'
#' @param object the object
#' @param ... further arguments
NULL

#' @rdname dendrogram
#' @export
as.dendrogram.hclust <- function(object, ...) {
  return(buildDendrogram(object))
}

#' @rdname dendrogram
#' @export
as.dendrogram.hdbscan <- function(object, ...) {
  return(buildDendrogram(object$hc))
}

#' @rdname dendrogram
#' @name reachability-coersion
#' @aliases as.reachability
NULL

#' @rdname dendrogram
#' @export
as.reachability <-
  function(object, ...)
    UseMethod("as.reachability")

#' @rdname dendrogram
#' @export
as.dendrogram.reachability <- function(object, ...) {
  if (length(which(object$reachdist == Inf)) > 1)
    stop(
      "Multiple Infinite reachability distances found. Reachability plots can only be converted if they contain enough information to fully represent the dendrogram structure. If using OPTICS, a larger eps value (such as Inf) may be needed in the parameterization."
    )
  #dup_x <- object
  c_order <- order(object$reachdist) - 1
  # dup_x$order <- dup_x$order - 1
  #q_order <- sapply(c_order, function(i) which(dup_x$order == i))
  res <- reach_to_dendrogram(object, c_order)
  # res <- dendrapply(res, function(leaf) { new_leaf <- leaf[[1]]; attributes(new_leaf) <- attributes(leaf); new_leaf })

  # add mid points for plotting
  res <- .midcache.dendrogram(res)

  res
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
  if (is.null(mp <- attr(x, "midpoint")))
    0
  else
    mp

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
