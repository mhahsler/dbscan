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

#' Plot Clusters
#'
#' This function produces a two-dimensional scatter plot of data points
#' and colors the data points according to a supplied clustering. Noise points
#' are marked as `x`. `hullplot()` also adds convex hulls to clusters.
#'
#' @name hullplot
#' @aliases hullplot clplot
#'
#' @param x a data matrix. If more than 2 columns are provided, then the data
#' is plotted using the first two principal components.
#' @param cl a clustering. Either a numeric cluster assignment vector or a
#' clustering object (a list with an element named `cluster`).
#' @param col colors used for clusters. Defaults to the standard palette.  The
#' first color (default is black) is used for noise/unassigned points (cluster
#' id 0).
#' @param pch a vector of plotting characters. By default `o` is used for
#'   points and `x` for noise points.
#' @param cex expansion factor for symbols.
#' @param hull_lwd,hull_lty line width and line type used for the convex hull.
#' @param main main title.
#' @param solid,alpha draw filled polygons instead of just lines for the convex
#' hulls? alpha controls the level of alpha shading.
#' @param ...  additional arguments passed on to plot.
#' @author Michael Hahsler
#' @keywords plot clustering
#' @examples
#' set.seed(2)
#' n <- 400
#'
#' x <- cbind(
#'   x = runif(4, 0, 1) + rnorm(n, sd = 0.1),
#'   y = runif(4, 0, 1) + rnorm(n, sd = 0.1)
#'   )
#' cl <- rep(1:4, time = 100)
#'
#'
#' ### original data with true clustering
#' clplot(x, cl, main = "True clusters")
#' hullplot(x, cl, main = "True clusters")
#' ### use different symbols
#' hullplot(x, cl, main = "True clusters", pch = cl)
#' ### just the hulls
#' hullplot(x, cl, main = "True clusters", pch = NA)
#' ### a version suitable for b/w printing)
#' hullplot(x, cl, main = "True clusters", solid = FALSE,
#'   col = c("grey", "black"), pch = cl)
#'
#'
#' ### run some clustering algorithms and plot the results
#' db <- dbscan(x, eps = .07, minPts = 10)
#' clplot(x, db, main = "DBSCAN")
#' hullplot(x, db, main = "DBSCAN")
#'
#' op <- optics(x, eps = 10, minPts = 10)
#' opDBSCAN <- extractDBSCAN(op, eps_cl = .07)
#' hullplot(x, opDBSCAN, main = "OPTICS")
#'
#' opXi <- extractXi(op, xi = 0.05)
#' hullplot(x, opXi, main = "OPTICSXi")
#'
#' # Extract minimal 'flat' clusters only
#' opXi <- extractXi(op, xi = 0.05, minimum = TRUE)
#' hullplot(x, opXi, main = "OPTICSXi")
#'
#' km <- kmeans(x, centers = 4)
#' hullplot(x, km, main = "k-means")
#'
#' hc <- cutree(hclust(dist(x)), k = 4)
#' hullplot(x, hc, main = "Hierarchical Clustering")
#' @export
hullplot <- function(x,
  cl,
  col = NULL,
  pch = NULL,
  cex = 0.5,
  hull_lwd = 1,
  hull_lty = 1,
  solid = TRUE,
  alpha = .2,
  main = "Convex Cluster Hulls",
  ...) {
  ### handle d>2 by using PCA
  if (ncol(x) > 2)
    x <- prcomp(x)$x

  ### extract clustering (keep hierarchical xICSXi structure)
  if (inherits(cl, "xics") || "clusters_xi" %in% names(cl)) {
    clusters_xi <- cl$clusters_xi
    cl_order <- cl$order
  } else
    clusters_xi <- NULL

  if (is.list(cl))
    cl <- cl$cluster
  if (!is.numeric(cl))
    stop("Could not get cluster assignment vector from cl.")

  #if(is.null(col)) col <- c("#000000FF", rainbow(n=max(cl)))
  if (is.null(col))
    col <- palette()

  # Note: We use the first color for noise points
  if (length(col) == 1L)
    col <- c(col, col)
  col_noise <- col[1]
  col <- col[-1]


  if (max(cl) > length(col)) {
    warning("Not enough colors. Some colors will be reused.")
    col <- rep(col, length.out = max(cl))
  }

  # mark noise points
  pch <- pch %||% ifelse(cl == 0L, 4L, 1L)

  plot(x[, 1:2],
    col = c(col_noise, col)[cl + 1L],
    pch = pch,
    cex = cex,
    main = main,
    ...)

  col_poly <- adjustcolor(col, alpha.f = alpha)
  border <- col

  ## no border?
  if (is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
    hull_lwd <- 1
    border <- NA
  }

  if (inherits(cl, "xics") || "clusters_xi" %in% names(cl)) {
    ## This is necessary for larger datasets: Ensure largest is plotted first
    clusters_xi <-
      clusters_xi[order(-(clusters_xi$end - clusters_xi$start)), ] # Order by size (descending)
    ci_order <- clusters_xi$cluster_id
  } else {
    ci_order <- 1:max(cl)
  }

  for (i in seq_along(ci_order)) {
    ### use all the points for xICSXi's hierarchical structure
    if (is.null(clusters_xi)) {
      d <- x[cl == i, , drop = FALSE]
    } else {
      d <-
        x[cl_order[clusters_xi$start[i]:clusters_xi$end[i]], , drop = FALSE]
    }

    ch <- chull(d)
    ch <- c(ch, ch[1])
    if (!solid) {
      lines(d[ch, ],
            col = border[ci_order[i]],
            lwd = hull_lwd,
            lty = hull_lty)
    } else {
      polygon(
        d[ch, ],
        col = col_poly[ci_order[i]],
        lwd = hull_lwd,
        lty = hull_lty,
        border = border[ci_order[i]]
      )
    }
  }
}

#' @rdname hullplot
#' @export
clplot <- function(x,
                   cl,
                   col = NULL,
                   pch = NULL,
                   cex = 0.5,
                   main = "Cluster Plot",
                   ...)
  hullplot(x, cl = cl, col = col, pch = pch, cex = cex, main = main,
          solid = FALSE, hull_lwd = NA)
