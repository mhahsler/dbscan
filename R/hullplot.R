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

hullplot <- function(x, cl, col = NULL,
  cex = 0.5, hull_lwd = 1, hull_lty = 1, solid = TRUE, alpha = .2,
  main = "Convex Cluster Hulls", ...) {

  ### handle d>2 by using PCA
  if(ncol(x)>2) x <- prcomp(x)$x

  ### extract clustering (keep hierarchical xICSXi structure)
  if(is(cl, "xics") || "clusters_xi" %in% names(cl)) {
    clusters_xi <- cl$clusters_xi
    cl_order <- cl$order
  } else clusters_xi <- NULL

  if(is.list(cl)) cl <- cl$cluster
  if(!is.numeric(cl)) stop("Could not get cluster assignment vector from cl.")

  #if(is.null(col)) col <- c("#000000FF", rainbow(n=max(cl)))
  if(is.null(col)) col <- palette()
  if(max(cl)+1L > length(col)) warning("Not enough colors. Some colors will be reused.")

  plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)
  col_poly <- adjustcolor(col, alpha.f = alpha)
  border <- col

  ## no border?
  if(is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
    hull_lwd <- 1
    border <- NA
  }

  if(is(cl, "xics") || "clusters_xi" %in% names(cl)) {
    ## This is necessary for larger datasets: Ensure largest is plotted first
    clusters_xi <- clusters_xi[order(-(clusters_xi$end-clusters_xi$start)),] # Order by size (descending)
    ci_order <- clusters_xi$cluster_id
  } else { ci_order <- 1:max(cl) }

  for(i in 1:length(ci_order)) {

    ### use all the points for xICSXi's hierarchical structure
    if(is.null(clusters_xi)) { d <- x[cl==i,]
    } else { d <- x[cl_order[clusters_xi$start[i] : clusters_xi$end[i]],] }

    ch <- chull(d)
    ch <- c(ch, ch[1])
    if(!solid) { lines(d[ch,], col = col[ci_order[i] %% length(col) +1L], lwd=hull_lwd, lty=hull_lty)
    } else { polygon(d[ch,], col = col_poly[ci_order[i] %% length(col_poly) +1L],
                     lwd=hull_lwd, lty=hull_lty, border = border[ci_order[i] %% length(col_poly) +1L]) }
  }
}
