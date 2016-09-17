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

  ### extract clustering (keep hierarchical OPTICSXi structure)
  if(is(cl, "optics") && !is.null(cl$clusters_xi)) clusters_xi <- cl
  else clusters_xi <- NULL

  if(is.list(cl)) cl <- cl$cluster
  if(!is.numeric(cl)) stop("could not get cluster assignment vector from cl.")

  #if(is.null(col)) col <- c("#000000FF", rainbow(n=max(cl)))
  if(is.null(col)) col <- palette()
  if(max(cl)+1L > length(col)) warning("not enought colors. some colors will be reused.")
  
  # Plot hierarchical clustering
  if (inherits(opt,what="optics") && opt$min == FALSE) { # Short circuiting should allow element access
    plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)
    if (!is.null(opt$clusters_xi) && nrow(opt$clusters_xi) > 0) {
      for(i in 1:nrow(opt$clusters_xi)) {
        rng <- opt$clusters_xi[i,1:2]
        d <- x[opt$order[rng$start : rng$end],]
        ch <- chull(d)
        ch <- c(ch, ch[1])
        lines(d[ch,], col = col[i %% length(col) +1L], lwd=lwd_hull)
      }
    }
  } else { # Simple 'flat' clustering 
    # Double-check: If minimum was manually set to true at some point but was OPTICS-XI was originally 
    # run with minimum set to false (default), need to recompute the clustering array
    if (inherits(cl,what="optics") && cl$min == TRUE) { cl <- extractXiClusters(opt, TRUE)$cluster }
  plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)

  col_poly <- adjustcolor(col, alpha.f = alpha)
  border <- col

  ## no border?
  if(is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
    hull_lwd <- 1
    border <- NA
  }


  for(i in 1:max(cl)) {

    ### use all the points for OPTICSXi's hierarchical structure
    if(is.null(clusters_xi)) d <- x[cl==i,]
    else d <- x[clusters_xi$order[clusters_xi$clusters_xi$start[i] : clusters_xi$clusters_xi$end[i]],]

    ch <- chull(d)
    ch <- c(ch, ch[1])
    if(!solid)
      lines(d[ch,], col = col[i %% length(col) +1L], lwd=hull_lwd, lty=hull_lty)
    else
      polygon(d[ch,], col = col_poly[i %% length(col_poly) +1L],
        lwd=hull_lwd, lty=hull_lty, border = border[i %% length(col_poly) +1L])
  }
  }
}
