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
  cex = 0.5, lwd_hull = 2, main = "Convex Cluster Hulls", ...) {

  if(is.list(cl)) cl <- cl$cluster
  if(!is.numeric(cl)) stop("could not get cluster assignment vector from cl.")

  if(is.null(col)) col <- palette()
  if(max(cl)+1L > length(col)) warning("not enought colors. some colors will be reused.")

  plot(x[,1:2], col = col[cl %% length(col) +1L], cex = cex, main = main, ...)

  for(i in 1:max(cl)) {
    d <- x[cl==i,]
    ch <- chull(d)
    ch <- c(ch, ch[1])
    lines(d[ch,], col = col[i %% length(col) +1L], lwd=lwd_hull)
  }
}
