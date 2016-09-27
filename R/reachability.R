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


# Add new Generic for other class extensions support reachability plotting 
as.reachability <- function(x, ...) UseMethod("as.reachability")

# Simple print method for reachability objects
print.reachability <- function(x) {
  min_reach <- na.omit(x$reachdist)[which(na.omit(x$reachdist) != Inf)]
  cat("Reachability collection for ", length(x$order), " objects.\n", 
      "Avg minimum reachability distance: ", mean(min_reach), "\n", 
      "Field Specializations: order, reachdist, cluster", sep="")
}

# Dendrogram --> Reachability conversion 
as.dendrogram.reachability <- function(x) {
  if(length(which(x$reachdist == Inf)) > 1) stop(cat("Multiple Infinite reachability distances found. Reachability plots can only be converted if they contain 
                                                     enough information to fully represent the dendrogram structure. If using OPTICS, a larger eps value 
                                                     (such as Inf) may be needed in the parameterization."))
  dup_x <- x
  c_order <- order(dup_x$reachdist) - 1
  dup_x$order <- dup_x$order - 1
  reach_to_dendrogram(dup_x, sort(dup_x$reachdist), c_order)
}

# Converting from dendrogram --> reachability plot 
as.reachability.dendrogram <- function(x) {
  if (!inherits(x, "dendrogram")) stop("The as.reachability method requires a dendrogram object.")
  fix_x <- dendrapply(x, function(leaf) { 
    new_leaf <- as.list(leaf); attributes(new_leaf) <- attributes(leaf); new_leaf 
  })
  res <- dendrogram_to_reach(fix_x)
  return(res); 
}

# Plotting method for reachability objects.
plot.reachability <- function(x, order_labels = FALSE, ...) {
  if (is.null(x$cluster) || is.na(x$cluster)) {
    plot(x$reachdist[x$order], type="h", xlab="Order", 
         ylab = "Reachability dist.", main = "Reachability Plot", ...)
    lines(x = c(1, 1), y = c(0, max(x$reachdist[x$reachdist != Inf])), lty=2)
    if (order_labels) { text(x = 1:length(x$order), y = x$reachdist[x$order], labels = x$order, pos=3) }
  } else {
    if (is(x$cluster, "xics") || all(c("start", "end", "cluster_id") %in% names(x$cluster))) {
      y_max <- max(x$reachdist[which(x$reachdist != Inf)])
      # Sort clusters by size
      hclusters <- x$cluster[order(x$cluster$end - x$cluster$start),]
      def.par <- par(no.readonly = TRUE)
      
      # .1 means to leave 15% for the cluster lines
      par(mar= c(2, 4, 4, 2) + 0.1, omd = c(0, 1, .15, 1))
      y_increments <- (y_max/0.85*.15)/(nrow(hclusters)+1L)
      
      top_level <- extractClusterLabels(x$cluster, x$order)
      plot(x$reachdist[x$order], type="h", col=top_level[x$order]+1L,
           ylab = "Reachability dist.", xlab = NA, xaxt = "n",
           main = "Reachability Plot", yaxs="i", ylim=c(0,y_max), xaxt='n', ...)
      if (order_labels) { text(x = 1:length(x$order), y = x$reachdist[x$order], labels = x$order, pos=3) }
      
      i <- 1:nrow(hclusters)
      segments(x0=hclusters$start[i], y0=-(y_increments*i),
               x1=hclusters$end[i], col=hclusters$cluster_id[i]+1L, lwd=2, xpd=NA)
      ## Restore previous settings
      par(def.par)
    } else if (is.numeric(x$cluster)) { # 'flat' clustering given 
      plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
           ylab = "Reachability dist.", xlab = "Order",
           main = "Reachability Plot", ...)
      if (order_labels) { text(x = 1:length(x$order), y = x$reachdist[x$order], labels = x$order, pos=3) }
    }
  }
}
