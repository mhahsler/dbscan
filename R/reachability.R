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
as.reachability <- function(object, ...) UseMethod("as.reachability")

# Simple print method for reachability objects
print.reachability <- function(x, ...) {
  avg_reach <- mean(x$reachdist[which(x$reachdist != Inf)], na.rm = T)
  cat("Reachability collection for ", length(x$order), " objects.\n", 
      "Avg minimum reachability distance: ", avg_reach, "\n", 
      "Available Fields: order, reachdist", sep="")
}

# Converting from dendrogram --> reachability plot 
as.reachability.dendrogram <- function(object, ...) {
  if (!inherits(object, "dendrogram")) stop("The as.reachability method requires a dendrogram object.")
  # Rcpp doesn't seem to import attributes well for vectors 
  fix_x <- dendrapply(object, function(leaf) { 
    new_leaf <- as.list(leaf); attributes(new_leaf) <- attributes(leaf); new_leaf 
  })
  res <- dendrogram_to_reach(fix_x)
  # Refix the ordering 
  res$reachdist <- res$reachdist[order(res$order)]
  return(res)
}

# Dendrogram --> Reachability conversion 
as.dendrogram.reachability <- function(object, ...) {
  if(length(which(object$reachdist == Inf)) > 1) stop("Multiple Infinite reachability distances found. Reachability plots can only be converted if they contain 
                                                     enough information to fully represent the dendrogram structure. If using OPTICS, a larger eps value 
                                                     (such as Inf) may be needed in the parameterization.")
  #dup_x <- object
  c_order <- order(object$reachdist) - 1
  # dup_x$order <- dup_x$order - 1
  #q_order <- sapply(c_order, function(i) which(dup_x$order == i))
  res <- reach_to_dendrogram(object, c_order)
  # res <- dendrapply(res, function(leaf) { new_leaf <- leaf[[1]]; attributes(new_leaf) <- attributes(leaf); new_leaf })
  res
}

# Plotting method for reachability objects.
plot.reachability <- function(x, order_labels = FALSE, xlab="Order", 
                              ylab = "Reachability dist.", main = "Reachability Plot", ...) {
  if (is.null(x$order) || is.null(x$reachdist)) stop("reachability objects need 'reachdist' and 'order' fields")
  plot(x$reachdist[x$order], xlab = xlab, ylab = ylab, main = main, type="h", ...)
  lines(x = c(1, 1), y = c(0, max(x$reachdist[x$reachdist != Inf])), lty=3)
  if (order_labels) { text(x = 1:length(x$order), y = x$reachdist[x$order], labels = x$order, pos=3) }
}
