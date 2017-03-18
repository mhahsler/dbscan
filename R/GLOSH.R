#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler, Matthew Piekenbrock

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

glosh <- function(x, k = 4, ...) {
  
  # get n
  if (is(x, "dist") || is(x, "matrix")){
    if(is(x, "dist")) n <- attr(x, "Size")
    else n <- nrow(x)
    # get k nearest neighbors + distances
    d <- kNN(x, k - 1, ...)
    x_dist <-  if(is(x, "dist")) x else dist(x, method = "euclidean") # copy since mrd changes by reference!
    mrd <- mrd(x_dist, d$dist[, k - 1])
    
    # need to assemble hclust object manually
    mst <- prims(mrd, n)
    hc <- hclustMergeOrder(mst, order(mst[, 3]))
  } else if (is(x, "hclust")){
    hc <- x
    n <- nrow(hc$merge) + 1
  }
  else stop("x needs to be a matrix, dist, or hclust object!")
  
  if(k < 2 || k >= n)
    stop("k has to be larger than 1 and smaller than the number of points")
  
  res <- computeStability(hc, k, compute_glosh = TRUE)
  
  # return 
  attr(res, "glosh")
}
