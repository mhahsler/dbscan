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

## Framework for Optimal Selection of Clusters
extractFOSC <- function(x, min_clsize = 2L, constraints = NA, alpha = 0){
  if (!class(x) %in% c("hclust")) stop("extractFOSC expects 'x' to be a valid hclust object.")
  if (!missing(constraints) && !class(constraints) %in% c("list", "integer", "numeric", "matrix")) { 
    stop("extractFOSC expects constraints to be either an adjacency list or adjacency matrix of constraints.")
  }
  if (!min_clsize >= 2L) stop("min_clsize must be at least 2.") 
  if (alpha < 0 || alpha > 1) stop("alpha can only takes values between [0, 1].")
  n <- nrow(x$merge) + 1L
  
  ## Make sure distances are non-zero to prevent NaN from divide by zero
  if (any(x$height == 0.0)) { x$height[x$height == 0.0] <- .Machine$double.eps }
    
  ## First step for both unsupervised and semisupervised - compute stability scores
  cl_tree <- computeStability(x, min_clsize)
  
  ## Unsupervised Extraction
  if (missing(constraints)){
    cl_tree <- extractUnsupervised(cl_tree)
  } 
  ## Semi-supervised Extraction
  else {
    ## If given as adjacency-list form
    if (is.list(constraints)) { cl_tree <- extractSemiSupervised(cl_tree, constraints, alpha) } 
    ## Adjacency matrix given (probably from dist object), retrieve adjacency list form
    else if (is.vector(constraints)){
      constraints_list <-  distToAdjacency(constraints, n)
      cl_tree <- extractSemiSupervised(cl_tree, constraints_list, alpha)
    } 
    ## Full nxn adjacency-matrix given, give warning and retrieve adjacency list form
    else if (is.matrix(constraints)){
      constraints <- constraints[lower.tri(constraints)]
      constraints_list <- distToAdjacency(constraints, n)
      cl_tree <- extractSemiSupervised(cl_tree, constraints_list, alpha)
    } else {
      stop(paste("'extractFOSC' doesn't know how to handle constraints of type", class(constraints)))
    }
  }
  ## Maximum possible tree stability; acts as a normalizing constant, if needed
  max_stability <- if (is.null(attr(cl_tree, "max_stability"))) 1 else attr(cl_tree, "max_stability")
  
  ## CL Tracker records the cluster at each merge step, going up 
  cl_track <- attr(cl_tree, "cl_tracker")
  
  ## Each merge step produces:
  ## 1. A stability score (unsupervised)
  ## 2. A constraint score (semi-supervised)
  ## 3. A total score, which may be identical to one of the above, or a convex combination thereof (if alpha is specified)
  stability_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$stability))
  constraint_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$vscore))
  total_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$score))
  
  ## Add scores to the existing hclust hierarchy 
  out <- append(x, list( #"cluster"=cl_track,
                        "stability"=stability_score,
                        "constraint"=constraint_score,
                        "total"=total_score))
  class(out) <- "hclust"
  
  ## Return the new merge information, and the resulting flat clustering 
  res <- list(cluster=attr(cl_tree, "cluster"), hc=out)
  # attr(res, ".internal") <- cl_tree # TODO: make better diagnostic and visualization tools
  return(res)
}