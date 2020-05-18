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
extractFOSC <- function(x, constraints = NA, alpha = 0, minPts = 2L, prune_unstable = FALSE, validate_constraints = FALSE){
  if (!class(x) %in% c("hclust")) stop("extractFOSC expects 'x' to be a valid hclust object.")
  if (!missing(constraints) && !class(constraints) %in% c("list", "integer", "numeric", "matrix")) {
    stop("extractFOSC expects constraints to be either an adjacency list or adjacency matrix of constraints.")
  }
  if (!minPts >= 2) stop("minPts must be at least 2.")
  if (alpha < 0 || alpha > 1) stop("alpha can only takes values between [0, 1].")
  n <- nrow(x$merge) + 1L

  ## First step for both unsupervised and semisupervised - compute stability scores
  cl_tree <- computeStability(x, minPts)

  ## Unsupervised Extraction
  if (missing(constraints)){
    cl_tree <- extractUnsupervised(cl_tree, prune_unstable)
  }
  ## Semi-supervised Extraction
  else {
    ## If given as adjacency-list form
    if (is.list(constraints)) {
      ## Checks for proper indexing, symmetry of constraints, etc.
      if (validate_constraints) {
        is_valid <- max(as.integer(names(constraints))) < n
        is_valid <- is_valid && all(sapply(constraints, function(ilc) all(ilc <= n)))
        if (!is_valid){ stop("Detected constraint indices not in the interval [1, n]") }
        constraints <- validateConstraintList(constraints, n)
      }
      cl_tree <- extractSemiSupervised(cl_tree, constraints, alpha, prune_unstable)
    }
    ## Adjacency matrix given (probably from dist object), retrieve adjacency list form
    else if (is.vector(constraints)){
      if (!all(constraints %in% c(-1, 0, 1))){ stop("'extractFOSC' only accepts instance-level constraints. See ?extractFOSC for more details.") }
      ## Checks for proper integer labels, symmetry of constraints, length of vector, etc.
      if (validate_constraints) {
        is_valid <- length(constraints) == choose(n, 2)
        constraints_list <- validateConstraintList(distToAdjacency(constraints, n), n)
      } else {
        constraints_list <-  distToAdjacency(constraints, n)
      }
      cl_tree <- extractSemiSupervised(cl_tree, constraints_list, alpha, prune_unstable)
    }
    ## Full nxn adjacency-matrix given, give warning and retrieve adjacency list form
    else if (is.matrix(constraints)){
      if (!all(constraints %in% c(-1, 0, 1))){ stop("'extractFOSC' only accepts instance-level constraints. See ?extractFOSC for more details.") }
      if (!all(dim(constraints) == c(n, n))) { stop("Given matrix is not square.") }
      warning("Full nxn matrix given; extractFOCS does not support asymmetric relational constraints. Using lower triangular.")

      constraints <- constraints[lower.tri(constraints)]

      ## Checks for proper integer labels, symmetry of constraints, length of vector, etc.
      if (validate_constraints) {
        is_valid <- length(constraints) == choose(n, 2)
        constraints_list <- validateConstraintList(distToAdjacency(constraints, n), n)
      } else {
        constraints_list <- distToAdjacency(constraints, n)
      }
      cl_tree <- extractSemiSupervised(cl_tree, constraints_list, alpha, prune_unstable)
    } else {
      stop(paste("'extractFOSC' doesn't know how to handle constraints of type", class(constraints)))
    }
  }
  total_stab <- if (is.null(attr(cl_tree, "total_stability"))) 1 else attr(cl_tree, "total_stability")
  cl_track <- attr(cl_tree, "cl_tracker")
  stability_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$stability))
  constraint_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$vscore))
  total_score <- unlist(sapply(cl_track, function(cid) cl_tree[[as.character(cid)]]$vscore))
  out <- append(x, list("cluster"=cl_track,
                        "stability"=stability_score,
                        "constraint"=constraint_score,
                        "total"=total_score))
  extraction_type <- ifelse(missing(constraints), "(w/ stability-based extraction)",
                            ifelse(alpha == 0, "(w/ constraint-based extraction)", "(w/ mixed-objective extraction)"))
  substrs <- unlist(strsplit(x$method, split = " \\(w\\/"))
  out[["method"]] <- if (length(substrs) > 1) paste(substrs[[1]], extraction_type) else paste(out[["method"]], extraction_type)
  class(out) <- "hclust"
  return(list(cluster=attr(cl_tree, "cluster"), hc=out))
}
