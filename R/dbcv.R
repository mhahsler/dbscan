#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2018 Michael Hahsler, Matt Piekenbrock

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

## Macro to allow indexing into lower-triangular. Useful for indexing distances in 'dist' objects.
INDEX_TF <- function(N,to,from) (N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1) ## 0-based!

dbcv <- function(x, cl, xdist = NULL, squared=TRUE){
  if (!.matrixlike(x)) { stop("'dbcv' expects x needs to be a matrix to calculate distances.") }
  if (!missing(xdist) && !is(xdist, 'dist')) { stop("'dbcv' expects xdist to be a dist object, if given.") }
  if (!missing(xdist) && attr(xdist, "method") != "euclidean") { stop("Only euclidean distance is supported.") }
  
  ## In DBCV, singletons are ambiguously defined. However, the cannot be considered valid clusters, for reasons 
  ## listed in section 4 of the original paper. To ensure coverage, they are assigned into the noise category. 
  cl_freq <- table(cl)
  cl[which(cl %in% as.integer(names(cl_freq[cl_freq == 1])))] <- 0L
  if (all(cl == 0)){ return(0) }
  
  ## Basic info
  if (missing(xdist) || is.null(xdist)) { 
    xdist <- dist(x, method = "euclidean")^(ifelse(squared, 2, 1)) ## use squared if flagged
  }
  n <- nrow(x)                    # number of coordinates
  dim_x <- ncol(x)                # dimensionality of space
  cl_ids <- unique(cl)            # all cluster ids
  cl_valid <- cl_ids[cl_ids != 0] # valid cluster indices (non-noise)
  n_cl <- length(cl_valid)        # number of clusters
  
  ## 1 or 0 clusters results in worst score + a warning
  if (n_cl <= 1){
    warning("DBCV is undefined for less than 2 non-noise clusters passed in.")
    return(-1L) 
  } 
  
  ## Indexes
  cl_ids_idx <- lapply(cl_valid, function(id) sort(which(cl == id))) ## the sort is important for indexing purposes
  
  ## Get the all-points-core-distance for each point, within each cluster
  all_pts_cd <- all_pts_core(x, cl_ids_idx, squared)
  
  ## Given an original point index 'idx', retrieve it's all point core distance
  all_cl_ids <- unlist(cl_ids_idx)
  cl_dist <- dist_subset(xdist, all_cl_ids)
  cl_mrd <- structure(mrd(cl_dist, unlist(all_pts_cd)), class = "dist", Size = length(all_cl_ids))
  
  ## Mutual reachability MSTs
  mrd_graphs <- lapply(cl_ids_idx, function(idx){
    rel_idx <- match(idx, all_cl_ids)
    mst <- prims(x_dist = dist_subset(cl_mrd, rel_idx), n = length(rel_idx))
    matrix(mst[order(mst[, 3]),], ncol = 3) # return mst ordered by edge weight
  })
  
  ## Get the indices of the points making up the internal nodes of the MSTs
  internal_nodes <- lapply(mrd_graphs, function(mst){
    node_deg <- table(c(mst[, 1], mst[, 2]))
    as.integer(names(node_deg)[which(node_deg > 1)])
  })
  
  ## Density Sparseness of a Cluster (DSC) [per cluster]
  dsc <- mapply(function(mst, int_idx){
    int_edge_idx <- (mst[, 1] %in% int_idx) & (mst[, 2] %in% int_idx)
    if (length(which(int_edge_idx)) == 0){
      # Exception case: if there are no internal edges (e.g. when there are 1-3 vertices in the cluster), density sparseness is undefined. 
      # Since this is a border case, assume the validation index isn't affected and treat all edges as internal.
      return(as.numeric(max(mst[, 3]))) 
    } 
    as.numeric(max(mst[which(int_edge_idx), 3])) # largest edge weight between internal nodes
  }, mrd_graphs, internal_nodes)
  
  ## Density Separation of a Pair of Clusters (DSPC) [per cluster pair]
  ## Using the internal nodes only, extract the MRD subsets corresponding to all of the pairwise clusters 
  ## Algorithmically, this uses a fast MST to find shortest path between them
  dspc_mat <- dspc(cl_ids_idx, internal_nodes, unlist(cl_ids_idx), cl_mrd) ## Computationally intensive

  ## Validity index [Vc(C_i)]; callable for a given cluster
  v_c <- function(i) {
    idx <- which(i == dspc_mat[, 1] | i == dspc_mat[, 2]) ## all cluster pairs involving i
    min_separation <- min(dspc_mat[idx, 3]) 
    (min_separation - dsc[i])/(max(c(min_separation, dsc[i])))
  }
  
  ## Final Density-Based Cluster Validation score
  res <- sum(sapply(1:n_cl, function(i) { (length(cl_ids_idx[[i]])/n) * v_c(i) }))
  return(res)
}