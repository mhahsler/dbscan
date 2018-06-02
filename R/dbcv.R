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
  if (is.numeric(x) && is.null(dim(x))){ x <- matrix(x) }
  if (!.matrixlike(x)) { stop("'dbcv' expects x needs to be a matrix to calculate distances.") }
  if (!missing(xdist) && !is(xdist, 'dist')) { stop("'dbcv' expects xdist to be a dist object, if given.") }
  if (!missing(xdist) && attr(xdist, "method") != "euclidean") { stop("Only euclidean distance is supported.") }
  
  ## Basic info
  if (missing(xdist) || is.null(xdist)) { 
    xdist <- dist(x, method = "euclidean")^(ifelse(squared, 2, 1)) ## use squared if flagged
  }
  n <- nrow(x)                    # number of coordinates
  dim_x <- ncol(x)                # dimensionality of space
  
  if (is.numeric(cl)){
    cl_ids_idx <- getClusterIdList(cl)
    all_pts_cd <- all_pts_core(x, cl_ids_idx, squared) ## Get the all-points-core-distance for each point, within each cluster
    res <- computeDBCV(x, xdist, cl_ids_idx, all_pts_cd, n, d = dim_x, squared = squared)
  } else if (is.list(cl)){    ## Multiple cluster version 
    x_dist_sorted <- apply(as.matrix(xdist), 1, sort) ## get sorted matrix-form of distances
    res <- vector(mode = "numeric", length = length(cl))
    for (i in 1:length(cl)){
      cl_ids_idx <- getClusterIdList(cl[[i]])
      all_pts_cd <- dbscan:::all_pts_core_sorted_dist(x_dist_sorted, cl = cl_ids_idx, d = dim_x, squared = squared)
      res[i] <- computeDBCV(x, xdist, cl_ids_idx, all_pts_cd, n, d = dim_x, squared = squared)
    }
  } else { stop("'dbcv' expects 'cl' to either be an integer vector representing cluster membership or a list of integer vectors.") }
  return(res)
}

computeDBCV <- function(x, xdist, cl_ids_idx, all_pts_cd, n, d, squared){
  n_cl <- length(cl_ids_idx)
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
    idx <- as.integer(names(node_deg)[which(node_deg > 1)])
    int_edge_idx <- (mst[, 1] %in% idx) & (mst[, 2] %in% idx)
    if (length(which(int_edge_idx)) == 0){
      # Exception case: if there are no internal edges (e.g. when there are 1-3 vertices in the cluster), density sparseness is undefined. 
      # Since this is a border case, assume the validation index isn't affected and treat all edges as internal.
      return(as.integer(names(node_deg)))
    } else {
      return(idx)
    }
  })
  
  ## Density Sparseness of a Cluster (DSC) [per cluster]
  dsc <- mapply(function(mst, int_idx){
    int_edge_idx <- (mst[, 1] %in% int_idx) & (mst[, 2] %in% int_idx)
    if (length(which(int_edge_idx)) == 0){
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

getClusterIdList <- function(cl){
  ## In DBCV, singletons are ambiguously defined. However, the cannot be considered valid clusters, for reasons 
  ## listed in section 4 of the original paper. To ensure coverage, they are assigned into the noise category. 
  cl_freq <- table(cl)
  cl[which(cl %in% as.integer(names(cl_freq[cl_freq == 1])))] <- 0L
  if (all(cl == 0)){ return(0) }
  
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
  return(cl_ids_idx)
}