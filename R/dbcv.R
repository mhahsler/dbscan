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

dbcv <- function(x, cl, xdist = NULL){
  if (!.matrixlike(x)) { stop("'dbcv' expects x needs to be a matrix to calculate distances.") }
  if (!missing(xdist) && !is(xdist, 'dist')) { stop("'dbcv' expects xdist to be a dist object, if given.") }
  if (!missing(xdist) && attr(xdist, "method") != "euclidean") { stop("Only euclidean distance is supported.") }
  if (all(cl == 0)){ return(0) }
  
  ## Basic info
  xdist <- if (!missing(xdist)) xdist else dist(x, method = "euclidean")^2 # Note the use of squared euclidean
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
  all_pts_cd <- lapply(cl_ids_idx, function(cl_idx){
    cl_x <- x[cl_idx,] # cluster point coordinates
    cl_n <- nrow(cl_x) # number of points in the cluster
    # cl_knn <- kNN(cl_x, k = cl_n - 1) # KNN distances 
    # cl_knn_dist <- cl_knn$dist^2 # KNN dist (NOTE: squared euclidean)
    dist_cl <- as.matrix(dist(cl_x, method = "euclidean")^2)
    apply(dist_cl, 1, function(i_all_knn) {
      (sum((1/i_all_knn[i_all_knn != 0])^dim_x)/(cl_n-1))^(-1/dim_x) # Note the squared euclidean distance
    })
  })
  
  ## Given an original point index 'idx', retrieve it's all point core distance
  acp_dist_map <- structure(unlist(all_pts_cd), names = as.character(unlist(cl_ids_idx)))
  
  ## Mutual Reachability distances 
  ## To avoid computing the full mutual reachability matrix, cluster indices are used to get subsets 
  ## of the 'dist' object corresponding to the given clustering, and then that is combined with the all points core 
  ## distance to compute the resulting mutual reachability distances for each point
  mrd_dist <- mapply(function(cl_idx, cl_cd){ 
    idx <- combn(length(cl_idx), 2L)
    idx <- cbind(cl_idx[idx[1,]], cl_idx[idx[2,]])
    ## You can verify with all(dist(x[cl_idx,])^2 == xdist[(INDEX_TF(n, idx[, 1] - 1, idx[, 2] - 1))+1])
    mrd(dm = xdist[(INDEX_TF(n, idx[, 1] - 1, idx[, 2] - 1))+1], cl_cd) 
  }, cl_ids_idx, all_pts_cd, SIMPLIFY = FALSE)
  
  ## Mutual reachability MSTs
  mrd_graphs <- mapply(function(cl_idx, cl_mrd){ 
    mst <- prims(x_dist = cl_mrd, n = length(cl_idx))
    mst[order(mst[, 3]),] # return mst ordered by edge weight
  }, cl_ids_idx, mrd_dist, SIMPLIFY = FALSE)
  
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
  node_ids <- mapply(function(idx, i) cl_ids_idx[[i]][idx], internal_nodes, 1:n_cl, SIMPLIFY = FALSE)
  config <- list(n = n, ncl = n_cl, n_pairs = choose(n_cl, 2), node_ids = node_ids, acp = acp_dist_map, xdist = xdist)
  dspc_mat <- dspc(config = config) ## Computationally intensive

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