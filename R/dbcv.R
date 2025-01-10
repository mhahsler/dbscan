#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2024 Michael Hahsler, Matt Piekenbrock

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


#' Density-Based Clustering Validation Index (DBCV)
#'
#' Calculate the Density-Based Clustering Validation Index (DBCV)  for a
#' clustering.
#'
#' DBCV (Moulavi et al, 2014) computes a score based on the density sparseness of each cluster
#' and the density separation of each pair of clusters.
#'
#' The density sparseness of a cluster (DSC) is deﬁned as the maximum edge weight of
#' a minimal spanning tree for the points of the cluster using the mutual
#' reachability distance based on the all-points-core-distance.
#'
#' The density separation of a pair of clusters (DSPC)
#' is deﬁned as the minimum reachability distance between the internal nodes of
#' the spanning trees of the two clusters.
#'
#' The validity index for a cluster is calculated using these measures and aggregated
#' to a validity index for the whole clustering using a weighted average.
#'
#' The index is in the range \eqn{[-1,1]}. If the cluster density compactness is better
#' than the density separation, a positive value is returned. The actual value depends
#' on the separability of the data. In general, greater values
#' of the measure indicating a better density-based clustering solution.
#'
#' Noise points are included in the calculation only in the weighted average,
#' therefore clustering with more noise points will get a lower index.
#'
#' **Performance note:** This implementation calculates a distance matrix and thus
#' can only be used for small or sampled datasets.
#'
#' @aliases dbcv DBCV
#' @family Evaluation Functions
#'
#' @param x a data matrix or a dist object.
#' @param cl a clustering (e.g., a integer vector)
#' @param d dimensionality of the original data if a dist object is provided.
#' @param metric distance metric used. The available metrics are the methods
#'        implemented by `dist()` plus `"sqeuclidean"` for the squared
#'        Euclidean distance used in the original DBCV implementation.
#' @param sample sample size used for large datasets.
#'
#' @return A list with the DBCV `score` for the clustering,
#'   the density sparseness of cluster (`dsc`) values,
#'   the density separation of pairs of clusters (`dspc`) distances,
#'   and the validity indices of clusters (`c_c`).
#'
#' @author Matt Piekenbrock and Michael Hahsler
#' @references Davoud Moulavi and Pablo A. Jaskowiak and
#' Ricardo J. G. B. Campello and Arthur Zimek and Jörg Sander (2014).
#' Density-Based Clustering Validation. In
#' _Proceedings of the 2014 SIAM International Conference on Data Mining,_
#' pages 839-847
#' \doi{10.1137/1.9781611973440.96}
#'
#' Pablo A. Jaskowiak (2022). MATLAB implementation of DBCV.
#' \url{https://github.com/pajaskowiak/dbcv}
#' @examples
#' # Load a test dataset
#' data(Dataset_1)
#' x <- Dataset_1[, c("x", "y")]
#' class <- Dataset_1$class
#'
#' clplot(x, class)
#'
#' # We use MinPts 3 and use the knee at eps = .1 for dbscan
#' kNNdistplot(x, minPts = 3)
#'
#' cl <- dbscan(x, eps = .1, minPts = 3)
#' clplot(x, cl)
#'
#' dbcv(x, cl)
#'
#' # compare to the DBCV index on the original class labels and
#' # with a random partitioning
#' dbcv(x, class)
#' dbcv(x, sample(1:4, replace = TRUE, size = nrow(x)))
#'
#' # find the best eps using dbcv
#' eps_grid <- seq(.01,.2, by = .01)
#' cls <- lapply(eps_grid, FUN = function(e) dbscan(x, eps = e, minPts = 3))
#' dbcvs <- sapply(cls, FUN = function(cl) dbcv(x, cl)$score)
#'
#' plot(eps_grid, dbcvs, type = "l")
#'
#' eps_opt <- eps_grid[which.max(dbcvs)]
#' eps_opt
#'
#' cl <- dbscan(x, eps = eps_opt, minPts = 3)
#' clplot(x, cl)
#' @export
dbcv <- function(x, cl, d, metric = "euclidean", sample = NULL) {
  # a clustering with a cluster element
  if (is.list(cl)) {
    cl <- cl$cluster
  }

  if (inherits(x, "dist")) {
    xdist <- x
    if (missing(d))
      stop("d needs to be specified if a distance matrix is supplied!")

  } else if (.matrixlike(x)) {
    if (!is.null(sample)) {
      take <- sample(nrow(x), size = sample)
      x <- x[take, ]
      cl <- cl[take]
    }

    x <- as.matrix(x)
    if (!missing(d) && d != ncol(x))
      stop("d does not match the number of columns in x!")
    d <- ncol(x)

    if (pmatch(metric, "sqeuclidean"))
      xdist <- dist(x, method = "euclidean")^2
    else
      xdist <- dist(x, method = metric)

  } else
    stop("'dbcv' expects x needs to be a matrix to calculate distances.")

  .check_dist(xdist)
  n <- attr(xdist, "Size")

  # in case we get a factor
  cl <- as.integer(cl)

  if (length(cl) != n)
    stop("cl does not match the number of rows in x!")

  ## getClusterIdList removes noise points and singleton clusters.
  ## returns indices reorder by cluster
  cl_ids_idx <- getClusterIdList(cl)
  all_cl_ids <- unlist(cl_ids_idx)
  n_cl <- length(cl_ids_idx)

  ## reordered distances w/o noise
  cl_dist <- dist_subset(xdist, all_cl_ids)

  ## 1. Calculate all-points-core-distance
  ## Calculate the all-points-core-distance for each point, within each cluster
  ## Note: this needs the dimensionality of the data d
  all_pts_core_dist <- lapply(
    cl_ids_idx,
    FUN = function(ids) {
      dists <- (rowSums(as.matrix((
        1 / dist_subset(xdist, ids)
      ) ^ d)) / (length(ids) - 1)) ^ (-1 / d)
      #dists <- ((rowSums(1 / as.matrix(xdist)[ids, ids]) ^ d) / (length(ids) - 1)) ^ (-1 / d)
    }
  )
  all_pts_core_dist <- unlist(all_pts_core_dist)


  ## 2. Create for each cluster a mutual reachability MSTs based on
  ## the all_pts_core_dist
  cl_mrd <- structure(mrd(cl_dist, all_pts_core_dist),
                      class = "dist",
                      Size = length(all_cl_ids))
  ## Noise points are removed, but the index is affected by dividing by the
  ## total number of objects including the noise points (n)!

  ## mst is a matrix with columns: from to and weight
  mrd_graphs <- lapply(cl_ids_idx, function(idx) {
    # is this necessary? cl_ids_idx is already in the same order as all_cl_ids
    #rel_idx <- match(idx, all_cl_ids)
    rel_idx <- idx
    mst(x_dist = dist_subset(cl_mrd, rel_idx),
                     n = length(rel_idx))
  })

  ## 3. Density Sparseness of a Cluster (DSC):
  ## The maximum edge weight of the internal edges in the cluster's
  ## mutual reachability MST.

  ## find internal nodes for DSC and DSPC. Internal nodes have a degree > 1
  internal_nodes <- lapply(mrd_graphs, function(mst) {
    node_deg <- table(c(mst[, 1], mst[, 2]))
    idx <- as.integer(names(node_deg)[which(node_deg > 1)])
    idx
  })

  dsc <- mapply(function(mst, int_idx) {
    # find internal edges
    int_edge_idx <- (mst[, 1] %in% int_idx) &
      (mst[, 2] %in% int_idx)
    if (length(which(int_edge_idx)) == 0) {
      return(as.numeric(max(mst[, 3])))
    }
    as.numeric(max(mst[int_edge_idx, 3]))
  }, mrd_graphs, internal_nodes)


  ## 4. Density Separation of a Pair of Clusters (DSPC):
  ## The minimum reachability distance between the internal nodes of the
  ## internal nodes of a pair of MST_MRD's of clusters Ci and Cj
  dspc_dist <- dspc(cl_ids_idx, internal_nodes, unlist(cl_ids_idx), cl_mrd)
  # returns a matrix with Ci, Cj, dist

  # make it into a full distance matrix
  dspc_dist <- dspc_dist[, 3L]
  class(dspc_dist) <- "dist"
  attr(dspc_dist, "Size") <- n_cl
  attr(dspc_dist, "Diag") <- FALSE
  attr(dspc_dist, "Upper") <- FALSE

  dspc_mm <- as.matrix(dspc_dist)
  diag(dspc_mm) <- NA

  ## 5. Validity index of a cluster:
  min_separation <- apply(dspc_mm, MARGIN = 1, min, na.rm = TRUE)
  v_c <- (min_separation - dsc) / pmax(min_separation, dsc)


  ## 5. Validity index for whole clustering
  res <- sum(lengths(cl_ids_idx) / n * v_c)

  return(list(
    score = res,
    n = n,
    n_c = lengths(cl_ids_idx),
    d = d,
    dsc = dsc,
    dspc = dspc_dist,
    v_c = v_c
  ))
}


getClusterIdList <- function(cl) {
  ## In DBCV, singletons are ambiguously defined. However, they cannot be
  ## considered valid clusters, for reasons listed in section 4 of the
  ## original paper. To ensure coverage, they are assigned into the noise category.
  cl_freq <- table(cl)
  cl[which(cl %in% as.integer(names(cl_freq[cl_freq == 1])))] <- 0L
  if (all(cl == 0)) {
    return(0)
  }

  cl_ids <- unique(cl)            # all cluster ids
  cl_valid <- cl_ids[cl_ids != 0] # valid cluster indices (non-noise)
  n_cl <- length(cl_valid)        # number of clusters

  ## 1 or 0 clusters results in worst score + a warning
  if (n_cl <= 1) {
    warning("DBCV is undefined for less than 2 non-noise clusters passed in.")
    return(-1L)
  }

  ## Indexes
  cl_ids_idx <- lapply(cl_valid, function(id)
    sort(which(cl == id))) ## the sort is important for indexing purposes
  return(cl_ids_idx)
}
