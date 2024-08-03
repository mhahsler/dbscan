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

#' Hierarchical DBSCAN (HDBSCAN)
#'
#' Fast C++ implementation of the HDBSCAN (Hierarchical DBSCAN) and its related
#' algorithms.
#'
#' This fast implementation of HDBSCAN (Campello et al., 2013) computes the
#' hierarchical cluster tree representing density estimates along with the
#' stability-based flat cluster extraction. HDBSCAN essentially computes the
#' hierarchy of all DBSCAN* clusterings, and
#' then uses a stability-based extraction method to find optimal cuts in the
#' hierarchy, thus producing a flat solution.
#'
#' HDBSCAN performs the following steps:
#'
#' 1. Compute mutual reachability distance mrd between points
#'    (based on distances and core distances).
#' 2. Use mdr as a distance measure to construct a minimum spanning tree.
#' 3. Prune the tree using stability.
#' 4. Extract the clusters.
#'
#' Additional, related algorithms including the "Global-Local Outlier Score
#' from Hierarchies" (GLOSH; see section 6 of Campello et al., 2015)
#' is available in function [glosh()]
#' and the ability to cluster based on instance-level constraints (see
#' section 5.3 of Campello et al. 2015) are supported. The algorithms only need
#' the parameter `minPts`.
#'
#' Note that `minPts` not only acts as a minimum cluster size to detect,
#' but also as a "smoothing" factor of the density estimates implicitly
#' computed from HDBSCAN.
#'
#' When using the optional parameter `cluster_selection_epsilon`,
#' a combination between DBSCAN* and HDBSCAN* can be achieved
#' (see Malzer & Baum 2020). This means that part of the
#' tree is affected by `cluster_selection_epsilon` as if
#' running DBSCAN* with `eps` = `cluster_selection_epsilon`.
#' The remaining part (on levels above the threshold) is still
#' processed by HDBSCAN*'s stability-based selection algorithm
#' and can therefore return clusters of variable densities.
#' Note that there is not always a remaining part, especially if
#' the parameter value is chosen too large, or if there aren't
#' enough clusters of variable densities. In this case, the result
#' will be equal to DBSCAN*.
# `cluster_selection_epsilon` is especially useful for cases
#' where HDBSCAN* produces too many small clusters that
#' need to be merged, while still being able to extract clusters
#' of variable densities at higher levels.
#'
#' `coredist()`: The core distance is defined for each point as
#' the distance to the `MinPts`'s neighbor. It is a density estimate.
#'
#' `mrdist()`: The mutual reachability distance is defined between two points as
#' `mrd(a, b) = max(coredist(a), coredist(b), dist(a, b))`. This distance metric is used by
#' HDBSCAN. It has the effect of increasing distances in low density areas.
#'
#' `predict()` assigns each new data point to the same cluster as the nearest point
#' if it is not more than that points core distance away. Otherwise the new point
#' is classified as a noise point (i.e., cluster ID 0).
#' @aliases hdbscan HDBSCAN print.hdbscan
#'
#' @family HDBSCAN functions
#' @family clustering functions
#'
#' @param x a data matrix (Euclidean distances are used) or a [dist] object
#' calculated with an arbitrary distance metric.
#' @param minPts integer; Minimum size of clusters. See details.
#' @param cluster_selection_epsilon double; a distance threshold below which
#  no clusters should be selected (see Malzer & Baum 2020)
#' @param gen_hdbscan_tree logical; should the robust single linkage tree be
#' explicitly computed (see cluster tree in Chaudhuri et al, 2010).
#' @param gen_simplified_tree logical; should the simplified hierarchy be
#' explicitly computed (see Campello et al, 2013).
#' @param verbose report progress.
#' @param ...  additional arguments are passed on.
#' @param scale integer; used to scale condensed tree based on the graphics
#' device. Lower scale results in wider trees.
#' @param gradient character vector; the colors to build the condensed tree
#' coloring with.
#' @param show_flat logical; whether to draw boxes indicating the most stable
#' clusters.
#' @param coredist numeric vector with precomputed core distances (optional).
#'
#' @return `hdbscan()` returns object of class `hdbscan` with the following components:
#' \item{cluster }{A integer vector with cluster assignments. Zero indicates
#' noise points.}
#' \item{minPts }{ value of the `minPts` parameter.}
#' \item{cluster_scores }{The sum of the stability scores for each salient
#' (flat) cluster. Corresponds to cluster IDs given the in `"cluster"` element.
#' }
#' \item{membership_prob }{The probability or individual stability of a
#' point within its clusters. Between 0 and 1.}
#' \item{outlier_scores }{The GLOSH outlier score of each point. }
#' \item{hc }{An [hclust] object of the HDBSCAN hierarchy. }
#'
#' `coredist()` returns a vector with the core distance for each data point.
#'
#' `mrdist()` returns a [dist] object containing pairwise mutual reachability distances.
#'
#' @author Matt Piekenbrock
#' @author Claudia Malzer (added cluster_selection_epsilon)
#'
#' @references
#' Campello RJGB, Moulavi D, Sander J (2013). Density-Based Clustering Based on
#' Hierarchical Density Estimates. Proceedings of the 17th Pacific-Asia
#' Conference on Knowledge Discovery in Databases, PAKDD 2013, _Lecture Notes
#' in Computer Science_ 7819, p. 160.
#' \doi{10.1007/978-3-642-37456-2_14}
#'
#' Campello RJGB, Moulavi D, Zimek A, Sander J (2015). Hierarchical density
#' estimates for data clustering, visualization, and outlier detection.
#' _ACM Transactions on Knowledge Discovery from Data (TKDD),_ 10(5):1-51.
#' \doi{10.1145/2733381}
#'
#' Malzer, C., & Baum, M. (2020). A Hybrid Approach To Hierarchical
#' Density-based Cluster Selection.
#' In 2020 IEEE International Conference on Multisensor Fusion
#' and Integration for Intelligent Systems (MFI), pp. 223-228.
#' \doi{10.1109/MFI49285.2020.9235263}
#' @keywords model clustering hierarchical
#' @examples
#' ## cluster the moons data set with HDBSCAN
#' data(moons)
#'
#' res <- hdbscan(moons, minPts = 5)
#' res
#'
#' plot(res)
#' plot(moons, col = res$cluster + 1L)
#'
#' ## cluster the moons data set with HDBSCAN using Manhattan distances
#' res <- hdbscan(dist(moons, method = "manhattan"), minPts = 5)
#' plot(res)
#' plot(moons, col = res$cluster + 1L)
#'
#' ## DS3 from Chameleon
#' data("DS3")
#'
#' res <- hdbscan(DS3, minPts = 50)
#' res
#'
#' ## Plot the simplified tree, highlight the most stable clusters
#' plot(res, show_flat = TRUE)
#'
#' ## Plot the actual clusters (noise has cluster id 0 and is shown in black)
#' plot(DS3, col = res$cluster + 1L, cex = .5)
#'
#' ## Example for HDBSCAN(e) using cluster_selection_epsilon
#' # data with clusters of various densities.
#' X <- data.frame(
#'  x = c(
#'   0.08, 0.46, 0.46, 2.95, 3.50, 1.49, 6.89, 6.87, 0.21, 0.15,
#'   0.15, 0.39, 0.80, 0.80, 0.37, 3.63, 0.35, 0.30, 0.64, 0.59, 1.20, 1.22,
#'   1.42, 0.95, 2.70, 6.36, 6.36, 6.36, 6.60, 0.04, 0.71, 0.57, 0.24, 0.24,
#'   0.04, 0.04, 1.35, 0.82, 1.04, 0.62, 0.26, 5.98, 1.67, 1.67, 0.48, 0.15,
#'   6.67, 6.67, 1.20, 0.21, 3.99, 0.12, 0.19, 0.15, 6.96, 0.26, 0.08, 0.30,
#'   1.04, 1.04, 1.04, 0.62, 0.04, 0.04, 0.04, 0.82, 0.82, 1.29, 1.35, 0.46,
#'   0.46, 0.04, 0.04, 5.98, 5.98, 6.87, 0.37, 6.47, 6.47, 6.47, 6.67, 0.30,
#'   1.49, 3.21, 3.21, 0.75, 0.75, 0.46, 0.46, 0.46, 0.46, 3.63, 0.39, 3.65,
#'   4.09, 4.01, 3.36, 1.43, 3.28, 5.94, 6.35, 6.87, 5.60, 5.99, 0.12, 0.00,
#'   0.32, 0.39, 0.00, 1.63, 1.36, 5.67, 5.60, 5.79, 1.10, 2.99, 0.39, 0.18
#'   ),
#'  y = c(
#'   7.41, 8.01, 8.01, 5.44, 7.11, 7.13, 1.83, 1.83, 8.22, 8.08,
#'   8.08, 7.20, 7.83, 7.83, 8.29, 5.99, 8.32, 8.22, 7.38, 7.69, 8.22, 7.31,
#'   8.25, 8.39, 6.34, 0.16, 0.16, 0.16, 1.66, 7.55, 7.90, 8.18, 8.32, 8.32,
#'   7.97, 7.97, 8.15, 8.43, 7.83, 8.32, 8.29, 1.03, 7.27, 7.27, 8.08, 7.27,
#'   0.79, 0.79, 8.22, 7.73, 6.62, 7.62, 8.39, 8.36, 1.73, 8.29, 8.04, 8.22,
#'   7.83, 7.83, 7.83, 8.32, 8.11, 7.69, 7.55, 7.20, 7.20, 8.01, 8.15, 7.55,
#'   7.55, 7.97, 7.97, 1.03, 1.03, 1.24, 7.20, 0.47, 0.47, 0.47, 0.79, 8.22,
#'   7.13, 6.48, 6.48, 7.10, 7.10, 8.01, 8.01, 8.01, 8.01, 5.99, 8.04, 5.22,
#'   5.82, 5.14, 4.81, 7.62, 5.73, 0.55, 1.31, 0.05, 0.95, 1.59, 7.99, 7.48,
#'   8.38, 7.12, 2.01, 1.40, 0.00, 9.69, 9.47, 9.25, 2.63, 6.89, 0.56, 3.11
#'  )
#' )
#'
#' ## HDBSCAN splits one cluster
#' hdb <- hdbscan(X, minPts = 3)
#' plot(hdb, show_flat = TRUE)
#' hullplot(X, hdb, main = "HDBSCAN")
#'
#' ## DBSCAN* marks the least dense cluster as outliers
#' db <- dbscan(X, eps = 1, minPts = 3, borderPoints = FALSE)
#' hullplot(X, db, main = "DBSCAN*")
#'
#' ## HDBSCAN(e) mixes HDBSCAN AND DBSCAN* to find all clusters
#' hdbe <- hdbscan(X, minPts = 3, cluster_selection_epsilon = 1)
#' plot(hdbe, show_flat = TRUE)
#' hullplot(X, hdbe, main = "HDBSCAN(e)")
#' @export
hdbscan <- function(
    x,
    minPts,
    cluster_selection_epsilon = 0.0,
    gen_hdbscan_tree = FALSE,
    gen_simplified_tree = FALSE,
    verbose = FALSE) {
  if (!inherits(x, "dist") && !.matrixlike(x)) {
    stop("hdbscan expects a numeric matrix or a dist object.")
  }

  ## 1. Calculate the mutual reachability between points
  if (verbose) {
    cat("Calculating core distances...\n")
  }
  coredist <- coredist(x, minPts)
  if (verbose) {
    cat("Calculating the mutual reachability matrix distances...\n")
  }
  mrd <- mrdist(x, minPts,
    coredist = coredist
  )
  n <- attr(mrd, "Size")

  ## 2. Construct a minimum spanning tree and convert to RSL representation
  if (verbose) {
    cat("Constructing the minimum spanning tree...\n")
  }
  mst <- prims(mrd, n)
  hc <- hclustMergeOrder(mst, order(mst[, 3]))
  hc$call <- match.call()

  ## 3. Prune the tree
  ## Process the hierarchy to retrieve all the necessary info needed by HDBSCAN
  if (verbose) {
    cat("Tree pruning...\n")
  }
  res <- computeStability(hc, minPts, compute_glosh = TRUE)
  res <- extractUnsupervised(res, cluster_selection_epsilon = cluster_selection_epsilon)
  cl <- attr(res, "cluster")

  ## 4. Extract the clusters
  if (verbose) {
    cat("Extract clusters...\n")
  }
  sl <- attr(res, "salient_clusters")

  ## Generate membership 'probabilities' using core distance as the measure of density
  prob <- rep(0, length(cl))
  for (cid in sl) {
    max_f <- max(coredist[which(cl == cid)])
    pr <- (max_f - coredist[which(cl == cid)]) / max_f
    prob[cl == cid] <- pr
  }

  ## Match cluster assignments to be incremental, with 0 representing noise
  if (any(cl == 0)) {
    cluster <- match(cl, c(0, sl)) - 1
  } else {
    cluster <- match(cl, sl)
  }
  cl_map <-
    structure(sl, names = unique(cluster[hc$order][cluster[hc$order] != 0]))

  ## Stability scores
  ## NOTE: These scores represent the stability scores -before- the hierarchy traversal
  cluster_scores <-
    sapply(sl, function(sl_cid) {
      res[[as.character(sl_cid)]]$stability
    })
  names(cluster_scores) <- names(cl_map)

  ## Return everything HDBSCAN does
  attr(res, "cl_map") <-
    cl_map # Mapping of hierarchical IDS to 'normalized' incremental ids
  out <- structure(
    list(
      cluster = cluster,
      minPts = minPts,
      coredist = coredist,
      cluster_scores = cluster_scores,
      # (Cluster-wide cumulative) Stability Scores
      membership_prob = prob,
      # Individual point membership probabilities
      outlier_scores = attr(res, "glosh"),
      # Outlier Scores
      hc = hc # Hclust object of MST (can be cut for quick assignments)
    ),
    class = "hdbscan",
    hdbscan = res
  ) # hdbscan attributes contains actual HDBSCAN hierarchy

  ## The trees don't need to be explicitly computed, but they may be useful if the user wants them
  if (gen_hdbscan_tree) {
    out$hdbscan_tree <- buildDendrogram(hc)
  }
  if (gen_simplified_tree) {
    out$simplified_tree <- simplifiedTree(res)
  }
  return(out)
}

#' @rdname hdbscan
#' @export
print.hdbscan <- function(x, ...) {
  writeLines(c(
    paste0("HDBSCAN clustering for ", nobs(x), " objects."),
    paste0("Parameters: minPts = ", x$minPts),
    paste0(
      "The clustering contains ",
      ncluster(x),
      " cluster(s) and ",
      nnoise(x),
      " noise points."
    )
  ))

  print(table(x$cluster))
  cat("\n")
  writeLines(strwrap(paste0(
    "Available fields: ", toString(names(x))
  ), exdent = 18))
}

#' @rdname hdbscan
#' @export
plot.hdbscan <-
  function(
      x,
      scale = "suggest",
      gradient = c("yellow", "red"),
      show_flat = FALSE,
      ...) {
    ## Logic checks
    if (!(scale == "suggest" ||
      scale > 0.0)) {
      stop("scale parameter must be greater than 0.")
    }

    ## Main information needed
    hd_info <- attr(x, "hdbscan")
    dend <- x$simplified_tree %||% simplifiedTree(hd_info)
    coords <-
      node_xy(hd_info, cl_hierarchy = attr(hd_info, "cl_hierarchy"))

    ## Variables to help setup the scaling of the plotting
    nclusters <- length(hd_info)
    npoints <- length(x$cluster)
    nleaves <-
      length(all_children(
        attr(hd_info, "cl_hierarchy"),
        key = 0,
        leaves_only = T
      ))
    scale <-
      ifelse(scale == "suggest", nclusters, nclusters / as.numeric(scale))

    ## Color variables
    col_breaks <-
      seq(0, length(x$cluster) + nclusters, by = nclusters)
    gcolors <-
      grDevices::colorRampPalette(gradient)(length(col_breaks))

    ## Depth-first search to recursively plot rectangles
    eps_dfs <- function(dend, index, parent_height, scale) {
      coord <- coords[index, ]
      cl_key <- as.character(attr(dend, "label"))

      ## widths == number of points in the cluster at each eps it was alive
      widths <-
        sapply(sort(hd_info[[cl_key]]$eps, decreasing = TRUE), function(eps) {
          length(which(hd_info[[cl_key]]$eps <= eps))
        })
      if (length(widths) > 0) {
        widths <-
          unlist(c(
            widths + hd_info[[cl_key]]$n_children,
            rep(hd_info[[cl_key]]$n_children, hd_info[[cl_key]]$n_children)
          ))
      } else {
        widths <-
          as.vector(rep(hd_info[[cl_key]]$n_children, hd_info[[cl_key]]$n_children))
      }

      ## Normalize and scale widths to length of x-axis
      normalize <- function(x) {
        (nleaves) * (x - 1) / (npoints - 1)
      }
      xleft <- coord[[1]] - normalize(widths) / scale
      xright <- coord[[1]] + normalize(widths) / scale

      ## Top is always parent height, bottom is when the points died
      ## Minor adjustment made if at the root equivalent to plot.dendrogram(edge.root=T)
      if (cl_key == "0") {
        ytop <-
          rep(
            hd_info[[cl_key]]$eps_birth + 0.0625 * hd_info[[cl_key]]$eps_birth,
            length(widths)
          )
        ybottom <- rep(hd_info[[cl_key]]$eps_death, length(widths))
      } else {
        ytop <- rep(parent_height, length(widths))
        ybottom <-
          c(
            sort(hd_info[[cl_key]]$eps, decreasing = TRUE),
            rep(hd_info[[cl_key]]$eps_death, hd_info[[cl_key]]$n_children)
          )
      }

      ## Draw the rectangles
      rect_color <-
        gcolors[.bincode(length(widths), breaks = col_breaks)]
      graphics::rect(
        xleft = xleft,
        xright = xright,
        ybottom = ybottom,
        ytop = ytop,
        col = rect_color,
        border = NA,
        lwd = 0
      )

      ## Highlight the most 'stable' clusters returned by the default flat cluster extraction
      if (show_flat) {
        salient_cl <- attr(hd_info, "salient_clusters")
        if (as.integer(attr(dend, "label")) %in% salient_cl) {
          x_adjust <-
            (max(xright) - min(xleft)) * 0.10 # 10% left/right border
          y_adjust <-
            (max(ytop) - min(ybottom)) * 0.025 # 2.5% above/below border
          graphics::rect(
            xleft = min(xleft) - x_adjust,
            xright = max(xright) + x_adjust,
            ybottom = min(ybottom) - y_adjust,
            ytop = max(ytop) + y_adjust,
            border = "red",
            lwd = 1
          )
          n_label <-
            names(which(attr(hd_info, "cl_map") == attr(dend, "label")))
          text(
            x = coord[[1]],
            y = min(ybottom),
            pos = 1,
            labels = n_label
          )
        }
      }

      ## Recurse in depth-first-manner
      if (is.leaf(dend)) {
        return(index)
      } else {
        left <-
          eps_dfs(
            dend[[1]],
            index = index + 1,
            parent_height = attr(dend, "height"),
            scale = scale
          )
        right <-
          eps_dfs(
            dend[[2]],
            index = left + 1,
            parent_height = attr(dend, "height"),
            scale = scale
          )
        return(right)
      }
    }

    ## Run the recursive plotting
    plot(
      dend,
      edge.root = TRUE,
      main = "HDBSCAN*",
      ylab = "eps value",
      leaflab = "none",
      ...
    )
    eps_dfs(dend,
      index = 1,
      parent_height = 0,
      scale = scale
    )
    return(invisible(x))
  }

#' @rdname hdbscan
#' @export
coredist <- function(x, minPts) {
  k <- minPts - 1
  kNN(x, k = k, sort = TRUE)$dist[, k]
}

#' @rdname hdbscan
#' @export
mrdist <- function(x, minPts, coredist = NULL) {
  if (inherits(x, "dist")) {
    .check_dist(x)
    x_dist <- x
  } else {
    x_dist <- dist(x,
      method = "euclidean",
      diag = FALSE,
      upper = FALSE
    )
  }

  if (is.null(coredist)) {
    coredist <- coredist(x, minPts)
  }
  mr_dist <- mrd(x_dist, coredist)

  class(mr_dist) <- "dist"
  attr(mr_dist, "Size") <- attr(x_dist, "Size")
  attr(mr_dist, "Diag") <- FALSE
  attr(mr_dist, "Upper") <- FALSE
  attr(mr_dist, "method") <- paste0(
    "mutual reachability (",
    attr(x_dist, "method"), ")"
  )
  mr_dist
}
