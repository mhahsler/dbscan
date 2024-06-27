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
#' @export
hdbscan <- function(x,
  minPts,
  gen_hdbscan_tree = FALSE,
  gen_simplified_tree = FALSE,
  verbose = FALSE) {
  if (!inherits(x, "dist") && !.matrixlike(x))
    stop("hdbscan expects a numeric matrix or a dist object.")

  ## 1. Calculate the mutual reachability between points
  if (verbose)
    cat("Calculating core distances...\n")
  coredist <- coredist(x, minPts)
  if (verbose)
    cat("Calculating the mutual reachability matrix distances...\n")
  mrd <- mrdist(x, minPts, coredist
    = coredist)
  n <- attr(mrd, "Size")

  ## 2. Construct a minimum spanning tree and convert to RSL representation
  if (verbose)
    cat("Constructing the minimum spanning tree...\n")
  mst <- prims(mrd, n)
  hc <- hclustMergeOrder(mst, order(mst[, 3]))
  hc$call <- match.call()

  ## 3. Prune the tree
  ## Process the hierarchy to retrieve all the necessary info needed by HDBSCAN
  if (verbose)
    cat("Tree pruning...\n")
  res <- computeStability(hc, minPts, compute_glosh = TRUE)
  res <- extractUnsupervised(res)
  cl <- attr(res, "cluster")

  ## 4. Extract the clusters
  if (verbose)
    cat("Extract clusters...\n")
  sl <- attr(res, "salient_clusters")

  ## Generate membership 'probabilities' using core distance as the measure of density
  prob <- rep(0, length(cl))
  for (cid in sl) {
    ccl <- res[[as.character(cid)]]
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
    sapply(sl, function(sl_cid)
      res[[as.character(sl_cid)]]$stability)
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
    out$hdbscan_tree = buildDendrogram(hc)
  }
  if (gen_simplified_tree) {
    out$simplified_tree = simplifiedTree(res)
  }
  return(out)
}

#' @rdname hdbscan
#' @export
print.hdbscan <- function(x, ...) {
  cl <- unique(x$cluster)
  cl <- length(cl[cl != 0L])
  writeLines(c(
    paste0("HDBSCAN clustering for ", length(x$cluster), " objects."),
    paste0("Parameters: minPts = ", x$minPts),
    paste0(
      "The clustering contains ",
      cl,
      " cluster(s) and ",
      sum(x$cluster == 0L),
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
  function(x,
    scale = "suggest",
    gradient = c("yellow", "red"),
    show_flat = FALSE,
    ...) {
    ## Logic checks
    if (!(scale == "suggest" ||
        scale > 0.0))
      stop("scale parameter must be greater than 0.")

    ## Main information needed
    hd_info <- attr(x, "hdbscan")
    dend <-
      if (is.null(x$simplified_tree))
        simplifiedTree(hd_info)
    else
      x$simplified_tree
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
        sapply(sort(hd_info[[cl_key]]$eps, decreasing = TRUE), function(eps)
          length(which(hd_info[[cl_key]]$eps <= eps)))
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
      normalize <- function(x)
        (nleaves) * (x - 1) / (npoints - 1)
      xleft <- coord[[1]] - normalize(widths) / scale
      xright <- coord[[1]] + normalize(widths) / scale

      ## Top is always parent height, bottom is when the points died
      ## Minor adjustment made if at the root equivalent to plot.dendrogram(edge.root=T)
      if (cl_key == "0") {
        ytop <-
          rep(hd_info[[cl_key]]$eps_birth + 0.0625 * hd_info[[cl_key]]$eps_birth,
            length(widths))
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
      scale = scale)
    return(invisible(x))
  }

#' @rdname hdbscan
#' @export
coredist <- function(x, minPts) {
  k <- minPts - 1
  kNN(x, k = k , sort = TRUE)$dist[, k]
}

#' @rdname hdbscan
#' @export
mrdist <- function(x, minPts, coredist = NULL) {
  if (inherits(x, "dist")) {
    .check_dist(x)
    x_dist <- x
  } else{
    x_dist <- dist(x,
      method = "euclidean",
      diag = FALSE,
      upper = FALSE)
  }

  if (is.null(coredist))
    coredist <- coredist(x, minPts)
  mr_dist <- mrd(x_dist, coredist)

  class(mr_dist) <- "dist"
  attr(mr_dist, "Size") <- attr(x_dist, "Size")
  attr(mr_dist, "Diag") <- FALSE
  attr(mr_dist, "Upper") <- FALSE
  attr(mr_dist, "method") <- paste0("mutual reachability (",
                                    attr(x_dist, "method"),")")
  mr_dist
}
