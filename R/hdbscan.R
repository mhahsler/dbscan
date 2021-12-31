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
#' This fast implementation of HDBSCAN (Hahsler et al, 2019) computes the
#' hierarchical cluster tree representing density estimates along with the
#' stability-based flat cluster extraction proposed by Campello et al. (2013).
#' HDBSCAN essentially computes the hierarchy of all DBSCAN* clusterings, and
#' then uses a stability-based extraction method to find optimal cuts in the
#' hierarchy, thus producing a flat solution.
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
#' [predict()] assigns each new data point to the same cluster as the point
#' in the original data that has the smallest mutual reachability distance.
#'
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
#' @param ...  additional arguments are passed on.
#' @param scale integer; used to scale condensed tree based on the graphics
#' device. Lower scale results in wider trees.
#' @param gradient character vector; the colors to build the condensed tree
#' coloring with.
#' @param show_flat logical; whether to draw boxes indicating the most stable
#' clusters.
#'
#' @return A object of class `hdbscan` with the following components:
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
#' @author Matt Piekenbrock
#' @references Hahsler M, Piekenbrock M, Doran D (2019). dbscan: Fast
#' Density-Based Clustering with R.  _Journal of Statistical Software_,
#' 91(1), 1-30.
#' \doi{10.18637/jss.v091.i01}
#'
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
#'
#' ## Predict cluster for new data points.
#' newdata <- data.frame(
#'   X = runif(n = 10, min(DS3[, 1]), max(DS3[, 1])),
#'   Y = runif(n = 10, min(DS3[, 2]), max(DS3[, 2]))
#'   )
#'
#' plot(DS3, col = res$cluster + 1L, pch = ".")
#' points(newdata, col = "red", pch = 4, lwd = 2)
#'
#' pred_label <- predict(res, newdata, data = DS3)
#' pred_label
#'
#' points(newdata, col = pred_label + 1L, pch = 1, lwd = 2, cex = 2)
#' @export hdbscan
hdbscan <- function(x,
  minPts,
  gen_hdbscan_tree = FALSE,
  gen_simplified_tree = FALSE) {
  if (.matrixlike(x) && !inherits(x, "dist")) {
    x <- as.matrix(x)
    if (!is.numeric(x))
      stop("hdbscan expects numerical data")
    xdist <- dist(x, method = "euclidean")
  } else if (inherits(x, "dist")) {
    ## this is for non-euclidean distances (note: this is slower for kNNdist)
    xdist <- x
  } else{
    stop(
      "hdbscan expects a matrix-coercible object of numerical data, and xdist to be a 'dist' object (or not supplied)."
    )
  }

  core_dist <- kNNdist(x, k = minPts - 1)

  ## At this point, xdist should be a dist object.
  n <- attr(xdist, "Size")

  ## Mutual reachability matrix
  mrd <- mrd(xdist, core_dist)

  ## Get MST, convert to RSL representation
  mst <- prims(mrd, n)
  hc <- hclustMergeOrder(mst, order(mst[, 3]))
  hc$call <- match.call()

  ## Process the hierarchy to retrieve all the necessary info needed by HDBSCAN
  res <- computeStability(hc, minPts, compute_glosh = TRUE)
  res <- extractUnsupervised(res)
  cl <- attr(res, "cluster")
  sl <- attr(res, "salient_clusters")

  ## Generate membership 'probabilities' using core distance as the measure of density
  prob <- rep(0, length(cl))
  for (cid in sl) {
    ccl <- res[[as.character(cid)]]
    max_f <- max(core_dist[which(cl == cid)])
    pr <- (max_f - core_dist[which(cl == cid)]) / max_f
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
    "Available fields: ", paste(names(x), collapse = ", ")
  ), exdent = 18))
}

#' @rdname hdbscan
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
    col_breaks <- seq(0, length(x$cluster) + nclusters, by = nclusters)
    gcolors <-
      grDevices::colorRampPalette(gradient)(length(col_breaks))

    ## Depth-first search to recursively plot rectangles
    eps_dfs <- function(dend, index, parent_height, scale) {
      coord <- coords[index, ]
      cl_key <- as.character(attr(dend, "label"))

      ## widths == number of points in the cluster at each eps it was alive
      widths <-
        sapply(sort(hd_info[[cl_key]]$eps, decreasing = T), function(eps)
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
            sort(hd_info[[cl_key]]$eps, decreasing = T),
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
