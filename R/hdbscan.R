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

hdbscan <- function(x, minPts = 5, xdist=NA, gen_rsl_tree = F) {
  ## Calculate Core distance using kNN
  if (missing(xdist) && (is(x, "data.table") || is(x, "data.frame") || is(x, "matrix"))) {
    euc_dist <- dist(x, method = "euclidean")
    core_dist <- kNNdist(x, k = minPts - 1)[, minPts - 1]
    n <- nrow(x)
  } else if (!missing(xdist) && is(xdist, "dist")) {
    euc_dist <- xdist
    core_dist <- kNNdist(x, k = minPts - 1)[, minPts - 1]
    n <- nrow(x)
  } else {
    stop("hdbscan expects a matrix-conformable object for x")
  }

  ## Mutual Reachability matrix
  mrd <- mrd(euc_dist, core_dist)

  ## Get MST, convert to RSL representation
  mst <- dbscan:::prims(mrd, n)
  hc <- dbscan:::hclustMergeOrder(mst, order(mst[, 3]))

  ## Process the hierarchy to retrieve all the necessary info needed by HDBSCAN
  res <- hdbscan_fast(hc, minPts)
  cl <- attr(res, "cluster")
  sl <- attr(res, "salient_clusters")
  
  ## Generate membership 'probabilities' using core distance as the measure of density
  prob <- rep(0, length(cl))
  for (cid in sl){
    ccl <- res[[as.character(cid)]]
    max_f <- max(core_dist[which(cl == cid)])
    pr <- (max_f - core_dist[which(cl == cid)])/max_f
    prob[cl == cid] <- pr
  }

  ## Match cluster assignments to be incremental, with 0 representing noise
  if (any(cl == 0)) {
    cluster <- match(cl, c(0, sl))-1
  } else { cluster <- match(cl, sl) }

  ## Final stability scores
  ## NOTE: These scores represent the stability scores before the hierarchy traversal
  cluster_scores <- sapply(sl, function(cid) res[[as.character(cid)]]$score)

  ## Return just a vector of cluster assignments by default, else return everything
  if (gen_rsl_tree){
    out <- structure(list(cluster=cluster, minPts=minPts,
                          cluster_scores=cluster_scores, # Stability Scores
                          membership_prob=prob, # Individual point membership probabilities
                          outlier_scores=attr(res, "glosh"), # Outlier Scores
                          mst=mst, # Minimum Spanning Tree (for debugging)
                          hc=hc, # Hclust object of MST (for debugging / maybe useful to user?)
                          rsl_tree=dbscan:::buildDendrogram(hcl)
    ), class="hdbscan") # runtime contains all the information used to build HDBSCAN
  } else {
    out <- structure(list(cluster=cluster, minPts=minPts,
                          cluster_scores=cluster_scores, # Stability Scores
                          membership_prob=prob, # Individual point membership probabilities
                          outlier_scores=attr(res, "glosh"), # Outlier Scores
                          mst=mst, # Minimum Spanning Tree (for debugging)
                          hc=hc # Hclust object of MST (for debugging)
    ), class="hdbscan") # runtime contains all the information used to build HDBSCAN
  }
  return(out)
}

print.hdbscan <- function(x, ...) {
  cl <- unique(x$cluster)
  cl <- length(cl[cl!=0L])
  writeLines(c(
    paste0("HDBSCAN clustering for ", length(x$cluster), " objects."),
    paste0("Parameters: minPts = ", x$minPts),
    paste0("The clustering contains ", cl, " cluster(s) and ",
           sum(x$cluster==0L), " noise points.")
  ))
  
  print(table(x$cluster))
  cat("\n")
  
  writeLines(strwrap(paste0("Available fields: ", paste(names(x), collapse = ", ")), exdent = 18))
}

plot.hdbscan <- function(x, scale=20, gradient=c("yellow", "red"), show_flat = F, ...){
  dend <- stats:::midcache.dendrogram(x$condensed_tree)
  runtime <- attr(x, "runtime")
  yBottom <- ceiling(max(sapply(runtime[[".cluster_info"]], function(cl) cl$height))) + 1

  dbscan:::all_children(attr(res, "cl_hierarchy"), 0, F)
  plot(dend, ylim=c(yBottom, 0), edge.root = F, main="HDBSCAN*", ylab="Lambda value")
  coords <- dendextend::get_nodes_xy(dend)

  ## Variables to help setup plotting
  nclusters <- length(runtime[[".cid"]])
  col_breaks <- seq(0, length(cl$cluster)+nclusters, by=nclusters)
  gcolors <- colorRampPalette(gradient)(length(col_breaks))
  normalize <- function(x) (nclusters)*(x - 1) / (length(cl$cluster) - 1)

  ## Depth-first search to recursively plot rectangles
  eps_dfs <- function(dend, index, parent_height, scale){
    coord <- coords[index,]
    widths <- sapply(sort(attr(dend, "eps"), decreasing = T), function(eps) length(which(attr(dend, "eps") <= eps)))
    xleft <- coord[[1]] - normalize(widths)/scale
    xright <- coord[[1]] + normalize(widths)/scale
    ytop <- rep(parent_height, length(widths))
    ybottom <- sort(1/attr(dend, "eps"), decreasing = F)
    rect(xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop,
         col = gcolors[.bincode(length(attr(dend, "contains")), breaks = col_breaks)],
         border=NA, lwd=0)
    if (show_flat){
      salient_cl <- names(attr(x, "runtime")[[".cluster"]])
      if (attr(dend, "label") %in% salient_cl){
        max_width <- max(normalize(widths)/scale)
        rect(xleft = min(xleft) - max_width, xright = max(xright) + max_width,
             ybottom = max(ybottom) + max_width, ytop = min(ytop) - max_width,
             border = "red", lwd=1)
      }
    }
    if (is.leaf(dend)){
      return(dend)
    } else {
      for (i in 1:length(dend)){
        eps_dfs(dend[[i]], as.integer(attr(dend[[i]], "label"))+1, parent_height = coords[index, 2], scale=scale)
      }
    }
  }
  eps_dfs(dend, index = 1, parent_height = 0, scale = scale)
}
