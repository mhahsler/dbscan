#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2015 Michael Hahsler

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

optics <- function(x, eps, minPts = 5, eps_cl, xi, search = "kdtree",
  bucketSize = 10, splitRule = "suggest", approx = 0) {

  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if(is.na(search)) stop("Unknown NN search type!")

  ### dist search
  if(search == 3) {
    if(!is(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  ## for dist we provide the R code with a frNN list and no x
  frNN <- list()
  if(is(x, "dist")) {
    frNN <- frNN(x, eps)
    ## add self match and use C numbering
    frNN$id <- lapply(1:length(frNN$id),
      FUN = function(i) c(i-1L, frNN$id[[i]]-1L))
    frNN$dist <- lapply(1:length(frNN$dist),
      FUN = function(i) c(0, frNN$dist[[i]])^2)

    x <- matrix()
    storage.mode(x) <- "double"

  }else{
    if(!.matrixlike(x)) stop("x needs to be a matrix")
    ## make sure x is numeric
    x <- as.matrix(x)
    if(storage.mode(x) == "integer") storage.mode(x) <- "double"
    if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")
  }

  if(length(frNN) == 0 && any(is.na(x))) stop("data/distances cannot contain NAs for optics (with kd-tree)!")

  ret <- optics_int(as.matrix(x), as.double(eps), as.integer(minPts),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx), frNN)

  ret$minPts <- minPts
  ret$eps <- eps
  ret$eps_cl <- NA
  ret$xi <- NA
  class(ret) <- "optics"
  
  ### find clusters
  if(!missing(eps_cl)) ret <- optics_cut(ret, eps_cl)
  if(!missing(xi)) ret <- opticsXi(ret, xi)
    
  ret
}

### extract clusters
optics_cut <- function(x, eps_cl) {
  reachdist <- x$reachdist[x$order]
  coredist <- x$coredist[x$order]
  n <- length(x$order)
  cluster <- integer(n)

  clusterid <- 0L         ### 0 is noise
  for(i in 1:n) {
    if(reachdist[i] > eps_cl) {
      if(coredist[i] <= eps_cl) {
        clusterid <- clusterid + 1L
        cluster[i] <- clusterid
      }else{
        cluster[i] <- 0L  ### noise
      }
    }else{
      cluster[i] <- clusterid
    }
  }

  x$eps_cl <- eps_cl
  ### fix the order so cluster is in the same order as the rows in x
  cluster[x$order] <- cluster
  x$cluster <- cluster

  x
}

print.optics <- function(x, ...) {
  cat("OPTICS clustering for ", length(x$order), " objects.", "\n", sep = "")
  cat("Parameters: ", "minPts = ", x$minPts,
    ", eps = ", x$eps,
    ", eps_cl = ", x$eps_cl,
    ", xi = ", x$xi, 
    "\n", sep = "")
  if(!is.null(x$cluster)) {
    cl <- unique(x$cluster)
    cl <- length(cl[cl!=0L])
    if(is.na(x$xi)) { 
      cat("The clustering contains ", cl, " cluster(s) and ",
          sum(x$cluster==0L), " noise points.",
          "\n", sep = "")
      print(table(x$cluster))
    } else {
      cat("The clustering contains ", nrow(x$clusters_xi), " cluster(s) and ",
          sum(x$cluster==0L), " noise points.",
          "\n", sep = "")
    }
    cat("\n")
  }
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}

plot.optics <- function(x, y=NULL, cluster = TRUE, dendrogram = FALSE, prop = 0.75, ...) {
    if(!is.null(x$cluster) && cluster && !dendrogram) {
      if(is.null(x$clusters_xi)) { 
        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
          ylab = "Reachability dist.", xlab = "OPTICS order",
          main = "Reachability Plot", ...)
      } else {
        y_max <- max(x$reachdist[which(x$reachdist != Inf)])
        # Sort clusters by size 
        hclusters <- x$clusters_xi[order(x$clusters_xi$end-x$clusters_xi$start),]
        colors <- palette() 
        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
             ylab = "Reachability dist.", xlab = NA, xaxt = "n",
             main = "Reachability Plot", yaxs="i", ylim=c(0,y_max), xaxt='n', ...)
        bottom_inches <- par("fin")[2] * par("plt")[3]
        middle_inches <- par("fin")[2] * abs(par("plt")[4] - par("plt")[3])
        y_increments <- ((bottom_inches/middle_inches)*y_max / nrow(hclusters) * prop)
        i <- 1:nrow(hclusters)
        segments(x0=hclusters$start[i], y0=-(y_increments*i), x1=hclusters$end[i], col=palette()[(i+1) %% 8 + 1], lwd=1, xpd=T)
      }
    } else if (dendrogram) {
      d_x <- as.dendrogram(x)
      if (!is.null(x$cluster) && cluster) {
        if(requireNamespace("dendextend", quietly = TRUE)) { # Use dendextend plotting abilities 
          if (!is.na(x$eps_cl)) { # Plotting according to epsilon threshold
            d_x <- dendextend::set(d_x, "labels_col", value = c(x$cluster + 1)[x$order])
          } else if (!is.na(x$xi)) { # Plotting according to xi threshold
            hier_asc <- x$clusters_xi[order(-(x$clusters_xi$end - x$clusters_xi$start)),]
            minimal <- sort(unique(extractXiClusters(x, T)$cluster))[-1]
            for(cl in 1:nrow(hier_asc)) {
              cl_indices <- unlist(hier_asc[cl,])
              cl_labels <- x$order[cl_indices[1]:cl_indices[2]]
              # rule <- if (cl_indices[3] %in% minimal) "all" else "any"
              d_x <- dendextend::branches_attr_by_labels(d_x, labels=cl_labels, type = "all", 
                                             TF_values = c(cl_indices[3]+1,Inf))
            }
          }
          plot(d_x)
        }
      } else { plot(d_x) }
    } else {
      plot(x$reachdist[x$order], type="h",
        ylab = "Reachability dist.", xlab = "OPTICS order",
        main = "Reachability Plot", ...)
    }
}

predict.optics <- function (object, data, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$cluster)

  nn <- frNN(rbind(data, newdata), eps = object$eps_cl,
    sort = TRUE)$id[-(1:nrow(data))]
  sapply(nn, function(x) {
    x <- x[x<=nrow(data)]
    x <- object$cluster[x][x>0][1]
    x[is.na(x)] <- 0L
    x
  })
}

# Add new Generic for other class extensions support reachability plotting 
as.reachability <- function(x, ...) UseMethod("as.reachability")

as.dendrogram.optics <- function(x) {
  if(!"optics" %in% class(x)) stop("'as.dendrogram' from the dbscan package requires an OPTICS object.")
  if(x$minPts > length(x$order)) { stop("'minPts' should be less or equal to the points in the dataset.") } 
  if(length(which(x$reachdist == Inf)) > 1) stop(cat("The current eps value is not large enough to capture the complete hiearchical structure of the dataset.\nPlease use a large eps value (such as Inf)."))
    
  # Start with 0-height leaves 
  leaves <- lapply(1:length(x$order), function(cid) {
    structure(x$order[cid], members=1, height=0, label=as.character(x$order[cid]),
        leaf=T, class="dendrogram")
  })
  
  # Get sorted reachability distances and indices (ascending)
  pl <- sort(x$reachdist)  
  pl_order <- order(x$reachdist) 
  
  # For all the reachability distances excluding Inf 
  for(i in which(pl != Inf)) {
    # Get the index of the point with next smallest reach dist and its neighbor
    p <- pl_order[i] 
    q <- x$order[which(x$order == p) - 1]
    
    # Get the actual index fo the branch(es) containing the p and q 
    hier_indices <- lapply(leaves, labels)
    p_leaf <- which(unlist(lapply(hier_indices, function(b) p %in% b)))
    q_leaf <- which(unlist(lapply(hier_indices, function(b) q %in% b)))
    
    # Create the new branch containing the leaves that were removed (midpoint filled in later)
    children <- list(leaves[[q_leaf]], leaves[[p_leaf]])
    cl_members <- as.integer(sum(sapply(children, function(b) attr(b, which="members"))))
    N <- structure(children, members=cl_members, midpoint=NA, 
                   height=pl[i], class="dendrogram")
    
    leaves <- append(leaves[-c(p_leaf, q_leaf)], list(N))
  }
  
  # Always shift down one; rev actually populates the midpoint calculation 
  return(rev(rev(leaves[[1]])))
}

as.reachability.dendrogram <- function(x) {
  if (!inherits(x, "dendrogram")) stop("The as.reachability method requires a dendrogram object.")
  DFS <- function(d, rp, pnode = 0, stack=c()) {
    if(is.null(attr(d, "height"))) {
      return(pnode)
    } else {
      stack <- c(attr(d, "height"), stack)
      pnode <- DFS(d[[1]], rp, pnode, stack)
      if (length(d) > 1) { 
        for (sub_b in 2:length(d))  { pnode <- DFS(d[[sub_b]], rp, pnode, stack) }
      }
      # If at a leaf node, compare to previous node
      if ((!is.null(attr(d, "leaf"))) && attr(d, "leaf")) { 
        rp[[labels(d)]] <- stack # Record the ancestors reachability values 
        if(is.null(rp[[as.character(pnode)]])) { # 1st point 
          new_reach <- Inf
        } else { # Smallest Common Ancestor 
          new_reach <- min(intersect(rp[[as.character(pnode)]][-1], stack[-1])) # Exclude height 0
        } 
        rp[["reachdist"]] <- c(rp[["reachdist"]], new_reach) 
        rp[["order"]] <- c(rp[["order"]], as.integer(labels(d)))
        return(as.integer(labels(d)))
      }
      return(pnode)
    }
  }
  # Keep track using environment to save progress by reference 
  rp <- new.env(parent = emptyenv())
  DFS(x, rp=rp)
  return(list(reachdist=rp[["reachdist"]], order=rp[["order"]]))
}

# getColors <- function(x) { 
#   if (!inherits(x, "optics")) { stop("'getColors' is only meant ot be used with OPTICS objects.") }
#   if (is.null(x$clusters_xi)) { stop("'getColors' is only meant ot be used in the context of plotting with the Xi method") }
#   leaves <- dbscan:::extractXiClusters(x, minimum=T)$cluster
#   leaf_i <- which(x$clusters_xi$cluster_id %in% sort(unique(leaves))[-1])
#   minimal <- lapply(leaf_i, function(leaf) structure(leaf, height=0))
#   hier <- x$clusters_xi[-leaf_i, ]
#   hier <- hier[order(hier$end - hier$start), ]
#   for (i in 1:nrow(hier)) {
#     c_cluster <- hier[i,]$start:hier[i,]$end
#     new_leaves <- leaves
#     new_leaves[x$order[c_cluster]] <- hier[i,]$cluster_id
#     replaced_clusters <- setdiff(leaves, new_leaves)
#     max_child_height <- max(sapply(which(unlist(minimal) %in% replaced_clusters), function(ind) attr(minimal[[ind]], "height")))
#     minimal <<- append(minimal, list(structure(hier[i,]$cluster_id, height=max_child_height+1)))
#     leaves <- new_leaves
#   }
#   res <- data.frame(cid=unlist(minimal), height=sapply(minimal, function(b) attr(b, "height")))
#   res$color <- rep("#FFFFFFFF", nrow(res))
#   cheight <- max(res$height)
#   for ()
#   outer <- which(res$height %in% cheight)
#   inner <- which(res$height %in% (cheight-1))
#   cl_colors <- rainbow(length(outer) + length(inner))
#   res[inner, ]$color <- cl_colors[order(inner)]
# }

