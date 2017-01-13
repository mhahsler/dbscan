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

hdbscan <- function(x, minPts = 5) {
  ## Calculate Core distance using kNN
  if (is(x, "data.table") || is(x, "data.frame") || is(x, "matrix")) {
    euc_dist <- as.matrix(dist(x, method = "euclidean"))
    core_dist <- kNNdist(x, k = minPts - 1)[, minPts - 1]
  } else if (is(x, "dist")) {
    euc_dist <- as.matrix(x)
    core_dist <- x[, minPts - 1]
  }
  
  ## Mutual Reachability matrix
  mrd <- mrd(euc_dist, core_dist)
  
  ## Convert MST to Robust Single Linkage dendrogram
  ## NOTE: hclust conversion apparently doesn't support n-ary trees, so for now we use our own! 
  mst_spantree <- vegan::spantree(mrd)
  mst <- cbind(from=mst_spantree$kid, to=2:nrow(mrd), dist=mst_spantree$dist)
  mst <- mst[order(-mst[,3]),]
  
  ## Workaround 
  hcl <- as.hclust(mst_spantree)
  rsl_tree <- as.dendrogram(hcl)
  
  ## n-ary method
  # mst <- mst[order(mst[, 3]),]
  # mst <- cbind(mst[, 1:2] - 1, mst[, 3]) 
  # rsl_tree <- mst_to_dendrogram(mst)
  
  ## This store midpoints for nicer default dendrogram plotting. Computing the midpoints is very slow.
  ## NOTE: Midpoints are only partially implemented for non-binary trees
  rsl_tree <- stats:::midcache.dendrogram(rsl_tree, quiet = T)
  #hcl <- as.hclust(mst_spantree) 
  #rsl_tree$merge <- hcl$merge
  #rsl_tree$order <- hcl$order
  #rsl_tree$height <- mst[, 3]
  
    
  ## Adjust recursion stack depth if necessary 
  if (nrow(mrd) > 1000) { options(expressions = nrow(mrd) + round(nrow(mrd) * log(nrow(mrd), base = 2))) }

  ## Compute `condensed` hierarchy + lambda values
  store <- new.env(parent = emptyenv())
  store[[".id"]] <- 0
  store[["ctree"]] <- structure(list(), class="dendrogram") # do not remove class
  buildHDBSCAN(rsl_tree, store, minPts = minPts)
  store[["ctree"]] <- store[["ctree"]][[1]]
  
  ## Extract most saliently separable clusters
  store[[".cluster"]] <- list()
  computeSalientClusters(store, dend = store[["ctree"]])
  
  ## Unfold cluster assignments into vector
  cluster_assignments <- lapply(store[[".cluster"]], function(cl) as.vector(cl$contains))
  cluster <- rep(0, nrow(x))
  for (cid in ls(cluster_assignments)) {
    cluster[cluster_assignments[[cid]]] <- as.integer(cid)
  }
  
  ## "Normalize" cluster assignments to be incremental, with 0 representing noise 
  if (any(cluster == 0)) {
    cluster <- match(cluster, as.integer(ls(cluster_assignments))) - 1
  } else { cluster <- match(cluster, as.integer(ls(cluster_assignments))) }

  ## Generate normalized membership probabilities 
  prob <- rep(0, length(cluster))
  for (cid in ls(cluster_assignments)){
    cl <- store[[".cluster_info"]][[cid]]
    ## Below is what Python is doing 
    # pr <- 1 - ((max(1/cl$eps) - 1/cl$eps) / max(1/cl$eps))
    
    ## But this seems to make more sense to me
    pr <- (1/cl$eps - 1/cl$eps_birth)/(max(1/cl$eps) - 1/cl$eps_birth)
    prob[cl$contains] <- pr
  }
  
  ## Final stability scores 
  ## NOTE: These scores represent the stability scores before the hierarchy traversal
  cluster_scores <- sapply(ls(store[[".cluster_info"]]), function(cid) {
    cl <- store[[".cluster_info"]][[as.character(cid)]]
    sum( 1/cl$eps - 1/cl$eps_birth )
  })
  
  ## Return just a vector of cluster assignments by default, else return everything 
  structure(list(cluster=cluster, minPts=minPts, 
                 rsl_tree=rsl_tree, # Robust Single Linkage Tree
                 condensed_tree=store[["ctree"]], # Condensed tree
                 cluster_scores=cluster_scores, # Stability Scores
                 membership_prob=prob, # Individual point membership probabilities
                 mst=mst, # Minimum Spanning Tree (for debugging)
                 hcl=hcl, # Hclust object of MST (for debugging)
                 mrd=mrd # Mutual Reachability matrix (for debugging)
                 ), class="hdbscan", runtime=store) # runtime contains all the information used to build HDBSCAN
}

print.hdbscan <- function(x, ...) {
  cat("HDBSCAN clustering for ", length(x$cluster), " objects.", "\n", sep = "")
  cat("Parameters: minPts = ", x$minPts, "\n", sep = "")
  cl <- unique(x$cluster)
  cl <- length(cl[cl!=0L])
  cat("The clustering contains ", cl, " cluster(s) and ", sum(x$cluster==0L),
      " noise points.",
      "\n", sep = "")
  print(table(x$cluster))
  cat("\nAvailable fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}

plot.hdbscan <- function(x, scale=20, gradient=c("yellow", "red"), show_flat = F, ...){
  dend <- stats:::midcache.dendrogram(x$condensed_tree)
  runtime <- attr(x, "runtime")
  yBottom <- ceiling(max(sapply(runtime[[".cluster_info"]], function(cl) cl$height))) + 1
  
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

## Traverse the RSL hierarchy to recursively store distance values that points dropped off clusters at
## This function builds the condensed hierarchical tree, and records all of the meta-information needed 
## to compute the optimal flat clustering immediately after returning 
## Notation: 
## rsl == robust single linkage tree 
## store == environment for storing the condensed tree (ctree), allows updating by reference
## eps == delta == epsilon == DBSCAN eps distance 
## lambda == 1/epsilon
buildHDBSCAN <- function(rsl, store, minPts, cid = 0, path=1, parent_height = Inf) {

  ## If current ID doesn't exist, start a new one, NA's purposefully put in as a logical check
  ## members attribute filled in recursively 
  if (!cid %in% store[[".cid"]]){
    store[[".cid"]] <- append(store[[".cid"]], cid)
    members <- unlist(rsl)
    store[["ctree"]][[path]] <- structure(list(), members=NA, # number of leaves; TBD
                                children=NA,  # child cluster IDs; TBD
                                height=1/parent_height, # Dendrogram height
                                label=as.character(cid), # cluster label
                                eps_birth=parent_height, # distance cluster was formed
                                contains=members, # points cluster contains
                                eps=rep(NA, length(members)), # distance thresholds cluster existed
                                class="dendrogram") # dendrogram class required for path indexing to work
  }
  
  ## Always compute current cluster ID's recorded epsilon distances, and retrieve the point indices
  cdelta <- attr(store[["ctree"]][[path]], "eps")
  current_index_set <- attr(store[["ctree"]][[path]], "contains")
  
  ## Save all branch sizes
  branch_sizes <- sapply(1:length(rsl), function(i)  attr(rsl[[i]], "members"))
  
  ## All branches are noise. Record the height that that happened at, 
  ## and update dendrogram heights with max lambda value 
  if(all(branch_sizes < minPts)) {
    ## Get the indices of the noise branches
    noise <- sapply(1:length(rsl), function(i) which(current_index_set %in% unlist(rsl[[i]])))
    noise <- unlist(noise)
    
    ## Assign noise
    cdelta[noise] <- attr(rsl, "height")
    attr(store[["ctree"]][[path]], "eps") <- cdelta
    
    ## Adjust dendrogram stuff for plotting
    attr(store[["ctree"]][[path]], "leaf") <- T
    attr(store[["ctree"]][[path]], "height") <- max(1/attr(store[["ctree"]][[path]], "eps"))
    
    ## Compute (initial) stability score 
    attr(store[["ctree"]][[path]], "score") <- sum(1/attr(store[["ctree"]][[path]], "eps") - 1/attr(store[["ctree"]][[path]], "eps_birth"))
    
    ## Leaves should have 0 children and 1 member
    attr(store[["ctree"]][[path]], "children") <- NULL 
    attr(store[["ctree"]][[path]], "members") <- 1
    
    ## Return cluster ID 
    return(attr(store[["ctree"]][[path]], "label"))
  } ## Points are falling off the cluster as noise. Record distances and recurse.  
  else if (length(which(branch_sizes >= minPts)) == 1)  {
    noise_branches <- which(branch_sizes < minPts)
    noise <- sapply(noise_branches, function(i) which(current_index_set %in% unlist(rsl[[i]])))
    noise <- unlist(noise)
    
    cdelta[noise] <- attr(rsl, "height")
    attr(store[["ctree"]][[path]], "eps") <- cdelta
    
    ## Recurse; Keep going through the hierarchy for significant branches
    leaves <- c()
    for (b_i in which(branch_sizes >= minPts)) {
      leaves <- unique(c(leaves, buildHDBSCAN(rsl[[b_i]], store, minPts, cid, path=path, parent_height=attr(store[["ctree"]][[path]], "height"))))
    }
    
    ## Store immediate children when not at leaf node
    if (!is.leaf(store[["ctree"]][[path]])){
      attr(store[["ctree"]][[path]], "children") <- sapply(1:length(store[["ctree"]][[path]]), function(b_i) attr(store[["ctree"]][[c(path, b_i)]], "label"))
    }
    
    ## Respect dendrogram structure definition and store number of leaves as the member attribute
    attr(store[["ctree"]][[path]], "members") <- length(leaves)
    
    ## Return cluster leaves
    return(leaves)
  } ## Both clusters are big enough to be considerd a true split, so increment cluster IDs, save the
    ## distance that the current cluster dissolved at and recurse, making child clusters
    ## inherit the current ID as their parent cluster ID 
  else {
    ## Only recurse over true clusters
    true_clusters <- which(branch_sizes >= minPts)

    ## Last delta value this cluster was alive
    attr(store[["ctree"]][[path]], "eps_death") <- attr(rsl, "height")
    
    ## If there are any unmarked points, they're part of condensed child clusters yet to be 
    ## computed, so set as smallest possible value
    if (any(is.na(cdelta))){
      cdelta[which(is.na(cdelta))] <- attr(rsl, "height")
      attr(store[["ctree"]][[path]], "eps") <- cdelta
    }
    
    ## Compute initial stability score
    attr(store[["ctree"]][[path]], "score") <- sum(1/attr(store[["ctree"]][[path]], "eps") - 1/attr(store[["ctree"]][[path]], "eps_birth"))
    
    ## Change dendrogram heights
    attr(store[["ctree"]][[path]], "height") <- max(1/attr(store[["ctree"]][[path]], "eps"))
    
    ## Recursively make new branches; 
    leaves <- c()
    for (b_i in 1:length(true_clusters)) {
      store[[".id"]] <- store[[".id"]] + 1
      leaves <- unique(c(leaves, buildHDBSCAN(rsl[[true_clusters[b_i]]], store, minPts, store[[".id"]], c(path, b_i), attr(rsl, "height"))))
    }
    
    ## Store immediate children when not at leaf node
    if (!is.leaf(store[["ctree"]][[path]])){
      attr(store[["ctree"]][[path]], "children") <- sapply(1:length(store[["ctree"]][[path]]), function(b_i) attr(store[["ctree"]][[c(path, b_i)]], "label"))
    }
    
    ## Respect dendrogram structure definition and store number of leaves as the member attribute
    attr(store[["ctree"]][[path]], "members") <- length(leaves)
    return(leaves)
  }
  return(NULL)
}

## Computes minPts-saliently separated clusters by computing the scores
computeSalientClusters <- function(store, dend) {
  ## Having information on every condensed cluster is useful 
  if (is.null(store[[".cluster_info"]][[attr(dend, "label")]])){
    store[[".cluster_info"]][[attr(dend, "label")]] <- attributes(dend)
  }
  
  ## Return stability score if at a leaf node
  if(is.leaf(dend)) {
    cluster_score <- sum(1/attr(dend, "eps") - 1/attr(dend, "eps_birth"))
    store[[".cluster"]][[attr(dend, "label")]] <- list(label=attr(dend, "label"), contains=attr(dend, "contains"))
    return(cluster_score)
  }
  ## Compute/retrieve children scores
  c_scores <- sapply(1:length(dend), function(b_i) computeSalientClusters(store, dend[[b_i]]))
  if (attr(dend, "label") != "0"){ # Do not ever consider root node
    score <- attr(dend, "score")
    children <- attr(dend, "children") 
    if (score > sum(c_scores)) {
      # Remove children and add parent
      store[[".cluster"]] <- Filter(function(cl) !cl$label %in% children, store[[".cluster"]])
      store[[".cluster"]][[attr(dend, "label")]] <- list(label=attr(dend, "label"), contains=attr(dend, "contains"))
      return(score)
    } else { # Keep children, update parent score
      return(sum(c_scores))
    }
  }
  return(0)
}

## Still needs a bit of work 
as.matrix.hdbscan <- function(x, ...) {
  # res <- matrix(0, nrow = length(x$cluster), ncol=length(x$cluster)-1)
  # for (eps_i in 1:nrow(x$mst)){
  #   res[, (nrow(x$mst)+1)-eps_i] <- unname(dendextend:::cutree.dendrogram(x$rsl_tree, h=x$mst[eps_i, 3]))
  # }
  # rownames(res) <- unlist(x$rsl_tree)
  # colnames(res) <- rev(sprintf("%.3f", x$mst[, 3]))
  # res
  mst <- cl$mst[order(-cl$mst[, 3]),]
  cuts <- cutree(cl$hcl, h = mst[,3])
  core_dist <- diag(cl$mrd)
  res <- sapply(1:nrow(mst), function(i) {
    noise <- which(core_dist > mst[i,3])
    if (length(noise) > 0){
      cuts[noise, i] <- 0
      return(match(cuts[, i], sort(unique(cuts[, i]))) - 1)
    }
    return(match(cuts[, i], sort(unique(cuts[, i]))))
  })
}


rise <- function(cl){
  nclusters <- which(mapply(function(x, y) x < 0 && y < 0, cl$merge[,1], cl$merge[,2]))
  # res <- matrix(0, nrow=nrow(cl$merge) +1, ncol=nrow(cl$merge))
  cc <- 1
  cc_indices <- new.env(parent = emptyenv()) # or rather 
  cc_other <- rep(0, nrow(cl$merge) + 1)
  for (i in 1:nrow(cl$merge)){
    m <- cl$merge[i,]
    if (all(m < 0)){
      cc_other[(-m)] <- cc
      #cc_indices[[as.character(cc)]] <- append(cc_indices[as.character(cc)], -m)
      cc <- cc + 1
    } else if (any(m < 0)){
      j <- -(m[which(m < 0)])
      #cc_indices[[as.character(cc)]] <- append(cc_indices[as.character(cc)], -m)
      prev_step <- cl$merge[m[which(m > 0)],]
      true_cc <- cc_other[prev_step[which(prev_step < 0)]]
      cc_other[j] <- true_cc
    } else {
      
    }
  }
}  


#   
#   ## New list idea: instead of using table, use armadillo vectors + subtraction to determine 
#   ## online when the update specific point differences 
#   
#   # added_noise <- mapply(function(i, j){
#   #   ifelse(is.na(table(res[, j])["0"]), -Inf, table(res[, j])["0"]) > ifelse(is.na(table(res[, i])["0"]), 0, table(res[, i])["0"])
#   # }, 1:(ncol(res)-1), 2:ncol(res))
#   # 
#   # ## Get the epsilon indices that produced 'true' cluster splits
#   # true_splits <- which(mapply(function(i, j) {
#   #   length(table(res[res[, i] != 0, i])) < length(table(res[res[, j] != 0, j]))
#   # }, 1:(ncol(res)-1), 2:ncol(res))) + 1
#   # length_changes <- sapply(apply(res, 2, table), length)
#   # true_splits <- true_splits[true_splits < max(which(length_changes == max(length_changes)))]
#   # 
#   # for (n in 2:nrow(mst)){
#   #   assignments <- res[, n][cl$hcl$order]
#   #   cids <- unique(assignments)
#   #   
#   #   if (n %in% true_splits) { ## true split
#   #     sapply(as.character(cids), function(cid) {
#   #       store[[cid]]$contains <- cl$hcl$order[res[, n][cl$hcl$order] == as.numeric(cid)]
#   #       store[[cid]]$eps <- eps_values[n]
#   #     })
#   #   } else if (added_noise[n]){ ## cluster is shrinking
#   #     sapply(as.character(cids), function(cid) {
#   #       index <- match(cl$hcl$order[res[, n][cl$hcl$order] == as.numeric(cid)], store[[cid]]$contains)
#   #       store[[cid]]$eps <- eps_values[n]
#   #     })
#   #   } else { ## Either a non-true-split multi-branch or all noise
#   #     
#   #   }
#   # }
#   # 
#   # for (i in 0:(nrow(mst)-1)){
#   #   j <- i + 1
#   #   new_clusters <- setdiff(res[, j], res[, i])
#   #   j_freq <- table(res[, j])
#   #   i_freq <- table(res[, i])
#   #   ## If the next tree cut down generated a new non-zero symbol, it was a true split
#   #   if (length(j_freq[names(j_freq) %in% "0"]) > length(i_freq[names(i_freq) %in% "0"])){
#   #     
#   #   } ## Else if noise points are falling off the cluster
#   #   else if (added_noise[i]) {
#   #     update_points <- which(res[, j] != res[, i])
#   #     sapply(update_points, function(cid_i){
#   #       cid <- as.character(res[cid_i, i])
#   #       store[[cid]]
#   #     })
#   #   } ## Else if only noise left
#   #   else {
#   #     
#   #   }
#   #   
#   #   ## Take advantage of the fact that it is known when the last true split happened
#   #   if (length(j_freq) < length(i_freq)){
#   #     only_noise_left <- T
#   #   }
#   #   
# #     ## Create entry in hash table
# #     # sapply(as.character(new_clusters), function(cid) assign(cid, list(), envir = store))
# #     
# #     ## Get the cluster ID that was replaced
# #     split_off_index <- res[, i][which((res[, j] - res[, i]) != 0)][1] #na.omit(unique(res[, i][which(res[, j] == new_clusters)]))
# #     
# #     ## Get
# #     new_branch <- j_freq[names(j_freq) %in% c(new_clusters, split_off_index)]
# #     if (all(new_branch > minPts)){
# #       if (j != 1){
# #         new_cids <- c(cid_count+1, cid_count+2)
# #         cid_count <- cid_count + 2
# #         to_replace <- res[, j] %in% c(new_clusters, split_off_index)
# #         res2[to_replace, j] <- (new_cids)[match(res[to_replace, j], sort(unique(res[to_replace, j])))]
# #         sapply(as.character(new_cids), function(cid) {
# #           store[[cid]]$contains <- unique(c(store[[cid]]$contains, which(res[, j] == as.numeric(cid))))
# #           store[[cid]]$eps <- rep(eps_values[j], length(store[[cid]]$contains))
# #         })
# #       } else {
# #         store[["1"]]$contains <- cl$hcl$order
# #         store[["1"]]$eps <- rep(eps_values[j], length(store[["1"]]$contains))
# #       }
# #     } else if (all(new_branch <= minPts)) {
# #       ## else just update noise
# #       store[[cid]]$eps[which(res[, j] != res[, i])] <- eps_values[j]
# #     } else {
# #       cid <- as.character(split_off_index)
# #       #store[[cid]]$contains <- unique(c(store[[cid]]$contains, which(res[, j] == as.numeric(cid))))
# #       store[[cid]]$eps[which(res[, j] != res[, i])] <- eps_values[j]
# #     }
# #   }
# # }
#   
#   
#   map_vals <- function(x, from, to){
#     
#   }
#   
#   
#   # mst <- cl$mst[order(-cl$mst[, 3]),]
#   # # Add microscopically small number for consistency between tree cuts and DBSCAN*
#   # cuts <- cutree(cl$hcl, h = mst[,3]+.Machine$double.eps) 
#   # core_dist <- diag(cl$mrd)
#   # 
#   # # res <- apply(cuts, 2, function(column){ column[core_dist > as.numeric()] <- 0; column })
#   # 
#   # apply(cuts, 2, function(branch) match(branch[hcl$order], unique(branch[hcl$order])))
#   # 
#   # 
#   # ## Fast cutree-based approach with a conditional using core distance 
#   # n <- nrow(mst)+1
#   # i_merge <- as.vector(hcl$merge)
#   # m_nr <- rep(0, n)
#   # j <- 0
#   # sing <- rep(T, n)
#   # 
#   # res <- matrix(0, nrow=nrow(mst)+1, ncol=nrow(mst))
#   # i <- 1
#   # for (dist in mst[1:20, 3]){
#   #   k <- n + 1L - apply(outer(c(hcl$height, Inf), dist, ">"), 2, which.max)
#   #   res[, i] <- my_cuttree(hcl$merge, k)
#   #   i <- i + 1
#   # }
# 
#   
#   # 
#   # # Magical k 
#   # k <- n + 1L - apply(outer(c(hcl$height, Inf), 0.25, ">"), 2, which.max)
#   # which <- k
#   # for (k in 1:n){
#   #   # m1 <- i_merge[k]
#   #   # m2 <- i_merge[n-1+k-1]
#   #   m <- hcl$merge[k,]
#   #   if (all(m < 0)){
#   #     m_nr[-m] <- k
#   #     sing[-m] <- F
#   #   } else if (any(m < 0)){
#   #     j <- ifelse(m[1] < 0, -m[1], m[2])
#   #     m_nr[m_nr == m[1]] <- k
#   #     m_nr[k] <- k
#   #     sing[j] <- F
#   #   } else {
#   #     m_nr[m_nr %in% m] <- k
#   #   }
#   # }
#   # 
#   # found_j <- F
#   # z <- rep(0, n)
#   # i_ans <- rep(0, n)
#   # 
#   # for (j in 1:length(which)){
#   #   if(which == n - k){
#   #     if (!found_j){
#   #       found_j <- T
#   #       nclust <- 0
#   #       mm <- j * n
#   #       m1 <- mm
#   #       for (l in 1:n){
#   #         if (sing[l]){
#   #           i_ans[m1] = nclust
#   #           nclust <- nclust + 1
#   #         }else {
#   #           if (z[m_nr[l]] == 0){
#   #             z[m_nr[l]] <- nclust 
#   #             nclust <- nclust + 1
#   #             i_ans[m1] <- z[m_nr[l]]
#   #           }
#   #         }
#   #         m1 <- m1 + 1
#   #       }
#   #     }
#   #   } else {
#   #     m1 <- j*n
#   #     m2 <- mm
#   #     for(l in 1:n){
#   #       i_ans[m1] = i_ans[m2]
#   #       m1 <- m1 + 1 
#   #       m2 <- m2 + 1
#   #     }
#   #   }
#   #   for(j in 1:length(which)){
#   #     if(which[j] == n){
#   #       m1 <- j*n
#   #       for(l in 1:n){
#   #         i_ans[m1] <- l
#   #         m1 <- m1 + 1
#   #       }
#   #     }
#   #   }
#   #   
#   # }
#   # 
#   
# 
#   
#   
#   
#   
#   # res <- sapply(1:nrow(mst), function(i) {
#   #   #cuts[core_dist > mst[i,3], i] <- 0
#   #   #cuts[,i]
#   # noise <- which(core_dist > mst[i,3])
#   # if (length(noise) > 0){
#   #   cuts[noise, i] <- 0
#   #   corder <- unique(cuts[, i][cl$hcl$order])
#   #   corder <- corder[corder != 0]
#   #   return(match(cuts[, i], c(0, unique(corder))) - 1)
#   # }
#   # #return(cuts[,i])
#   # return(match(cuts[, i], unique(cuts[, i][cl$hcl$order])))
#   # })
#   
#   # for (i in 1:ncol(cuts)){
#   #   labels <- cuts[, i]
#   #   noise <- which(core_dist > mst[i, 3])
#   #   if (length(noise) > 0) 
#   #     labels[core_dist > mst[i, 3]] <- 0 
#   #   }
#   # }
#   # 
#   # res2 <- res
#   # minPts <- 5
#   # store <- new.env(parent = emptyenv())
#   # eps_values <- mst[, 3]
#   # cid_count <- 1
#   # only_noise_left <- F
#   # 
#   # added_noise <- mapply(function(i, j){
#   #   ifelse(is.na(table(res[, j])["0"]), -Inf, table(res[, j])["0"]) > ifelse(is.na(table(res[, i])["0"]), 0, table(res[, i])["0"])
#   # }, 1:(ncol(res)-1), 2:ncol(res))
#   # 
#   # ## Get the epsilon indices that produced 'true' cluster splits
#   # true_splits <- which(mapply(function(i, j) {
#   #   length(table(res[res[, i] != 0, i])) < length(table(res[res[, j] != 0, j]))
#   #   }, 1:(ncol(res)-1), 2:ncol(res))) + 1
#   # length_changes <- sapply(apply(res, 2, table), length)
#   # true_splits <- true_splits[true_splits < max(which(length_changes == max(length_changes)))]
#   # 
#   # for (n in 2:nrow(mst)){
#   #   assignments <- res[, n][cl$hcl$order]
#   #   cids <- unique(assignments)
#   # 
#   #   if (n %in% true_splits) { ## true split
#   #   sapply(as.character(cids), function(cid) {
#   #     store[[cid]]$contains <- cl$hcl$order[res[, n][cl$hcl$order] == as.numeric(cid)]
#   #     store[[cid]]$eps <- eps_values[n]
#   #   })
#   #   } else if (added_noise[n]){ ## cluster is shrinking
#   #     sapply(as.character(cids), function(cid) {
#   #       index <- match(cl$hcl$order[res[, n][cl$hcl$order] == as.numeric(cid)], store[[cid]]$contains)
#   #       store[[cid]]$eps <- eps_values[n]
#   #     })
#   #   } else { ## Either a non-true-split multi-branch or all noise
#   # 
#   #   }
#   # }
#   # 
#   # for (i in 0:(nrow(mst)-1)){
#   #   j <- i + 1
#   #   new_clusters <- setdiff(res[, j], res[, i])
#   #   j_freq <- table(res[, j])
#   #   i_freq <- table(res[, i])
#   #   ## If the next tree cut down generated a new non-zero symbol, it was a true split
#   #   if (length(j_freq[names(j_freq) %in% "0"]) > length(i_freq[names(i_freq) %in% "0"])){
#   # 
#   #   } ## Else if noise points are falling off the cluster
#   #   else if (added_noise[i]) {
#   #     update_points <- which(res[, j] != res[, i])
#   #     sapply(update_points, function(cid_i){
#   #       cid <- as.character(res[cid_i, i])
#   #       store[[cid]]
#   #     })
#   #   } ## Else if only noise left
#   #   else {
#   # 
#   #   }
#   # 
#   #   ## Take advantage of the fact that it is known when the last true split happened
#   #   if (length(j_freq) < length(i_freq)){
#   #     only_noise_left <- T
#   #   }
#   # 
#   #     ## Create entry in hash table
#   #     # sapply(as.character(new_clusters), function(cid) assign(cid, list(), envir = store))
#   # 
#   #     ## Get the cluster ID that was replaced
#   #     split_off_index <- res[, i][which((res[, j] - res[, i]) != 0)][1] #na.omit(unique(res[, i][which(res[, j] == new_clusters)]))
#   # 
#   # ## Get
#   # new_branch <- j_freq[names(j_freq) %in% c(new_clusters, split_off_index)]
#   # if (all(new_branch > minPts)){
#   #   if (j != 1){
#   #     new_cids <- c(cid_count+1, cid_count+2)
#   #     cid_count <- cid_count + 2
#   #     to_replace <- res[, j] %in% c(new_clusters, split_off_index)
#   #     res2[to_replace, j] <- (new_cids)[match(res[to_replace, j], sort(unique(res[to_replace, j])))]
#   #     sapply(as.character(new_cids), function(cid) {
#   #       store[[cid]]$contains <- unique(c(store[[cid]]$contains, which(res[, j] == as.numeric(cid))))
#   #       store[[cid]]$eps <- rep(eps_values[j], length(store[[cid]]$contains))
#   #     })
#   #   } else {
#   #     store[["1"]]$contains <- cl$hcl$order
#   #     store[["1"]]$eps <- rep(eps_values[j], length(store[["1"]]$contains))
#   #   }
#   # } else if (all(new_branch <= minPts)) {
#   #   ## else just update noise
#   #   store[[cid]]$eps[which(res[, j] != res[, i])] <- eps_values[j]
#   # } else {
#   #   cid <- as.character(split_off_index)
#   #   #store[[cid]]$contains <- unique(c(store[[cid]]$contains, which(res[, j] == as.numeric(cid))))
#   #   store[[cid]]$eps[which(res[, j] != res[, i])] <- eps_values[j]
#   # }
#   #   }
#   # }
#   #cuts[which(core_dist > eps)] <- 0 # Use core distance to distinguish noise; minPts doesn't matter? 
#   #cut_cl <- match(cuts, sort(unique(cuts))) - 1 # 'Normalize' so cluster ID increments; sort ensures noise is 0
# }


