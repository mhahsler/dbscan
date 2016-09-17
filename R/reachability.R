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


# Add new Generic for other class extensions support reachability plotting 
as.reachability <- function(x, ...) UseMethod("as.reachability")

# Dendrogram --> Reachability conversion 
as.dendrogram.reachability <- function(x) {
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
  
  # Always shift down one; plot(*, center=T) needed since midpoint calculation expensive
  return(leaves[[1]])
}

# Converting from dendrogram --> reachability plot 
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
  return(structure(list(reachdist=rp[["reachdist"]], order=rp[["order"]])), class="reachability")
}

# Plotting method for reachability objects.
plot.reachability <- function(x){
  
}
