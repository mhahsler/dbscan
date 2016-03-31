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

opticsXi <- function(object, xi = 0.001, minimum=F)
{
  if (!"optics" %in% class(object)) stop("opticsXi only accepts objects resulting from dbscan::optics!")
  if (xi >= 1.0 || xi <= 0.0) stop("The Xi parameter must be (0, 1)")
  
  # Initial variables  
  object$ord_rd <- object$reachdist[object$order]
  object$ixi <- (1 - xi)
  SetOfSteepDownAreas <- list()
  SetOfClusters <- list()
  index <- 1
  mib <- 0
  sdaset <- list()
  while (index <= length(object$order))
  {
    mib <- max(mib, object$ord_rd[index]) 
    if (!valid(index+1, object)) break
    
    # Test if this is a steep down area 
    if (steepDown(index, object))
    {
      # Update mib values with current mib and filter 
      sdaset <- updateFilterSDASet(mib, sdaset, object$ixi)
      startval <- object$ord_rd[index]
      mib <- 0
      startsteep <- index; endsteep <- index + 1
      while(!is.na(object$order[index+1])) {
        index <- index + 1
        if (steepDown(index, object)) { endsteep <- index + 1; next }
        if (!steepDown(index, object, ixi=1.0) || index - endsteep > object$minPts) break
      }
      sda <- list(s=startsteep, e=endsteep, maximum=startval, mib=0)
      # print(paste("New steep down area:", toString(sda)))
      sdaset <- append(sdaset, list(sda))
      next
    }
    if (steepUp(index, object))
    {
      sdaset <- updateFilterSDASet(mib, sdaset, object$ixi)
      {
        startsteep <- index; endsteep <- index + 1
        mib <- object$ord_rd[index]
        esuccr <- if (!valid(index+1, object)) Inf else object$ord_rd[index+1]
        if (esuccr != Inf) {
          while(!is.na(object$order[index+1])) {
            index <- index + 1
            if (steepUp(index, object)) { 
              endsteep <- index + 1
              mib <- object$ord_rd[index]
              esuccr <- if (!valid(index+1, object)) Inf else object$ord_rd[index+1]
              if (esuccr == Inf) { endsteep <- endsteep - 1; break }
              next 
            }
            if (!steepUp(index, object, ixi=1.0) || index - endsteep > object$minPts) break
          }
        } else {
          endsteep <- endsteep - 1
          index <- index + 1
        }
        sua <- list(s=startsteep, e=endsteep, maximum=esuccr)
        # print(paste("New steep up area:", toString(sua)))
      }
      for (sda in rev(sdaset))
      {
        # Condition 3B 
        if (mib * object$ixi < sda$mib) next 
        
        # Default values 
        cstart <- sda$s
        cend <- sua$e
      
        # Condition 4
        {
          # Case b
          if (sda$maximum * object$ixi >= sua$maximum) {
            while(cstart < cend && object$ord_rd[cstart+1] > sua$maximum) cstart <- cstart + 1
          } 
          # Case c
          else if (sua$maximum * object$ixi >= sda$maximum) {
            while(cend > cstart && object$ord_rd[cend-1] > sda$maximum) cend <- cend - 1
          }
        }
        
        # obey minpts 
        if (cend - cstart + 1 < object$minPts) next
        SetOfClusters <- append(SetOfClusters, list(list(start=cstart, end=cend)))
        next
      }
    } else { index <- index + 1 }
  }
  # Remove aliases  
  object$ord_rd <- NULL
  object$ixi <- NULL

  # Keep hiearchical cluster data
  object$xi <- xi
  object$clusters_xi <- do.call(rbind, SetOfClusters)
  
  if (length(SetOfClusters) == 0) stop(paste("No clusters were found with threshold:", xi))
  else { 
    # Order cluster data 
    object$clusters_xi <- data.frame(start=unlist(object$clusters_xi[,1]), end=unlist(object$clusters_xi[,2]))
    object$clusters_xi <- object$clusters_xi[order(object$clusters_xi$start, object$clusters_xi$end),]
    row.names(object$clusters_xi) <- NULL
  }

  # Replace cluster attribute with XI results and return
  object <- extractXiClusters(object, minimum=minimum)
  object
}

# Removes obsolete steep areas 
updateFilterSDASet <- function(mib, sdaset, ixi) {
  sdaset <- Filter(function(sda) sda$maximum*ixi > mib, sdaset)
  lapply(sdaset, function(sda) { if (mib > sda$mib) sda$mib <- mib; sda })
}

steepUp <- function(i, object, ixi = object$ixi) {
  if(object$ord_rd[i] >= Inf) return(F)
  if(!valid(i+1, object)) return(T)
  return(object$ord_rd[i] <= object$ord_rd[i+1] * ixi)
}

steepDown <- function(i, object, ixi = object$ixi) {
  if(!valid(i+1, object)) return(F)
  if(object$ord_rd[i+1] >= Inf) return(F)
  return(object$ord_rd[i] * ixi >= object$ord_rd[i+1])
}

valid <- function(index, object) {
  return(!is.na(object$ord_rd[index]))
}

### Extract xi clusters (minimum == T extracts clusters that do not contain other clusters)
extractXiClusters <- function(object, minimum=F) {
  if (length(object$clusters_xi) == 0) stop("No Xi clusters detected. The steepness (xi) parameter might be too high.")
  
  # Add cluster_id to xi clusters 
  if (!"cluster_id" %in% names(object$clusters_xi)) object$clusters_xi <- cbind(object$clusters_xi, cluster_id=1:nrow(object$clusters_xi))
 
  # Copy the clusters and sort them based on minimum parameter value 
  clusters_xi <- object$clusters_xi
  if (!"cluster_size" %in% names(clusters_xi)) clusters_xi <- cbind(clusters_xi, list(cluster_size = (clusters_xi$end - clusters_xi$start)))
  clusters_xi <- if (minimum) { clusters_xi[order(clusters_xi$cluster_size),] } else { clusters_xi[order(-clusters_xi$cluster_size),] }
  
  # Fill in the matrix 
  clusters <- rep(0, length(object$order))
  for(cid in clusters_xi$cluster_id) {
    cluster <- clusters_xi[clusters_xi$cluster_id == cid,]
    if (minimum) {
      if (all(clusters[cluster$start:cluster$end] == 0)) { 
        clusters[cluster$start:cluster$end] <- cid
      }
    } else clusters[cluster$start:cluster$end] <- cid
  }
  
  # Fix the ordering
  clusters[object$order] <- clusters
  object$cluster <- clusters
  object
}