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

opticsXi <- function(optics_obj, xi = 0.001, minimum=F)
{
  if (!"optics" %in% class(optics_obj)) stop("opticsXi only accepts objects resulting from dbscan::optics!")
  
  # Initial variables  
  optics_obj$ord_rd <- optics_obj$reachdist[optics_obj$order]
  optics_obj$ixi <- (1 - xi)
  SetOfSteepDownAreas <- list()
  SetOfClusters <- list()
  index <- 1
  mib <- 0
  sdaset <- list()
  nocorrect <- F
  while (index <= length(optics_obj$order))
  {
    mib <- max(mib, optics_obj$ord_rd[index]) 
    if (!valid(index+1, optics_obj)) break
    
    # Test if this is a steep down area 
    if (steepDown(index, optics_obj))
    {
      # Update mib values with current mib and filter 
      sdaset <- updateFilterSDASet(mib, sdaset, optics_obj$ixi)
      startval <- optics_obj$ord_rd[index]
      mib <- 0
      startsteep <- index; endsteep <- index + 1
      while(!is.na(optics_obj$order[index+1])) {
        index <- index + 1
        if (steepDown(index, optics_obj)) { endsteep <- index + 1; next }
        if (!steepDown(index, optics_obj, ixi=1.0) || index - endsteep > optics_obj$minPts) break
      }
      sda <- list(s=startsteep, e=endsteep, maximum=startval, mib=0)
      # print(paste("New steep down area:", toString(sda)))
      sdaset <- append(sdaset, list(sda))
      next
    }
    if (steepUp(index, optics_obj))
    {
      sdaset <- updateFilterSDASet(mib, sdaset, optics_obj$ixi)
      {
        startsteep <- index; endsteep <- index + 1
        mib <- optics_obj$ord_rd[index]
        esuccr <- if (!valid(index+1, optics_obj)) Inf else optics_obj$ord_rd[index+1]
        if (esuccr != Inf) {
          while(!is.na(optics_obj$order[index+1])) {
            index <- index + 1
            if (steepUp(index, optics_obj)) { 
              endsteep <- index + 1
              mib <- optics_obj$ord_rd[index]
              esuccr <- if (!valid(index+1, optics_obj)) Inf else optics_obj$ord_rd[index+1]
              if (esuccr == Inf) { endsteep <- endsteep - 1; break }
              next 
            }
            if (!steepUp(index, optics_obj, ixi=1.0) || index - endsteep > optics_obj$minPts) break
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
        if (mib * optics_obj$ixi < sda$mib) next 
        
        # Default values 
        cstart <- sda$s
        cend <- sua$e
        
        # Credit to ELKI 
        if (!nocorrect) while(cend > cstart && optics_obj$ord_rd[cend] == Inf) { cend <- cend - 1 }
        
        # Condition 4
        {
          # Case b
          if (sda$maximum * optics_obj$ixi >= sua$maximum) {
            while(cstart < cend && optics_obj$ord_rd[cstart+1] > sua$maximum) cstart <- cstart + 1
          } 
          # Case c
          else if (sua$maximum * optics_obj$ixi >= sda$maximum) {
            while(cend > cstart && optics_obj$ord_rd[cend-1] > sda$maximum) cend <- cend - 1
          }
        }
        
        # Credit to ELKI 
        if (!nocorrect) {
          simplify <- F
          while (cend > cstart) {
            tmp2 <- optics_obj$order[cend]
            for (i in start:cend) {
              if (tmp2 == optics_obj$order[i]) simplify <- T
            }
            if (simplify) break
            cend <- cend - 1
          }
        }
        
        # obey minpts 
        if (cend - cstart + 1 < optics_obj$minPts) next
        SetOfClusters <- append(SetOfClusters, list(list(start=cstart, end=cend)))
        next
      }
    } else { index <- index + 1 }
  }
  # Remove aliases  
  optics_obj$ord_rd <- NULL
  optics_obj$ixi <- NULL

  # Keep hiearchical cluster data
  optics_obj$xi <- xi
  optics_obj$clusters_xi <- data.table::rbindlist(SetOfClusters)

  # Replace cluster attribute with XI results and return
  optics_obj <- extractXiClusters(optics_obj, minimum=minimum)
  optics_obj
}

# Removes obsolete steep areas 
updateFilterSDASet <- function(mib, sdaset, ixi) {
  sdaset <- Filter(function(sda) sda$maximum*ixi > mib, sdaset)
  lapply(sdaset, function(sda) { if (mib > sda$mib) sda$mib <- mib; sda })
}

steepUp <- function(i, optics_obj, ixi = optics_obj$ixi) {
  if(optics_obj$ord_rd[i] >= Inf) return(F)
  if(!valid(i+1, optics_obj)) return(T)
  return(optics_obj$ord_rd[i] <= optics_obj$ord_rd[i+1] * ixi)
}

steepDown <- function(i, optics_obj, ixi = optics_obj$ixi) {
  if(!valid(i+1, optics_obj)) return(F)
  if(optics_obj$ord_rd[i+1] >= Inf) return(F)
  return(optics_obj$ord_rd[i] * ixi >= optics_obj$ord_rd[i+1])
}

valid <- function(index, optics_obj) {
  return(!is.na(optics_obj$ord_rd[index]))
}

### Extract xi clusters (minimum == T extracts clusters that do not contain other clusters)
extractXiClusters <- function(optics_obj, minimum=F) {
  
  if (length(optics_obj$clusters_xi) == 0) stop("No Xi clusters detected. The steepness (xi) parameter might be too high.")
  
  # Prepare cluster list based on minimum parameter
  clusters_xi <- optics_obj$clusters_xi
  if (!"cluster_id" %in% names(clusters_xi)) clusters_xi <- cbind(clusters_xi, cluster_id=1:nrow(clusters_xi))
  if (!"cluster_size" %in% names(clusters_xi)) clusters_xi <- cbind(clusters_xi, clusters_xi[, list(cluster_size = (end - start))])
  clusters_xi <- if (minimum) { clusters_xi[order(cluster_size)] } else { clusters_xi[order(-cluster_size)] }
  
  # Fill in the matrix 
  clusters <- rep(0, length(optics_obj$order))
  for(cid in clusters_xi$cluster_id) {
    if (minimum) {
      if (clusters_xi[cluster_id == cid, all(clusters[start:end] == 0)]) { 
        clusters[clusters_xi[cluster_id == cid, start:end]] <- cid
      }
    } else clusters[clusters_xi[cluster_id == cid, start:end]] <- cid
  }
  
  # Fix the ordering
  clusters[optics_obj$order] <- clusters
  optics_obj$cluster <- clusters
  optics_obj
}