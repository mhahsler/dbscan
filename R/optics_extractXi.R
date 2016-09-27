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

extractXi <- function(x, xi, minimum = FALSE, correctPredecessors = TRUE)
{
  if (!"optics" %in% class(x)) stop("extractXi only accepts xs resulting from dbscan::optics!")
  if (xi >= 1.0 || xi <= 0.0) stop("The Xi parameter must be (0, 1)")

  # Initial variables
  x$ord_rd <- x$reachdist[x$order]
  x$ixi <- (1 - xi)
  SetOfSteepDownAreas <- list()
  SetOfClusters <- list()
  index <- 1
  mib <- 0
  sdaset <- list()
  while (index <= length(x$order))
  {
    mib <- max(mib, x$ord_rd[index])
    if (!valid(index+1, x)) break

    # Test if this is a steep down area
    if (steepDown(index, x))
    {
      # Update mib values with current mib and filter
      sdaset <- updateFilterSDASet(mib, sdaset, x$ixi)
      startval <- x$ord_rd[index]
      mib <- 0
      startsteep <- index; endsteep <- index + 1
      while(!is.na(x$order[index+1])) {
        index <- index + 1
        if (steepDown(index, x)) { endsteep <- index + 1; next }
        if (!steepDown(index, x, ixi=1.0) || index - endsteep > x$minPts) break
      }
      sda <- list(s=startsteep, e=endsteep, maximum=startval, mib=0)
      # print(paste("New steep down area:", toString(sda)))
      sdaset <- append(sdaset, list(sda))
      next
    }
    if (steepUp(index, x))
    {
      sdaset <- updateFilterSDASet(mib, sdaset, x$ixi)
      {
        startsteep <- index; endsteep <- index + 1
        mib <- x$ord_rd[index]
        esuccr <- if (!valid(index+1, x)) Inf else x$ord_rd[index+1]
        if (esuccr != Inf) {
          while(!is.na(x$order[index+1])) {
            index <- index + 1
            if (steepUp(index, x)) {
              endsteep <- index + 1
              mib <- x$ord_rd[index]
              esuccr <- if (!valid(index+1, x)) Inf else x$ord_rd[index+1]
              if (esuccr == Inf) { endsteep <- endsteep - 1; break }
              next
            }
            if (!steepUp(index, x, ixi=1.0) || index - endsteep > x$minPts) break
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
        if (mib * x$ixi < sda$mib) next

        # Default values
        cstart <- sda$s
        cend <- sua$e

        # Credit to ELKI
        if(!correctPredecessors) {
          while(cend > cstart && x$ord_rd[cend] == Inf) {
            cend <- cend - 1
          }
        }

        # Condition 4
        {
          # Case b
          if (sda$maximum * x$ixi >= sua$maximum) {
            while(cstart < cend && x$ord_rd[cstart+1] > sua$maximum) cstart <- cstart + 1
          }
          # Case c
          else if (sua$maximum * x$ixi >= sda$maximum) {
            while(cend > cstart && x$ord_rd[cend-1] > sda$maximum) cend <- cend - 1
          }
        }

        # This NOT in the original article - credit to ELKI for finding this.
        # Ensure that the predecessor is in the current cluster. This filter
        # removes common artifacts from the Xi method
        if(!correctPredecessors) {
          while(cend > cstart) {
            tmp2 <- x$predecessor[x$order[cend]]
            if (!is.na(tmp2) && any(x$order[cstart:(cend-1)] == tmp2, na.rm = TRUE))
               break
            # Not found.
            cend <- cend - 1
          }
        }

        # Ensure the last steep up point is not included if it's xi significant
        if (steepUp(index-1, x)) {
          cend <- cend - 1
        }

        # obey minpts
        if (cend - cstart + 1 < x$minPts) next
        SetOfClusters <- append(SetOfClusters, list(list(start=cstart, end=cend)))
        next
      }
    } else { index <- index + 1 }
  }
  # Remove aliases
  x$ord_rd <- NULL
  x$ixi <- NULL
  
  # Keep xi parameter
  x$xi <- xi
  
  # Zero-out clusters (only noise) if none found
  if (length(SetOfClusters) == 0) {
    warning(paste("No clusters were found with threshold:", xi))
    x$clusters_xi <- NULL
    x$cluster < rep(0, length(x$cluster))
    return(invisible(x))
  } else { # Cluster data exists; organize it by starting and ending index
    x$clusters_xi <- do.call(rbind, SetOfClusters)
    x$clusters_xi <- data.frame(start=unlist(x$clusters_xi[,1]), end=unlist(x$clusters_xi[,2]))
    x$clusters_xi <- x$clusters_xi[order(x$clusters_xi$start, x$clusters_xi$end),]
    row.names(x$clusters_xi) <- NULL
  }

  # Populate cluster vector with top-level cluster labels 
  x$cluster <- extractClusterLabels(x, minimum=minimum)
  
  # Remove non-local clusters if minimum was specified
  if (minimum) { x$clusters_xi <- x$clusters_xi[unique(x$cluster)[-1],] }
  
  class(x$cluster) <- append(class(x$cluster), "xics")
  class(x$clusters_xi) <- append(class(x$clusters_xi), "xics")
  x
}

# Removes obsolete steep areas
updateFilterSDASet <- function(mib, sdaset, ixi) {
  sdaset <- Filter(function(sda) sda$maximum*ixi > mib, sdaset)
  lapply(sdaset, function(sda) { if (mib > sda$mib) sda$mib <- mib; sda })
}

# Determines if the reachability distance at the current index 'i' is
# (xi) significantly lower than the next index
steepUp <- function(i, x, ixi = x$ixi) {
  if(x$ord_rd[i] >= Inf) return(FALSE)
  if(!valid(i+1, x)) return(TRUE)
  return(x$ord_rd[i] <= x$ord_rd[i+1] * ixi)
}

# Determines if the reachability distance at the current index 'i' is
# (xi) significantly higher than the next index
steepDown <- function(i, x, ixi = x$ixi) {
  if(!valid(i+1, x)) return(FALSE)
  if(x$ord_rd[i+1] >= Inf) return(FALSE)
  return(x$ord_rd[i] * ixi >= x$ord_rd[i+1])
}

# Determines if the reachability distance at the current index 'i' is a valid distance
valid <- function(index, x) {
  return(!is.na(x$ord_rd[index]))
}

### Extract xi clusters (minimum == T extracts clusters that do not contain other clusters)
extractClusterLabels <- function(cl, order, minimum = FALSE) {
  ## Add cluster_id to xi clusters
  if (!"cluster_id" %in% names(cl)) cl <- cbind(cl, cluster_id=1:nrow(cl))
  if (!all(c("start", "end") %in% names(cl))) stop("extractClusterLabels expects start and end references")
  
  ## Sort cl based on minimum parameter / cluster size
  if (!"cluster_size" %in% names(cl)) cl <- cbind(cl, list(cluster_size = (cl$end - cl$start)))
  cl <- if (minimum) { cl[order(cl$cluster_size),] } else { cl[order(-cl$cluster_size),] }

  ## Fill in the [cluster] vector with cluster IDs
  clusters <- rep(0, length(order))
  for(cid in cl$cluster_id) {
    cluster <- cl[cl$cluster_id == cid,]
    if (minimum) {
      if (all(clusters[cluster$start:cluster$end] == 0)) {
        clusters[cluster$start:cluster$end] <- cid
      }
    } else clusters[cluster$start:cluster$end] <- cid
  }

  # Fix the ordering
  clusters[order] <- clusters
  return(clusters)
}
