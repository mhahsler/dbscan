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

### extract DBSCAN clustering
extractDBSCAN <- function(object, eps_cl) {
  if (!"optics" %in% class(object)) stop("extractDBSCAN only accepts objects resulting from dbscan::optics!")

  reachdist <- object$reachdist[object$order]
  coredist <- object$coredist[object$order]
  n <- length(object$order)
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

  object$eps_cl <- eps_cl
  object$xi <- NA
  ### fix the order so cluster is in the same order as the rows in x
  cluster[object$order] <- cluster
  object$cluster <- cluster

  object
}
