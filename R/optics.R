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

optics <- function(x, eps, minPts = 5, eps_cl, search = "kdtree",
  bucketSize = 10, splitRule = "suggest", approx = 0) {

  ## make sure x is numeric
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")

  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR"))
  if(is.na(search)) stop("Unknown NN search type!")

  ret <- optics_int(as.matrix(x), as.double(eps), as.integer(minPts),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx))

  ### find clusters
  if(!missing(eps_cl)) ret <- optics_cut(ret, eps_cl)

  ret$eps <- eps
  ret$minPts <- minPts

  class(ret) <- "optics"
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
  cat("Parameters: eps = ", x$eps, ", minPts = ", x$minPts, "\n", sep = "")
  if(!is.null(x$cluster)) {
    cl <- unique(x$cluster)
    cl <- length(cl[cl!=0L])
    cat("The clustering contains ", cl, " cluster(s).",
      "\n", sep = "")
  }
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
  }

plot.optics <- function(x, y=NULL, cluster = TRUE, ...) {
    if(!is.null(x$cluster) && cluster) {
      plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
        ylab = "Reachability dist.", xlab = "OPTICS order",
        main = "Reachability Plot")
      # abline(h=x$eps_cl, col="gray", lty=2)
    }else{
      plot(x$reachdist[x$order], type="h",
        ylab = "Reachability dist.", xlab = "OPTICS order",
        main = "Reachability Plot")
    }
}
