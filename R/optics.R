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
    if(is.null(x$xi)) { 
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

plot.optics <- function(x, y=NULL, cluster = TRUE, ...) {
    if(!is.null(x$cluster) && cluster) {
      if(is.null(x$clusters_xi)) { 
        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
          ylab = "Reachability dist.", xlab = "OPTICS order",
          main = "Reachability Plot", ...)
      } else {
        y_max <- max(x$reachdist[which(x$reachdist != Inf)])
        hclusters <- x$clusters_xi[order(x$clusters_xi$end-x$clusters_xi$start),]
        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
             ylab = "Reachability dist.", xlab = NA, xaxt = "n",
             main = "Reachability Plot", yaxs="i", ylim=c(0,y_max), xaxt='n', ...)
        y_increments <- ((y_max/(par("plt")[4]-par("plt")[3]))*(par("plt")[3]))/(2*nrow(hclusters))
        i <- 1:nrow(hclusters)
        segments(x0=hclusters$start[i], y0=-(y_increments*i), x1=hclusters$end[i], col=hclusters$cluster_id[i]+1L, lwd=1, xpd=T)
      }
    }else{
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
