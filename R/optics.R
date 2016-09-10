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

optics <- function(x, eps, minPts = 5, ...) {

  ### extra contains settings for frNN
  ### search = "kdtree", bucketSize = 10, splitRule = "suggest", approx = 0
  extra <- list(...)
  args <- c("search", "bucketSize", "splitRule", "approx")
  m <- pmatch(names(extra), args)
  if(any(is.na(m))) stop("Unknown parameter: ",
    paste(names(extra)[is.na(m)], collapse = ", "))
  names(extra) <- args[m]

  search <- if(is.null(extra$search)) "kdtree" else extra$search
  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if(is.na(search)) stop("Unknown NN search type!")

  bucketSize <- if(is.null(extra$bucketSize)) 10L else
    as.integer(extra$bucketSize)

  splitRule <- if(is.null(extra$splitRule)) "suggest" else extra$splitRule
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  approx <- if(is.null(extra$approx)) 0L else as.integer(extra$approx)

  ### dist search
  if(search == 3) {
    if(!is(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  ## for dist we provide the R code with a frNN list and no x
  frNN <- list()
  if(is(x, "dist")) {
    frNN <- frNN(x, eps, ...)
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

  ret
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

plot.optics <- function(x, y=NULL, cluster = TRUE, ...) {
    if(!is.null(x$cluster) && cluster) {
      if(is.null(x$clusters_xi)) {
        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
          ylab = "Reachability dist.", xlab = "OPTICS order",
          main = "Reachability Plot", ...)
      } else {
        y_max <- max(x$reachdist[which(x$reachdist != Inf)])
        # Sort clusters by size
        hclusters <- x$clusters_xi[order(x$clusters_xi$end-x$clusters_xi$start),]

        def.par <- par(no.readonly = TRUE)

        # .1 means to leave 15% for the cluster lines
        par(mar= c(2, 4, 4, 2) + 0.1, omd = c(0, 1, .15, 1))
        y_increments <- (y_max/0.85*.15)/(nrow(hclusters)+1L)

        plot(x$reachdist[x$order], type="h", col=x$cluster[x$order]+1L,
             ylab = "Reachability dist.", xlab = NA, xaxt = "n",
             main = "Reachability Plot", yaxs="i", ylim=c(0,y_max), xaxt='n', ...)

        i <- 1:nrow(hclusters)
        segments(x0=hclusters$start[i], y0=-(y_increments*i),
          x1=hclusters$end[i], col=hclusters$cluster_id[i]+1L, lwd=2, xpd=NA)
        par(def.par)

        }
    }else{
      plot(x$reachdist[x$order], type="h",
        ylab = "Reachability dist.", xlab = "OPTICS order",
        main = "Reachability Plot", ...)
    }
}

predict.optics <- function (object, newdata = NULL, data, ...) {
  if (is.null(newdata)) return(object$cluster)
  if (is.null(object$cluster)) stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")

  nn <- frNN(rbind(data, newdata), eps = object$eps_cl,
    sort = TRUE, ...)$id[-(1:nrow(data))]
  sapply(nn, function(x) {
    x <- x[x<=nrow(data)]
    x <- object$cluster[x][x>0][1]
    x[is.na(x)] <- 0L
    x
  })
}
