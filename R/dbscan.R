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

dbscan <- function(x, eps, minPts = 5, weights = NULL, borderPoints = TRUE,
  ...) {

  if(is(x, "frNN") && missing(eps)) eps <- x$eps

  ### extra contains settings for frNN
  ### search = "kdtree", bucketSize = 10, splitRule = "suggest", approx = 0
  ### also check for MinPts for fpc compartibility (does not work for
  ### search method dist)
  extra <- list(...)
  args <- c("MinPts", "search", "bucketSize", "splitRule", "approx")
  m <- pmatch(names(extra), args)
  if(any(is.na(m))) stop("Unknown parameter: ",
    paste(names(extra)[is.na(m)], collapse = ", "))
  names(extra) <- args[m]

  if(!is.null(extra$MinPts)) {
    warning("converting argument MinPts (fpc) to minPts (dbscan)!")
    minPts <- extra$MinPts
  }
  extra$MinPts <- NULL

  search <- if(is.null(extra$search)) "kdtree" else extra$search
  search <- pmatch(toupper(search), c("KDTREE", "LINEAR", "DIST"))
  if(is.na(search)) stop("Unknown NN search type!")

  bucketSize <- if(is.null(extra$bucketSize)) 10L else
    as.integer(extra$bucketSize)

  splitRule <- if(is.null(extra$splitRule)) "suggest" else extra$splitRule
  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  approx <- if(is.null(extra$approx)) 0L else as.integer(extra$approx)

  ### do dist search
  if(search == 3) {
    if(!is(x, "dist"))
      if(.matrixlike(x)) x <- dist(x)
      else stop("x needs to be a matrix to calculate distances")
  }

  ## for dist we provide the R code with a frNN list and no x
  frNN <- list()
  if(is(x, "dist")) {
    frNN <- frNN(x, eps, ...)$id
    x <- matrix(0.0, nrow=0, ncol=0)
  }else if(is(x, "frNN")) {
    if(x$eps != eps) {
      eps <- x$eps
      warning("Using the eps of ", eps, " provided in the fixed-radius NN object.")
    }
    frNN <- x$id
    x <- matrix(0.0, nrow=0, ncol=0)

  }else{
    if(!.matrixlike(x)) stop("x needs to be a matrix")
    ## make sure x is numeric
    x <- as.matrix(x)
    if(storage.mode(x) == "integer") storage.mode(x) <- "double"
    if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")
  }

  if(length(frNN) == 0 && any(is.na(x)))
    stop("data/distances cannot contain NAs for dbscan (with kd-tree)!")

  ## add self match and use C numbering if frNN is used
  if(length(frNN) > 0)
    frNN <- lapply(1:length(frNN), FUN = function(i) c(i-1L, frNN[[i]]-1L))

  ret <- dbscan_int(x, as.double(eps), as.integer(minPts),
    as.double(weights), as.integer(borderPoints),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx), frNN)

  structure(list(cluster = ret, eps = eps, minPts = minPts),
    class = c("dbscan_fast", "dbscan"))
}


print.dbscan_fast <- function(x, ...) {
  cat("DBSCAN clustering for ", length(x$cluster), " objects.", "\n", sep = "")
  cat("Parameters: eps = ", x$eps, ", minPts = ", x$minPts, "\n", sep = "")
  cl <- unique(x$cluster)
  cl <- length(cl[cl!=0L])
  cat("The clustering contains ", cl, " cluster(s) and ", sum(x$cluster==0L),
    " noise points.",
      "\n", sep = "")
  print(table(x$cluster))
  message <- paste("\nAvailable fields: ", paste(names(x), collapse = ", "), "\n", sep="")
  cat(paste0(strwrap(message, width = 75), collapse = "\n                  "))
}

predict.dbscan_fast <- function (object, newdata = NULL, data, ...) {
  if (is.null(newdata)) return(object$cluster)

  nn <- frNN(rbind(data, newdata), eps = object$eps,
    sort = TRUE, ...)$id[-(1:nrow(data))]
  sapply(nn, function(x) {
    x <- x[x<=nrow(data)]
    x <- object$cluster[x][x>0][1]
    x[is.na(x)] <- 0L
    x
    })
}
