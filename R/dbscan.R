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

dbscan <- function(x, eps, minPts = 5, weights = NULL,
  borderPoints = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0) {

  ## make sure x is numeric
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")


  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR"))
  if(is.na(search)) stop("Unknown NN search type!")

  ret <- dbscan_int(x, as.double(eps), as.integer(minPts),
    as.double(weights), as.integer(borderPoints),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx))

  ret <- list(cluster = ret, eps = eps, minPts = minPts)
  class(ret) <- "dbscan"
  ret
}


print.dbscan <- function(x, ...) {
  cat("DBSCAN clustering for ", length(x$cluster), " objects.", "\n", sep = "")
  cat("Parameters: eps = ", x$eps, ", minPts = ", x$minPts, "\n", sep = "")
  cl <- unique(x$cluster)
  cl <- length(cl[cl!=0L])
  cat("The clustering contains ", cl, " cluster(s).",
      "\n", sep = "")
  cat("Available fields: ", paste(names(x), collapse = ", "), "\n", sep = "")
}

