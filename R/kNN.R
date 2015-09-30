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


.ANNsplitRule <- c("STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR", "SUGGEST")

kNN <- function(x, k, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0) {

  ## make sure x is numeric
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")


  k <- as.integer(k)

  if(k < 1) stop("Illegal k: needs to be k>=1!")
  if(k >= nrow(x)) stop("Not enought neighbors in data set!")

  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR"))
  if(is.na(search)) stop("Unknown NN search type!")

  ## returns NO self matches
  ret <- kNN_int(as.matrix(x), as.integer(k),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx))

  ### sort entries (by dist and id)?
  if(sort && k>1) {
    o <- sapply(1:nrow(ret$dist), FUN =
        function(i) order(ret$dist[i,], ret$id[i,], decreasing=FALSE))
    for(i in 1:ncol(o)) {
      ret$dist[i,] <- ret$dist[i,][o[,i]]
      ret$id[i,] <- ret$id[i,][o[,i]]
      }
  }

  dimnames(ret$dist) <- list(rownames(x), 1:k)
  dimnames(ret$id) <- list(rownames(x), 1:k)

  ret
}
