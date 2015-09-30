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

frNN <- function(x, eps, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0) {

  ## make sure x is numeric
  x <- as.matrix(x)
  if(storage.mode(x) == "integer") storage.mode(x) <- "double"
  if(storage.mode(x) != "double") stop("x has to be a numeric matrix.")

  splitRule <- pmatch(toupper(splitRule), .ANNsplitRule)-1L
  if(is.na(splitRule)) stop("Unknown splitRule!")

  search <- pmatch(toupper(search), c("KDTREE", "LINEAR"))
  if(is.na(search)) stop("Unknown NN search type!")

  ret <- frNN_int(as.matrix(x), as.double(eps),
    as.integer(search), as.integer(bucketSize),
    as.integer(splitRule), as.double(approx))

  if(sort) {
    o <- lapply(1:length(ret$dist), FUN =
        function(i) order(ret$dist[[i]], ret$id[[i]], decreasing=FALSE))
    ret$dist <- lapply(1:length(o), FUN = function(p) ret$dist[[p]][o[[p]]])
    ret$id <- lapply(1:length(o), FUN = function(p) ret$id[[p]][o[[p]]])
  }

  names(ret$dist) <- rownames(x)
  names(ret$id) <- rownames(x)

  ret$eps <- eps

  ret
}
