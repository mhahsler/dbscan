#######################################################################
# dbscan - Density Based Clustering of Applications with Noise
#          and Related Algorithms
# Copyright (C) 2017 Michael Hahsler

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

sNNclust <- function(x, k, minstr, minlinks, ...) {
  eps <- k - minstr
  minPts <- minlinks

  # find kNN
  if(is(x, "kNN")) {
    if(k > ncol(x$id))
      stop("Provided kNN object contains less than ", k, "neighbors!")
    if(k != ncol(x$id) && !x$sort)
      stop("Provided kNN object needs to be sorted! Use kNN() with sort = TRUE")

    nn <- x$id[,1:k]
  } else nn <- kNN(x, k, ..., sort = FALSE)$id

  # calculate JP clustering similarity and convert to distance
  d <- k - SNN_sim_int(nn)

  # convert into a frNN object
  remove <- d > eps
  nn_list <- lapply(seq(nrow(nn)), FUN = function(i) unname(nn[i, !remove[i,]]))
  snn <- structure(list(id = nn_list, eps = eps),
    class = c("NN", "frNN"))

  # run dbscan
  cl <- dbscan(snn, minPts = minPts, borderPoints = FALSE)
  cl$cluster
}
