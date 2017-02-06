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

sNN <- function(x, k, sort = TRUE, search = "kdtree", bucketSize = 10,
  splitRule = "suggest", approx = 0){

  if(is(x, "kNN")) {
    nn <- x
    if(missing(k)) k <- nn$k
    if(ncol(nn$id) < k) stop("kNN object does not contain enough neighbors!")
    if(ncol(nn$id) < k) {
      if(!nn$sort) nn <- sort(nn)
      nn$id <- nn$id[,1:k]
      nn$dist <- nn$dist[,1:k]
    }
  } else nn <- kNN(x, k, sort = FALSE, search = search, bucketSize = bucketSize,
    splitRule = splitRule, approx = approx)

  nn$shared <- dbscan:::SNN_sim_int(nn$id)
  class(nn) <- c("NN", "kNN", "sNN")

  if(sort) nn <- sort(nn)

  nn
}
