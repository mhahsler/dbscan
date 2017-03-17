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

sNNclust <- function(x, k, eps, minPts, borderPoints = TRUE, ...) {
  nn <- sNN(x, k=k, ...)

  # convert into a frNN object which already enforces eps
  nn_list <- lapply(seq(nrow(nn$id)),
    FUN = function(i) unname(nn$id[i, nn$shared[i,] >= eps]))
  snn <- structure(list(id = nn_list, eps = eps),
    class = c("NN", "frNN"))

  # run dbscan
  cl <- dbscan(snn, minPts = minPts, borderPoints = borderPoints)

  structure(list(cluster = cl$cluster,
    type = "SharedNN clustering",
    param = list(k = k, eps = eps, minPts = minPts, borderPoints = borderPoints)),
    class = c("general_clustering"))
}
