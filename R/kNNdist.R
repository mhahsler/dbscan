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

kNNdist <- function(x, k, ...) dbscan::kNN(x, k, sort = TRUE, ...)$dist

kNNdistplot <- function(x, k = 4, ...) {
  kNNdist <- sort(kNNdist(x, k ,...))
  plot(sort(kNNdist), type="l", ylab=paste(k, "-NN distance", sep=""),
    xlab = "Pointes (sample) sorted by distance")
}
