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

lof <- function(x, k = 4, ...) {

  # get n
  if(is(x, "dist")) n <- attr(x, "Size")
  else n <- nrow(x)

  if(is.null(n)) stop("x needs to be a matrix or a dist object!")

  if(k<1 || k>=n)
    stop("k has to be larger than 1 and smaller than the number of points")

  # get k nearest neighbors + distances
  d <- kNN(x, k, ...)

  # calculate local reachability density
  lrd <- numeric(n)
  for(i in 1:n) lrd[i] <- 1/(sum(apply(
    cbind(d$dist[d$id[i,], k], d$dist[i,]),
    1, max)) / k)

  # calculate lof
  lof <- numeric(n)
  for (i in 1:n) lof[i] <- sum(lrd[d$id[i,]])/k / lrd[i]

  # with more than MinPts duplicates lrd can become infinity
  # we define them not to be outliers
  lof[is.nan(lof)] <- 1

  lof
}
