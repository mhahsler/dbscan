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


jpclust <- function(x, k, kt, ...) {

  # Step 1
  if(is(x, "kNN")) {
    if(missing(k)) k <- nn$k
    nn <- x$id[,1:k]
  } else {
    nn <- kNN(x, k, sort = FALSE, ...)$id
  }

  if(length(kt) != 1 || kt < 1 || kt > k)
    stop("kt needs to be a threshold in range [1, k].")

  # Step 2
  cl <- JP_int(nn, kt = as.integer(kt))

  # Step 3
  as.integer(factor(cl))
}
