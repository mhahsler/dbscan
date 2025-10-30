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

.matrixlike <- function(x) {
  if  (is.null(dim(x)))
       return(FALSE)

  # check that there is at least one row and one column!
  if (nrow(x) < 1L) stop("the provided data has 0 rows!")
  if (ncol(x) < 1L) stop("the provided data has 0 columns!")

  TRUE
}
