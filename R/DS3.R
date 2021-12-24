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


#' DS3: Spatial data with arbitrary shapes
#'
#' Contains 8000 2-d points, with 6 "natural" looking shapes, all of which have
#' an sinusoid-like shape that intersects with each cluster.
#'
#' Originally used as a benchmark data set for the Chameleon clustering
#' algorithm (Karypis, Han and Kumar, 1999) to
#' illustrate the a data set containing arbitrarily shaped
#' spatial data surrounded by both noise and artifacts.
#'
#' @name DS3
#' @docType data
#' @format A data frame with 8000 observations on the following 2 variables.
#' \describe{
#' \item{X}{a numeric vector}
#' \item{Y}{a numeric vector} }
#'
#' @references Karypis, George, Eui-Hong Han, and Vipin Kumar (1999).
#' "Chameleon: Hierarchical clustering using dynamic modeling." _Computer_
#' 32(8): 68-75.
#' @source Obtained from \url{https://cs.joensuu.fi/sipu/datasets/}
#' @keywords datasets
#' @examples
#'
#' data(DS3)
#' plot(DS3, pch = 20, cex = 0.25)
#'
NULL
