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

#' DBCV Paper Datasets
#'
#' The four synthetic 2D datasets used in Moulavi et al (2014).
#'
#' @name DBCV_datasets
#' @aliases Dataset_1 Dataset_2 Dataset_3 Dataset_4
#' @docType data
#' @format Four data frames with the following 3 variables.
#' \describe{
#' \item{x}{a numeric vector}
#' \item{y}{a numeric vector}
#' \item{class}{an integer vector indicating the class label. 0 means noise.} }
#' @references Davoud Moulavi and Pablo A. Jaskowiak and
#' Ricardo J. G. B. Campello and Arthur Zimek and Jörg Sander (2014).
#' Density-Based Clustering Validation. In
#' _Proceedings of the 2014 SIAM International Conference on Data Mining,_
#' pages 839-847
#' \doi{10.1137/1.9781611973440.96}
#' @source https://github.com/pajaskowiak/dbcv
#' @keywords datasets
#' @examples
#' data("Dataset_1")
#' clplot(Dataset_1[, c("x", "y")], cl = Dataset_1$class)
#'
#' data("Dataset_2")
#' clplot(Dataset_2[, c("x", "y")], cl = Dataset_2$class)
#'
#' data("Dataset_3")
#' clplot(Dataset_3[, c("x", "y")], cl = Dataset_3$class)
#'
#' data("Dataset_4")
#' clplot(Dataset_4[, c("x", "y")], cl = Dataset_4$class)
NULL



