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

#' Moons Data
#'
#' Contains 100 2-d points, half of which are contained in two moons or
#' "blobs"" (25 points each blob), and the other half in asymmetric facing
#' crescent shapes. The three shapes are all linearly separable.
#'
#' This data was generated with the following Python commands using the
#' SciKit-Learn library:
#'
#' `> import sklearn.datasets as data`
#'
#' `> moons = data.make_moons(n_samples=50, noise=0.05)`
#'
#' `> blobs = data.make_blobs(n_samples=50, centers=[(-0.75,2.25), (1.0, 2.0)], cluster_std=0.25)`
#'
#' `> test_data = np.vstack([moons, blobs])`
#'
#' @name moons
#' @docType data
#' @format A data frame with 100 observations on the following 2 variables.
#' \describe{
#' \item{X}{a numeric vector}
#' \item{Y}{a numeric vector} }
#' @references Pedregosa, Fabian, Gael Varoquaux, Alexandre Gramfort,
#' Vincent Michel, Bertrand Thirion, Olivier Grisel, Mathieu Blondel et al.
#' Scikit-learn: Machine learning in Python. _Journal of Machine Learning
#' Research_ 12, no. Oct (2011): 2825-2830.
#' @source See the HDBSCAN notebook from github documentation:
#' \url{http://hdbscan.readthedocs.io/en/latest/how_hdbscan_works.html}
#' @keywords datasets
#' @examples
#' data(moons)
#' plot(moons, pch=20)
NULL



