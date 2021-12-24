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

#' Model Predictions for DBSCAN Clusterings
#'
#' Predict the membership of a new point given a clustering.
#'
#' For DBSCAN: If a new point is in the eps neighborhood of a point in the clustering, the
#' the label of that point is predicted, otherwise the prediction is `0` (noise).
#'
#' For OPTICS: Extracts a DBSCAN clustering and use predict for DBSCAN.
#'
#' For HDBSCAN: Find the k-nearest neighbors with the smallest mutual
#' reachability distance to predict the label
#'
#' @name predict
#' @aliases predict
#'
#' @param object clustering object.
#' @param data the data set used to create the clustering object.
#' @param newdata new data points for which the cluster membership should be
#' predicted.
#' @param ... further arguments.
NULL

#' @rdname predict
predict.dbscan_fast <- function (object, newdata, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

#' @rdname predict
predict.optics <- function (object, newdata, data, ...) {
  if (is.null(object$cluster))
    stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}


#' @rdname predict
predict.hdbscan <- function(object, newdata, data, ...) {
  k <- object$minPts - 1

  knn <- kNN(data, newdata, k = k , sort = TRUE)
  nn <- knn$id
  nn_dist <- knn$dist
  core_dist <- nn_dist[, k]

  nn_mrd <- mrd_m(nn_dist, core_dist)

  nn_id <- apply(nn_mrd, MARGIN = 1, which.min)
  object$cluster[nn[cbind(seq(nrow(nn)), nn_id)]]
}

.predict_frNN <- function(newdata, data, clusters, eps, ...) {
  if (is.null(newdata))
    return(clusters)

  # calculate the frNN between newdata and data (only keep entries for newdata)
  nn <- frNN(rbind(data, newdata), eps = eps,
    sort = TRUE, ...)$id[-(1:nrow(data))]

  # get rid of indices of newdata and pick first neighbor
  cl <- sapply(nn, function(x)
    clusters[x[x <= nrow(data) & x > 0][1]])
  cl[is.na(cl)] <- 0L
  cl
}
