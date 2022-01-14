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

#' @rdname dbscan
#' @param object clustering object.
#' @param data the data set used to create the clustering object.
#' @param newdata new data points for which the cluster membership should be
#' predicted.
predict.dbscan_fast <- function (object, newdata, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

#' @rdname optics
#' @param object clustering object.
#' @param data the data set used to create the clustering object.
#' @param newdata new data points for which the cluster membership should be
#' predicted.
predict.optics <- function (object, newdata, data, ...) {
  if (is.null(object$cluster) ||
      is.null(object$eps_cl) || is.na(object$eps_cl))
    stop("no extracted clustering available in object! run extractDBSCAN() first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}

#' @rdname hdbscan
#' @param object clustering object.
#' @param data the data set used to create the clustering object.
#' @param newdata new data points for which the cluster membership should be
#' predicted.
predict.hdbscan <- function(object, newdata, data, ...) {
  clusters <- object$cluster

  if (is.null(newdata))
    return(clusters)

  # don't use noise
  coredist <- object$coredist[clusters != 0]
  data <- data[clusters != 0,]
  clusters <- clusters[clusters != 0]

  # find minPts - 1 nearest neighbor
  nns <- kNN(data, query = newdata, k = 1)

  # choose cluster if dist <= coredist of that point
  drop(ifelse(nns$dist > coredist[nns$id], 0L, clusters[nns$id]))
}

## find the cluster id of the closest NN in the eps neighborhood or return 0 otherwise.
.predict_frNN <- function(newdata, data, clusters, eps, ...) {
  if (is.null(newdata))
    return(clusters)

  # don't use noise
  data <- data[clusters != 0,]
  clusters <- clusters[clusters != 0]

  # calculate the frNN between newdata and data (only keep entries for newdata)
  nn <- frNN(data,
    query = newdata,
    eps = eps,
    sort = TRUE,
    ...)

  sapply(nn$id, function(nns) {
    if (length(nns) == 0)
      0L
    else
      clusters[nns[1]]
  })
}
