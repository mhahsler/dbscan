.predict_frNN <- function(newdata, data, clusters, eps, ...) {
  if (is.null(newdata)) return(clusters)

  # calculate the frNN between newdata and data (only keep entries for newdata)
  nn <- frNN(rbind(data, newdata), eps = eps,
    sort = TRUE, ...)$id[-(1:nrow(data))]

  # get rid of indices of newdata and pick first neighbor
  cl <- sapply(nn, function(x) clusters[x[x<=nrow(data) & x>0][1]])
  cl[is.na(cl)] <- 0L
  cl
}

predict.dbscan_fast <- function (object, newdata, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

predict.optics <- function (object, newdata, data, ...) {
  if (is.null(object$cluster)) stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}

# Find the k-nearest Neighbor with the smallest mutual reachability distance to predict
# the label
predict.hdbscan <- function(object, newdata, data, ...) {
  k <- object$minPts - 1

  knn <- kNN(data, newdata, k = k , sort = TRUE)
  nn <- knn$id
  nn_dist <- knn$dist
  core_dist <- nn_dist[, k]

  nn_mrd <- mrd_m(nn_dist, core_dist)

  nn_id <- apply(nn_mrd, MARGIN = 1, which.min)
  object$cluster[nn[cbind(seq(nrow(nn)),nn_id)]]
}
