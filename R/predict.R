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

predict.dbscan_fast <- function (object, newdata = NULL, data, ...)
  .predict_frNN(newdata, data, object$cluster, object$eps, ...)

predict.optics <- function (object, newdata = NULL, data, ...) {
  if (is.null(object$cluster)) stop("no extracted clustering available in object! run extractDBSCAN or extractXi first.")
  .predict_frNN(newdata, data, object$cluster, object$eps_cl, ...)
}
